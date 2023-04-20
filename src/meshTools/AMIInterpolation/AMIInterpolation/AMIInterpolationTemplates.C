/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "profiling.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
namespace AMI
{
template<class T>
class MultiplyWeightOp
{
public:

    MultiplyWeightOp()
    {}

    void operator()
    (
        T& x,
        const UList<T>& values,
        const UList<label>& faces,
        const UList<scalar>& weights,
        const scalar /* unused sumWeight */
    ) const
    {
        x = Zero;
        forAll(faces, i)
        {
            const label facei = faces[i];
            const scalar w = weights[i];

            x += w*values[facei];
        }
    }
};


template<class T>
class StabiliseWeightOp
:
    public MultiplyWeightOp<T>
{
public:

    StabiliseWeightOp()
    {}

    void operator()
    (
        T& x,
        const UList<T>& values,
        const UList<label>& faces,
        const UList<scalar>& weights,
        const scalar sumWeight
    ) const
    {
        const T corr = (1 - sumWeight)*x;

        MultiplyWeightOp<T>::operator()(x, values, faces, weights, sumWeight);

        x += corr;
    }
};


template<class T>
class NormaliseWeightOp
:
    public MultiplyWeightOp<T>
{
public:

    NormaliseWeightOp()
    {}

    void operator()
    (
        T& x,
        const UList<T>& values,
        const UList<label>& faces,
        const UList<scalar>& weights,
        const scalar sumWeight
    ) const
    {
        MultiplyWeightOp<T>::operator()(x, values, faces, weights, sumWeight);

        x /= sumWeight;
    }
};


template<class T>
class AverageOp
:
    public MultiplyWeightOp<T>
{
public:

    AverageOp()
    {}

    void operator()
    (
        T& x,
        const UList<T>& values,
        const UList<label>& faces,
        const UList<scalar>& weights,
        const scalar sumWeight
    ) const
    {
        const T xOld = x;

        MultiplyWeightOp<T>::operator()(x, values, faces, weights, sumWeight);

        x = 0.5*(x + xOld);
    }
};


template<class T, class Cop>
class CombineOp
{
    // Combine operator, e.g.minOp, maxOp
    Cop cop_;


public:

    CombineOp(const Cop& cop)
    :
        cop_(cop)
    {}

    void operator()
    (
        T& x,
        const UList<T>& values,
        const UList<label>& faces,
        const UList<scalar>&, /* unused weights */
        const scalar /* unused sumWeight */
    ) const
    {
        x = values[faces[0]];
        for (label facei=1; facei<faces.size(); ++facei)
        {
            cop_(x, values[facei]);
        }
    }
};

} // End namespace AMI
} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToSource3
(
    const UList<Type>& tgtFld,
    List<Type>& srcFld,
    const CombineOp& cop
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToSource3");

    if (tgtFld.size() != tgtAddress_.size())
    {
        FatalErrorInFunction
            << "Target field size is not equal to target patch size" << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << tgtFld.size()
            << abort(FatalError);
    }

    if (srcFld.size() != srcAddress_.size())
    {
        FatalErrorInFunction
            << "Source field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    supplied field = " << srcFld.size()
            << abort(FatalError);
    }

    if (distributed())
    {
        List<Type> work(tgtFld);
        tgtMapPtr_->distribute(work);

        forAll(srcFld, facei)
        {
            if (srcWeightsSum_[facei] > lowWeightCorrection_)
            {
                const labelList& faces = srcAddress_[facei];
                const scalarList& weights = srcWeights_[facei];

                cop(srcFld[facei], work, faces, weights, srcWeightsSum_[facei]);
            }
        }
    }
    else
    {
        forAll(srcFld, facei)
        {
            if (srcWeightsSum_[facei] > lowWeightCorrection_)
            {
                const labelList& faces = srcAddress_[facei];
                const scalarList& weights = srcWeights_[facei];

                cop
                (
                    srcFld[facei],
                    tgtFld,
                    faces,
                    weights,
                    srcWeightsSum_[facei]
                );
            }
        }
    }
}


template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToTarget3
(
    const UList<Type>& srcFld,
    List<Type>& tgtFld,
    const CombineOp& cop
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToTarget3");

    if (srcFld.size() != srcAddress_.size())
    {
        FatalErrorInFunction
            << "Source field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    supplied field = " << srcFld.size()
            << abort(FatalError);
    }

    if (tgtFld.size() != tgtAddress_.size())
    {
        FatalErrorInFunction
            << "Target field size is not equal to target patch size" << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << tgtFld.size()
            << abort(FatalError);
    }

    if (distributed())
    {
        List<Type> work(srcFld);
        srcMapPtr_->distribute(work);

        forAll(tgtFld, facei)
        {
            if (tgtWeightsSum_[facei] > lowWeightCorrection_)
            {
                const labelList& faces = tgtAddress_[facei];
                const scalarList& weights = tgtWeights_[facei];

                cop(tgtFld[facei], work, faces, weights, tgtWeightsSum_[facei]);
            }
        }
    }
    else
    {
        forAll(tgtFld, facei)
        {
            if (tgtWeightsSum_[facei] > lowWeightCorrection_)
            {
                const labelList& faces = tgtAddress_[facei];
                const scalarList& weights = tgtWeights_[facei];

                cop
                (
                    tgtFld[facei],
                    srcFld,
                    faces,
                    weights,
                    tgtWeightsSum_[facei]
                );
            }
        }
    }
}


template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToTarget");

    if (fld.size() != srcAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != tgtAddress_.size())
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    target patch   = " << tgtAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(tgtAddress_.size());

    if (distributed())
    {
        const mapDistribute& map = srcMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, facei)
        {
            if (tgtWeightsSum_[facei] < lowWeightCorrection_)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& faces = tgtAddress_[facei];
                const scalarList& weights = tgtWeights_[facei];

                forAll(faces, i)
                {
                    cop(result[facei], facei, work[faces[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            if (tgtWeightsSum_[facei] < lowWeightCorrection_)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& faces = tgtAddress_[facei];
                const scalarList& weights = tgtWeights_[facei];

                forAll(faces, i)
                {
                    cop(result[facei], facei, fld[faces[i]], weights[i]);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
void Foam::AMIInterpolation::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    addProfiling(ami, "AMIInterpolation::interpolateToSource");

    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != srcAddress_.size())
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to source "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    source patch   = " << srcAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(srcAddress_.size());

    if (distributed())
    {
        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, facei)
        {
            if (srcWeightsSum_[facei] < lowWeightCorrection_)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& faces = srcAddress_[facei];
                const scalarList& weights = srcWeights_[facei];

                forAll(faces, i)
                {
                    cop(result[facei], facei, work[faces[i]], weights[i]);
                }
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            if (srcWeightsSum_[facei] < lowWeightCorrection_)
            {
                result[facei] = defaultValues[facei];
            }
            else
            {
                const labelList& faces = srcAddress_[facei];
                const scalarList& weights = srcWeights_[facei];

                forAll(faces, i)
                {
                    cop(result[facei], facei, fld[faces[i]], weights[i]);
                }
            }
        }
    }
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(srcAddress_.size(), Zero);

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), cop, defaultValues);
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    auto tresult = tmp<Field<Type>>::New(tgtAddress_.size(), Zero);

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), cop, defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>(), defaultValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::AMIInterpolation::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>(), defaultValues);
}


// ************************************************************************* //
