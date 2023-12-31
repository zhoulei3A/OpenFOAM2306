/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

Class
    Foam::GeometricField

Description
    Generic GeometricField class.

SourceFiles
    GeometricFieldI.H
    GeometricField.C
    GeometricFieldFunctions.H
    GeometricFieldFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GeometricField_H
#define Foam_GeometricField_H

#include "regIOobject.H"
#include "GeometricBoundaryField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);


/*---------------------------------------------------------------------------*\
                           Class GeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField
:
    public DimensionedField<Type, GeoMesh>
{
public:

    // Public Typedefs

        //- The mesh type for the GeometricField
        typedef typename GeoMesh::Mesh Mesh;

        //- The boundary mesh type for the boundary fields
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        //- The internal field type from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh> Internal;

        //- Type of boundary fields
        typedef GeometricBoundaryField<Type, PatchField, GeoMesh> Boundary;

        //- The patch field type for the GeometricBoundaryField
        typedef PatchField<Type> Patch;

        //- Component type of the field elements
        typedef typename Field<Type>::cmptType cmptType;


private:

    // Private Data

        //- Current time index.
        //  Used to trigger the storing of the old-time value
        mutable label timeIndex_;

        //- Pointer to old time field
        mutable GeometricField<Type, PatchField, GeoMesh>* field0Ptr_;

        //- Pointer to previous iteration (used for under-relaxation)
        mutable GeometricField<Type, PatchField, GeoMesh>* fieldPrevIterPtr_;

        //- Boundary field containing boundary field values
        Boundary boundaryField_;


    // Private Member Functions

        //- Read from file if it is present
        bool readIfPresent();

        //- Read old time field from file if it is present
        bool readOldTimeIfPresent();

        //- Read the field from the dictionary
        void readFields(const dictionary& dict);

        //- Read the field - create the field dictionary on-the-fly
        void readFields();

public:

    //- Runtime type information
    TypeName("GeometricField");


    // Static Member Functions

        //- Return a null geometric field
        inline static const GeometricField<Type, PatchField, GeoMesh>& null();


    // Constructors

        //- Construct given IOobject, mesh, dimensions and patch type.
        //  This allocates storage for the field but does not set values.
        //  Used only within this class to create TEMPORARY variables
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Construct given IOobject, mesh, dimensions and patch types.
        //  This allocates storage for the field but does not set values.
        //  Used only within this class to create TEMPORARY variables
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct given IOobject, mesh, dimensions, initial field value
        //- and given patch type.
        //  This assigns both dimensions and values.
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const Type& value,
            const dimensionSet& dims,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Construct given IOobject, mesh, dimensions, initial field value
        //- and given patch types.
        //  This assigns both dimensions and values.
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const Type& value,
            const dimensionSet& dims,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct given IOobject, mesh, dimensioned<Type> and patch type.
        //  This assigns both dimensions and values.
        //  The name of the dimensioned\<Type\> has no influence.
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Construct given IOobject, mesh, dimensioned<Type> and patch types.
        //  This assigns both dimensions and values.
        //  The name of the dimensioned\<Type\> has no influence.
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Copy construct from internal field and a patch list to clone
        GeometricField
        (
            const IOobject& io,
            const Internal& diField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Move construct from internal field and a patch list to clone
        GeometricField
        (
            const IOobject& io,
            Internal&& diField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Move construct from internal field and a patch list to clone
        GeometricField
        (
            const IOobject& io,
            const tmp<Internal>& tfield,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Copy construct from internal field and a patch list to clone
        GeometricField
        (
            const Internal& diField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Move construct from internal field and a patch list to clone
        GeometricField
        (
            Internal&& diField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Copy construct from internal field, with specified patch type
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const Field<Type>& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Move construct from internal field, with specified patch type
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            Field<Type>&& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Copy construct from components
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const Field<Type>& iField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Move construct from internal field and a patch list to clone
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            Field<Type>&& iField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Copy construct from components
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& dims,
            const tmp<Field<Type>>& tiField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Construct and read given IOobject
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const bool readOldTime = true
        );

        //- Construct from dictionary
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dictionary& dict
        );

        //- Copy construct
        GeometricField
        (
            const GeometricField<Type, PatchField, GeoMesh>& gf
        );

        //- Construct from tmp\<GeometricField\> deleting argument
        GeometricField
        (
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
        );

        //- Construct as copy resetting IO parameters
        GeometricField
        (
            const IOobject& io,
            const GeometricField<Type, PatchField, GeoMesh>& gf
        );

        //- Construct from tmp\<GeometricField\> resetting IO parameters
        GeometricField
        (
            const IOobject& io,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
        );

        //- Copy construct with a new name
        GeometricField
        (
            const word& newName,
            const GeometricField<Type, PatchField, GeoMesh>& gf
        );

        //- Construct with a new name from tmp\<GeometricField\>
        GeometricField
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
        );

        //- Construct as copy resetting IO parameters and patch type
        GeometricField
        (
            const IOobject& io,
            const GeometricField<Type, PatchField, GeoMesh>& gf,
            const word& patchFieldType
        );

        //- Construct as copy resetting IO parameters and boundary type
        //- for selected patchIDs
        GeometricField
        (
            const IOobject& io,
            const GeometricField<Type, PatchField, GeoMesh>& gf,
            const labelList& patchIDs,
            const word& patchFieldType
        );

        //- Construct as copy resetting IO parameters and boundary types
        GeometricField
        (
            const IOobject& io,
            const GeometricField<Type, PatchField, GeoMesh>& gf,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct as copy resetting IO parameters and boundary types
        GeometricField
        (
            const IOobject& io,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Clone
        tmp<GeometricField<Type, PatchField, GeoMesh>> clone() const;


    // Factory Methods

        //- Return tmp field from name, mesh, dimensions and patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& dims,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensions,
        //- copy of internal field, with specified patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& dims,
            const Field<Type>& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensions,
        //- moved internal field contents, with specified patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& dims,
            Field<Type>&& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensions, initial field value
        //- and patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const Type& value,
            const dimensionSet& dims,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensioned\<Type\>
        //- and patch types.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const Type& value,
            const dimensionSet& dims,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Return tmp field from name, mesh, dimensioned\<Type\>
        //- and patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensioned\<Type\>
        //- and patch types.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Return renamed tmp field
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
        );

        //- Rename tmp field and reset patch field type
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
            const word& patchFieldType
        );

        //- Rename tmp field and reset patch field types and return
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct tmp field based on mesh/registry information from
        //- an existing field.
        //  Created NO_READ, NO_WRITE, NO_REGISTER, using the instance
        //  from the field
        template<class AnyType>
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const GeometricField<AnyType, PatchField, GeoMesh>& fld,
            const word& name,
            const dimensionSet& dims,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Construct tmp field based on mesh/registry information from
        //- an existing field and initialise with value.
        //  Created NO_READ, NO_WRITE, NO_REGISTER, using the instance
        //  from the field
        template<class AnyType>
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const GeometricField<AnyType, PatchField, GeoMesh>& fld,
            const word& name,
            const dimensioned<Type>& dt,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );


    //- Destructor
    virtual ~GeometricField();


    // Member Functions

        //- Return a const-reference to the dimensioned internal field.
        inline const Internal& internalField() const noexcept;

        //- Return a reference to the dimensioned internal field.
        //  \param updateAccessTime update event counter and check
        //      old-time fields
        //
        //  \note Should avoid using updateAccessTime = true within loops.
        Internal& internalFieldRef(const bool updateAccessTime = true);

        //- Same as internalFieldRef()
        Internal& ref(const bool updateAccessTime = true)
        {
            return this->internalFieldRef(updateAccessTime);
        }

        //- Return a const-reference to the dimensioned internal field
        //- of a "vol" field.
        //  Useful in the formulation of source-terms for FV equations
        //
        //  \note definition in finiteVolume/fields/volFields/volFieldsI.H
        inline const Internal& v() const;

        //- Return a const-reference to the internal field values.
        inline const typename Internal::FieldType& primitiveField()
        const noexcept;

        //- Return a reference to the internal field values.
        //  \param updateAccessTime update event counter and check
        //      old-time fields
        //
        //  \note Should avoid using updateAccessTime = true within loops.
        typename Internal::FieldType& primitiveFieldRef
        (
            const bool updateAccessTime = true
        );

        //- Return const-reference to the boundary field
        inline const Boundary& boundaryField() const noexcept;

        //- Return a reference to the boundary field
        //  \param updateAccessTime update event counter and check
        //      old-time fields
        //
        //  \note Should avoid using updateAccessTime = true within loops.
        Boundary& boundaryFieldRef(const bool updateAccessTime = true);

        //- Return the time index of the field
        inline label timeIndex() const noexcept;

        //- Write-access to the time index of the field
        inline label& timeIndex() noexcept;

        //- Store the old-time fields
        void storeOldTimes() const;

        //- Store the old-time field
        void storeOldTime() const;

        //- Return the number of old time fields stored
        label nOldTimes() const;

        //- Return old time field
        const GeometricField<Type, PatchField, GeoMesh>& oldTime() const;

        //- Return non-const old time field
        //  (Not a good idea but it is used for sub-cycling)
        GeometricField<Type, PatchField, GeoMesh>& oldTime();

        //- Store the field as the previous iteration value
        void storePrevIter() const;

        //- Return previous iteration field
        const GeometricField<Type, PatchField, GeoMesh>& prevIter() const;

        //- Correct boundary field
        void correctBoundaryConditions();

        //- Does the field need a reference level for solution
        bool needReference() const;

        //- Return a component of the field
        tmp<GeometricField<cmptType, PatchField, GeoMesh>> component
        (
            const direction
        ) const;

        //- WriteData member function required by regIOobject
        bool writeData(Ostream&) const;

        //- Return transpose (only if it is a tensor field)
        tmp<GeometricField<Type, PatchField, GeoMesh>> T() const;

        //- Relax field (for steady-state solution).
        //  alpha = 1 : no relaxation
        //  alpha < 1 : relaxation
        //  alpha = 0 : do nothing
        void relax(const scalar alpha);

        //- Relax field (for steady-state solution).
        //  alpha is read from controlDict
        void relax();

        //- Select the final iteration parameters if \c final is true
        //- by returning the field name + "Final"
        //- otherwise the standard parameters by returning the field name
        word select(bool final) const;

        //- Helper function to write the min and max to an Ostream
        void writeMinMax(Ostream& os) const;


    // Member Function *this Operators

        //- Negate the field inplace. See notes in Field
        void negate();

        //- Normalise the field inplace. See notes in Field
        void normalise();

        //- Replace specified field component with content from another field
        void replace
        (
            const direction d,
            const GeometricField<cmptType, PatchField, GeoMesh>& gcf
        );

        //- Replace specified field component with specified value
        void replace
        (
            const direction d,
            const dimensioned<cmptType>& ds
        );

        //- Impose lower (floor) clamp on the field values (in-place)
        void clamp_min(const Type& lower);

        //- Impose upper (ceiling) clamp on the field values (in-place)
        void clamp_max(const Type& upper);

        //- Impose lower (floor) clamp on the field values (in-place)
        //  No dimension checking
        void clamp_min(const dimensioned<Type>& lower);

        //- Impose upper (ceiling) clamp on the field values (in-place)
        //  No dimension checking
        void clamp_max(const dimensioned<Type>& upper);

        //- Clamp field values (in-place) to the specified range.
        //  Does not check if range is valid or not. No dimension checking.
        void clamp_range(const dimensioned<MinMax<Type>>& range);

        //- Clamp field values (in-place) to the specified range.
        //  Does not check if range is valid or not.
        void clamp_range(const Type& lower, const Type& upper);

        //- Clamp field values (in-place) to the specified range.
        //  Does not check if range is valid or not.
        void clamp_range(const MinMax<Type>& range);

        //- Clamp field values (in-place) to the specified range.
        //  Does not check if range is valid or not. No dimension checking.
        void clamp_range
        (
            const dimensioned<Type>& lower,
            const dimensioned<Type>& upper
        );


    // Member Operators

        //- Return a const-reference to the dimensioned internal field.
        //- Same as internalField().
        //  Useful in the formulation of source-terms for FV equations
        const Internal& operator()() const { return *this; }

        void operator=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);
        void operator=(const dimensioned<Type>&);

        void operator==(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);
        void operator==(const dimensioned<Type>&);

        void operator+=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator+=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        void operator-=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator-=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        void operator*=(const GeometricField<scalar, PatchField, GeoMesh>&);
        void operator*=(const tmp<GeometricField<scalar,PatchField,GeoMesh>>&);

        void operator/=(const GeometricField<scalar, PatchField, GeoMesh>&);
        void operator/=(const tmp<GeometricField<scalar,PatchField,GeoMesh>>&);

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator*=(const dimensioned<scalar>&);
        void operator/=(const dimensioned<scalar>&);


    // Ostream Operators

        friend Ostream& operator<< <Type, PatchField, GeoMesh>
        (
            Ostream&,
            const GeometricField<Type, PatchField, GeoMesh>&
        );

        friend Ostream& operator<< <Type, PatchField, GeoMesh>
        (
            Ostream&,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&
        );


    // Housekeeping

        //- Clamp field values (in-place) to the specified range.
        //  \deprecated(2023-01) prefer clamp_range() naming
        void clip(const dimensioned<MinMax<Type>>& range)
        {
            this->clamp_range(range);
        }

        //- Clamp field values (in-place) to the specified range.
        //  \deprecated(2023-01) prefer clamp_range() naming
        void clip(const dimensioned<Type>& lo, const dimensioned<Type>& hi)
        {
            this->clamp_range(lo.value(), hi.value());
        }

        //- Use minimum of the field and specified value. ie, clamp_max().
        //  This sets the \em ceiling on the field values
        //  \deprecated(2023-01) prefer clamp_max()
        void min(const dimensioned<Type>& upper) { this->clamp_max(upper); }

        //- Use maximum of the field and specified value. ie, clamp_min().
        //  This sets the \em floor on the field values
        //  \deprecated(2023-01) prefer clamp_min()
        void max(const dimensioned<Type>& lower) { this->clamp_min(lower); }

        //- Deprecated(2019-01) identical to clamp_range()
        //  \deprecated(2019-01) identical to clamp_range()
        FOAM_DEPRECATED_FOR(2019-01, "clamp_range() method")
        void maxMin(const dimensioned<Type>& lo, const dimensioned<Type>& hi)
        {
            return this->clamp_range(lo.value(), hi.value());
        }
};


template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const typename GeometricField<Type, PatchField, GeoMesh>::
    Boundary&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GeometricFieldI.H"

#ifdef NoRepository
    #include "GeometricField.C"
#endif

#include "GeometricFieldFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
