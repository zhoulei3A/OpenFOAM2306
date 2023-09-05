/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "masterOFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::masterOFstream::checkWrite
(
    const fileName& fName,
    const UList<char>& charData
)
{
    if (charData.empty())
    {
        // Can probably skip all of this if there is nothing to write
        return;
    }

    Foam::mkDir(fName.path());

    OFstream os
    (
        atomic_,
        fName,
        IOstreamOption(IOstreamOption::BINARY, version(), compression_),
        append_
    );
    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Could not open file " << fName << nl
            << exit(FatalIOError);
    }

    // Use writeRaw() to write characters directly
    os.writeRaw(charData.cdata(), charData.size_bytes());

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Failed writing to " << fName << nl
            << exit(FatalIOError);
    }
}


void Foam::masterOFstream::commit()
{
    // Retrieve serialized content
    List<char> charData;
    OCharStream::swap(charData);

    if (UPstream::parRun())
    {
        // Discard content if not writing.
        // Ensures nothing will be communicated
        if (!writeOnProc_)
        {
            charData.clear();
        }

        List<fileName> filePaths(UPstream::nProcs(comm_));
        filePaths[UPstream::myProcNo(comm_)] = pathName_;
        Pstream::gatherList(filePaths, UPstream::msgType(), comm_);

        bool uniform =
        (
            UPstream::master(comm_)
         && fileOperation::uniformFile(filePaths)
        );

        Pstream::broadcast(uniform, comm_);

        if (uniform)
        {
            if (UPstream::master(comm_) && writeOnProc_)
            {
                checkWrite(pathName_, charData);
            }
            return;
        }

        // Different files.
        // Would be slightly easier if we didn't support writeOnProc
        // but instead of using listGatherValues<bool>() to manage
        // things on the master, we simply check based on the received
        // size and only write non-zero file content.
        // This avoids an additional reduction.

        const label startOfRequests = UPstream::nRequests();

        if (UPstream::is_subrank(comm_))
        {
            // Send to master. When (!writeOnProc_) it is zero-sized.
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                UPstream::masterNo(),
                charData.cdata(),
                charData.size_bytes(),
                comm_
            );
        }
        else if (UPstream::master(comm_))
        {
            List<List<char>> procData(UPstream::nProcs(comm_));

            const UPstream::rangeType recvProcs = UPstream::subProcs(comm_);

            for (const int proci : recvProcs)
            {
                List<char>& procSlice = procData[proci];

                // Probe the incoming message source/size
                std::pair<int, int> probed =
                    UPstream::probeMessage
                    (
                        UPstream::commsTypes::blocking,
                        proci,
                        UPstream::msgType(),
                        comm_
                    );

                procSlice.resize_nocopy(probed.second);

                // Receive content (can also be zero-sized)
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    proci,
                    procSlice.data(),
                    procSlice.size_bytes(),
                    UPstream::msgType(),
                    comm_
                );
            }

            if (writeOnProc_)
            {
                // Write master data
                checkWrite(filePaths[UPstream::masterNo()], charData);
            }

            // Poll for completed receive requests and dispatch
            DynamicList<int> indices(recvProcs.size());
            while
            (
                UPstream::waitSomeRequests
                (
                    startOfRequests,
                    recvProcs.size(),
                   &indices
                )
            )
            {
                for (const int idx : indices)
                {
                    const int proci = recvProcs[idx];

                    if (!procData[proci].empty())
                    {
                        // Write non-empty sub-proc data
                        checkWrite(filePaths[proci], procData[proci]);
                    }

                    // Eager cleanup?
                    // procData[proci].clear();
                }
            }
        }

        UPstream::waitRequests(startOfRequests);
    }
    else
    {
        checkWrite(pathName_, charData);
    }

    // This method is only called once (internally)
    // so no need to clear/flush old buffered data
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterOFstream::masterOFstream
(
    IOstreamOption::atomicType atomic,
    const label comm,
    const fileName& pathName,
    IOstreamOption streamOpt,
    IOstreamOption::appendType append,
    const bool writeOnProc
)
:
    OCharStream(streamOpt),
    pathName_(pathName),
    atomic_(atomic),
    compression_(streamOpt.compression()),
    append_(append),
    writeOnProc_(writeOnProc),
    comm_(comm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterOFstream::~masterOFstream()
{
    commit();
}


// ************************************************************************* //
