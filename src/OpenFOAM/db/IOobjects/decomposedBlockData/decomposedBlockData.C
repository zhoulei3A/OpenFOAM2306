/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "decomposedBlockData.H"
#include "OPstream.H"
#include "IPstream.H"
#include "Fstream.H"
#include "SpanStream.H"
#include "dictionary.H"
#include "objectRegistry.H"
#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(decomposedBlockData, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::decomposedBlockData::isCollatedType
(
    const word& objectType
)
{
    return
    (
        objectType == decomposedBlockData::typeName
    );
}


bool Foam::decomposedBlockData::isCollatedType
(
    const IOobject& io
)
{
    return decomposedBlockData::isCollatedType(io.headerClassName());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decomposedBlockData::decomposedBlockData
(
    const label comm,
    const IOobject& io,
    const UPstream::commsTypes commsType
)
:
    regIOobject(io),
    commsType_(commsType),
    comm_(comm),
    contentData_()
{
    // Temporary warning
    if (readOpt() == IOobjectOption::READ_MODIFIED)
    {
        WarningInFunction
            << "decomposedBlockData " << name()
            << " constructed with READ_MODIFIED"
            " but decomposedBlockData does not support automatic rereading."
            << endl;
    }
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        read();
    }
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

bool Foam::decomposedBlockData::readBlockEntry
(
    Istream& is,
    List<char>& charData
)
{
    // Handle any of these:

    // 0.  NCHARS (...)
    // 1.  List<char> NCHARS (...)
    // 2.  processorN  List<char> NCHARS (...) ;
    // 3.  processorN  NCHARS (...) ;

    is.fatalCheck(FUNCTION_NAME);
    token tok(is);
    is.fatalCheck(FUNCTION_NAME);

    // Dictionary format has primitiveEntry keyword:
    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if (!isDictFormat && tok.good())
    {
        is.putBack(tok);
    }
    charData.readList(is);

    if (isDictFormat)
    {
        is.fatalCheck(FUNCTION_NAME);
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return true;
}


bool Foam::decomposedBlockData::skipBlockEntry(Istream& is)
{
    // As per readBlockEntry but seeks instead of reading.
    // Internals like charList::readList

    // Handle any of these:
    // 0.  NCHARS (...)
    // 1.  List<char> NCHARS (...)
    // 2.  processorN  List<char> NCHARS (...) ;
    // 3.  processorN  NCHARS (...) ;

    if (!is.good()) return false;
    token tok(is);
    if (!is.good()) return false;

    // Dictionary format has primitiveEntry keyword:
    const bool isDictFormat = (tok.isWord() && !tok.isCompound());

    if (isDictFormat)
    {
        is >> tok;
        if (!is.good()) return false;
    }


    bool handled = false;

    // Like charList::readList
    if (tok.isCompound())
    {
        handled = true;
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..) or just a plain '0'

        const label len = tok.labelToken();

        // Binary, always contiguous

        if (len)
        {
            const auto oldFmt = is.format(IOstreamOption::BINARY);

            // read(...) includes surrounding start/end delimiters.

            // Note: nullptr to ignore instead of reading
            is.read(nullptr, std::streamsize(len));

            is.format(oldFmt);
        }

        handled = true;
    }
    else
    {
        // Incorrect token
        return false;
    }

    if (isDictFormat)
    {
        is.fatalCheck(FUNCTION_NAME);
        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        // Swallow trailing ';'
        if (tok.good() && !tok.isPunctuation(token::END_STATEMENT))
        {
            is.putBack(tok);
        }
    }

    return handled;
}


Foam::label Foam::decomposedBlockData::getNumBlocks
(
    Istream& is,
    const label maxNumBlocks
)
{
    label nBlocks = 0;

    // Handle OpenFOAM header if it is the first entry
    if (is.good())
    {
        token tok(is);

        if (is.good() && tok.isWord("FoamFile"))
        {
            dictionary headerDict(is);  // Read sub-dictionary content

            if (headerDict.readIfPresent("version", tok))
            {
                is.version(tok);
            }

            word formatName;
            if (headerDict.readIfPresent("format", formatName))
            {
                is.format(formatName);
            }

            //// Obtain number of blocks directly
            ///  This may not be reliable...
            //if (headerDict.readIfPresent("blocks", nBlocks))
            //{
            //    return nBlocks;
            //}
        }
        else if (tok.good())
        {
            is.putBack(tok);
        }
    }

    while (is.good() && skipBlockEntry(is))
    {
        ++nBlocks;

        if (maxNumBlocks == nBlocks)
        {
            break;
        }
    }

    return nBlocks;
}


bool Foam::decomposedBlockData::hasBlock(Istream& is, const label blockNumber)
{
    return
    (
        blockNumber >= 0
     && (blockNumber < getNumBlocks(is, blockNumber+1))
    );
}


std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    const label blocki,
    const UList<char>& charData
)
{
    // Offset to the beginning of this output
    // This should generally be OK for non-compressed streams
    // (eg, std::ofstream)

    std::streamoff blockOffset = os.stdStream().tellp();

    const word procName("processor" + Foam::name(blocki));

    // Write as commented content
    {
        os  << nl << "// " << procName << nl;
        charData.writeList(os) << nl;
    }
    // Write as primitiveEntry
    // {
    //     os << nl << procName << nl;
    //     charData.writeList(os);
    //     os.endEntry();
    // }

    return blockOffset;
}


std::streamoff Foam::decomposedBlockData::writeBlockEntry
(
    OSstream& os,
    IOstreamOption streamOptData,
    const regIOobject& io,
    const label blocki,
    const bool withLocalHeader
)
{
    // Serialize content to write
    List<char> contentChars;
    {
        OCharStream buf(streamOptData);

        bool ok = true;

        // Generate FoamFile header on master, without comment banner
        if (withLocalHeader)
        {
            const bool old = IOobject::bannerEnabled(false);

            ok = io.writeHeader(buf);

            IOobject::bannerEnabled(old);
        }

        // Write the data to the Ostream
        ok = ok && io.writeData(buf);

        if (!ok)
        {
            return std::streamoff(-1);
        }

        // Retrieve the output content
        buf.swap(contentChars);
    }

    return decomposedBlockData::writeBlockEntry(os, blocki, contentChars);
}


Foam::autoPtr<Foam::ISstream>
Foam::decomposedBlockData::readBlock
(
    const label blocki,
    ISstream& is,
    IOobject& headerIO
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlock:"
            << " stream:" << is.name() << " attempt to read block " << blocki
            << endl;
    }

    // Extracted header information
    IOstreamOption streamOptData;
    unsigned labelWidth = is.labelByteSize();
    unsigned scalarWidth = is.scalarByteSize();

    autoPtr<ISstream> realIsPtr;

    // Read master for header
    List<char> data;
    decomposedBlockData::readBlockEntry(is, data);

    if (blocki == 0)
    {
        realIsPtr.reset(new ICharStream(std::move(data)));
        realIsPtr->name() = is.name();

        {
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*realIsPtr))
            {
                FatalIOErrorInFunction(*realIsPtr)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
        }
    }
    else
    {
        {
            // Read header from first block
            ISpanStream headerStream(data);
            if (!headerIO.readHeader(headerStream))
            {
                FatalIOErrorInFunction(headerStream)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
            streamOptData = static_cast<IOstreamOption>(headerStream);
            labelWidth = headerStream.labelByteSize();
            scalarWidth = headerStream.scalarByteSize();
        }

        // Skip intermediate blocks
        for (label i = 1; i < blocki; ++i)
        {
            decomposedBlockData::skipBlockEntry(is);
        }

        // Read the block of interest
        decomposedBlockData::readBlockEntry(is, data);

        realIsPtr.reset(new ICharStream(std::move(data)));
        realIsPtr->name() = is.name();

        // Apply stream settings
        realIsPtr().format(streamOptData.format());
        realIsPtr().version(streamOptData.version());
        realIsPtr().setLabelByteSize(labelWidth);
        realIsPtr().setScalarByteSize(scalarWidth);
    }

    return realIsPtr;
}


bool Foam::decomposedBlockData::readBlocks
(
    const label comm,
    autoPtr<ISstream>& isPtr,
    List<char>& data,
    const UPstream::commsTypes  /* unused */
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "invalid")
            << " non-blocking comm:" << comm << endl;
    }

    // Read data on master and transmit. Always non-blocking

    bool ok = false;
    List<List<char>> procData;

    const label startOfRequests = UPstream::nRequests();

    if (UPstream::master(comm))
    {
        auto& is = *isPtr;
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, data);

        // Read proc data and setup non-blocking sends
        procData.resize(UPstream::nProcs(comm));
        for (const int proci : UPstream::subProcs(comm))
        {
            List<char>& procSlice = procData[proci];
            decomposedBlockData::readBlockEntry(is, procSlice);

            // Send content (non-blocking)
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                procSlice.cdata(),
                procSlice.size_bytes(),
                UPstream::msgType(),
                comm
            );
        }

        ok = is.good();
    }
    else if (UPstream::is_subrank(comm))
    {
        List<char>& procSlice = data;

        // Probe the incoming message source/size
        std::pair<int, int> probed =
            UPstream::probeMessage
            (
                UPstream::commsTypes::blocking,
                UPstream::masterNo(),
                UPstream::msgType(),
                comm
            );

        procSlice.resize_nocopy(probed.second);

        // Receive content (can also be zero-sized)
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::masterNo(),
            procSlice.data(),
            procSlice.size_bytes(),
            UPstream::msgType(),
            comm
        );
    }

    UPstream::waitRequests(startOfRequests);
    procData.clear();

    // Sync the status
    Pstream::broadcast(ok, comm);

    return ok;
}


Foam::autoPtr<Foam::ISstream> Foam::decomposedBlockData::readBlocks
(
    const label comm,
    const fileName& fName,
    autoPtr<ISstream>& isPtr,
    IOobject& headerIO,
    const UPstream::commsTypes  /* unused */
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::readBlocks:"
            << " stream:" << (isPtr ? isPtr->name() : "invalid")
            << " non-blocking" << endl;
    }

    // Read data on master and transmit. Always non-blocking
    bool ok = false;
    List<char> localData;
    List<List<char>> procData;
    autoPtr<ISstream> realIsPtr;

    const label startOfRequests = UPstream::nRequests();

    if (UPstream::master(comm))
    {
        auto& is = *isPtr;
        is.fatalCheck(FUNCTION_NAME);

        // Read master data
        decomposedBlockData::readBlockEntry(is, localData);

        // Move block data into a stream
        realIsPtr.reset(new ICharStream(std::move(localData)));
        realIsPtr->name() = fName;

        {
            // Read header from first block,
            // advancing the stream position
            if (!headerIO.readHeader(*realIsPtr))
            {
                FatalIOErrorInFunction(*realIsPtr)
                    << "Problem while reading object header "
                    << is.relativeName() << nl
                    << exit(FatalIOError);
            }
        }

        // Read proc data and setup non-blocking sends
        procData.resize(UPstream::nProcs(comm));
        for (const int proci : UPstream::subProcs(comm))
        {
            List<char>& procSlice = procData[proci];
            decomposedBlockData::readBlockEntry(is, procSlice);

            // Send content - non-blocking mode
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                proci,
                procSlice.cdata(),
                procSlice.size_bytes(),
                UPstream::msgType(),
                comm
            );
        }

        ok = is.good();
    }
    else if (UPstream::is_subrank(comm))
    {
        List<char>& procSlice = localData;

        // Probe the incoming message source/size
        std::pair<int, int> probed =
            UPstream::probeMessage
            (
                UPstream::commsTypes::blocking,
                UPstream::masterNo(),
                UPstream::msgType(),
                comm
            );

        procSlice.resize_nocopy(probed.second);

        // Receive content (can also be zero-sized)
        UIPstream::read
        (
            UPstream::commsTypes::nonBlocking,
            UPstream::masterNo(),
            procSlice.data(),
            procSlice.size_bytes(),
            UPstream::msgType(),
            comm
        );
    }

    UPstream::waitRequests(startOfRequests);
    procData.clear();

    if (UPstream::is_subrank(comm))
    {
        // Move block data into a stream
        realIsPtr.reset(new ICharStream(std::move(localData)));
        realIsPtr->name() = fName;
    }

    // Sync information

    // // Not needed...
    // Pstream::broadcast(ok, comm);

    // Broadcast master header info,
    // set stream properties from realIsPtr on master

    int verValue;
    int fmtValue;
    unsigned labelWidth;
    unsigned scalarWidth;
    word headerName(headerIO.name());

    if (UPstream::master(comm))
    {
        auto& is = realIsPtr();
        verValue = is.version().canonical();
        fmtValue = static_cast<int>(is.format());
        labelWidth = is.labelByteSize();
        scalarWidth = is.scalarByteSize();
    }

    Pstream::broadcasts
    (
        comm,
        verValue,
        fmtValue,
        labelWidth,
        scalarWidth,
        headerName,
        headerIO.headerClassName(),
        headerIO.note()
        // Unneeded: headerIO.instance()
        // Unneeded: headerIO.local()
    );

    if (realIsPtr)
    {
        auto& is = realIsPtr();
        is.version(IOstreamOption::versionNumber::canonical(verValue));
        is.format(IOstreamOption::streamFormat(fmtValue));
        is.setLabelByteSize(labelWidth);
        is.setScalarByteSize(scalarWidth);
    }

    headerIO.rename(headerName);

    if (debug)
    {
        Info<< "reading ok:" << ok << endl;
    }

    return realIsPtr;
}


void Foam::decomposedBlockData::gatherProcData
(
    const label comm,
    const UList<char>& localData,
    const labelUList& recvSizes,

    const labelRange& fromProcs,

    List<int>& sliceOffsets,
    DynamicList<char>& recvData,
    const UPstream::commsTypes commsType
)
{
    const label myRank = UPstream::myProcNo(comm);
    const label nProcs = UPstream::nProcs(comm);

    int nSendBytes = 0;
    recvData.clear();

    // On master, calculate sizing/offsets and resize the recv buffer
    List<int> sliceSizes;
    if (UPstream::master(comm))
    {
        sliceSizes.resize_nocopy(nProcs);
        sliceSizes = 0;
        sliceOffsets.resize_nocopy(nProcs+1);
        sliceOffsets = 0;

        int totalSize = 0;
        for (const label proci : fromProcs)
        {
            sliceSizes[proci] = int(recvSizes[proci]);
            sliceOffsets[proci] = totalSize;
            totalSize += sliceSizes[proci];
        }

        // One beyond the end of the range
        const label endProci = fromProcs.end_value();

        sliceOffsets[endProci] = totalSize;
        recvData.resize_nocopy(totalSize);
    }
    else if (fromProcs.contains(myRank) && !localData.empty())
    {
        // Note: UPstream::gather limited to int
        nSendBytes = int(localData.size_bytes());
    }


    if (UPstream::commsTypes::nonBlocking == commsType)
    {
        if (UPstream::master(comm))
        {
            for (const label proci : fromProcs)
            {
                SubList<char> procSlice
                (
                    recvData,
                    sliceOffsets[proci+1]-sliceOffsets[proci],
                    sliceOffsets[proci]
                );

                if (procSlice.empty())
                {
                    continue;
                }
                else if (proci == UPstream::masterNo())
                {
                    // No self-communication, although masterNo is usually
                    // not in the fromProcs list anyhow.
                    std::copy
                    (
                        localData.cbegin(),
                        localData.cbegin(procSlice.size()),
                        procSlice.begin()
                    );
                }
                else
                {
                    // Receive non-zero content
                    UIPstream::read
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        procSlice.data(),
                        procSlice.size_bytes(),
                        UPstream::msgType(),
                        comm
                    );
                }
            }
        }
        else if (fromProcs.contains(myRank) && !localData.empty())
        {
            // Send non-zero content
            UOPstream::write
            (
                UPstream::commsTypes::nonBlocking,
                UPstream::masterNo(),
                localData.cdata(),
                localData.size_bytes(),
                UPstream::msgType(),
                comm
            );
        }
    }
    else
    {
        // This is MPI_Gatherv()
        UPstream::gather
        (
            localData.cdata(),
            nSendBytes,

            recvData.data(),
            sliceSizes,
            sliceOffsets,
            comm
        );
    }
}


bool Foam::decomposedBlockData::writeBlocks
(
    const label comm,
    autoPtr<OSstream>& osPtr,
    List<std::streamoff>& blockOffset,
    const UList<char>& localData,

    const labelUList& recvSizes,
    const UList<stdFoam::span<char>>& procData,

    const UPstream::commsTypes commsType,
    const bool syncReturnState
)
{
    if (debug)
    {
        Pout<< "decomposedBlockData::writeBlocks:"
            << " stream:" << (osPtr ? osPtr->name() : "none")
            << " data:" << localData.size()
            << " (master only) procData:" << procData.size()
            << " commsType:" << UPstream::commsTypeNames[commsType] << endl;
    }

    const label nProcs = UPstream::nProcs(comm);

    bool ok = true;

    // Write master data
    if (UPstream::master(comm))
    {
        if (notNull(blockOffset))
        {
            // Recovery of blockOffset is optional
            blockOffset.resize(nProcs);
        }

        OSstream& os = osPtr();

        std::streamoff currOffset =
            decomposedBlockData::writeBlockEntry
            (
                os,
                UPstream::masterNo(),
                localData
            );

        if (UPstream::masterNo() < blockOffset.size())
        {
            blockOffset[UPstream::masterNo()] = currOffset;
        }

        // Write all pre-gathered proc data.
        if (procData.size())
        {
            for (label proci = 1; proci < nProcs; ++proci)
            {
                std::streamoff currOffset =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        procData[proci]
                    );

                if (proci < blockOffset.size())
                {
                    blockOffset[proci] = currOffset;
                }
            }
        }

        ok = os.good();
    }

    if (!procData.empty())
    {
        // Had pre-gathered proc data
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        if (UPstream::master(comm))
        {
            // Master data already written ...
            OSstream& os = osPtr();

            // Receive and write proc data
            label maxNonLocalSize = 0;
            for (label proci = 1; proci < nProcs; ++proci)
            {
                maxNonLocalSize = max(maxNonLocalSize, recvSizes[proci]);
            }

            DynamicList<char> recvData(maxNonLocalSize);
            for (label proci = 1; proci < nProcs; ++proci)
            {
                recvData.resize_nocopy(recvSizes[proci]);

                if (!recvData.empty())
                {
                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        proci,
                        recvData.data(),
                        recvData.size_bytes(),
                        UPstream::msgType(),
                        comm
                    );
                }

                std::streamoff currOffset =
                    decomposedBlockData::writeBlockEntry
                    (
                        os,
                        proci,
                        recvData
                    );

                if (proci < blockOffset.size())
                {
                    blockOffset[proci] = currOffset;
                }
            }

            ok = os.good();
        }
        else if (UPstream::is_subrank(comm) && !localData.empty())
        {
            UOPstream::write
            (
                UPstream::commsTypes::scheduled,
                UPstream::masterNo(),
                localData.cdata(),
                localData.size_bytes(),
                UPstream::msgType(),
                comm
            );
        }
    }
    else
    {
        // Master data already written ...

        // Find out how many ranks can be received into
        // maxMasterFileBufferSize

        const off_t maxBufferSize
        (
            fileOperations::masterUncollatedFileOperation::
            maxMasterFileBufferSize
        );

        List<int> sliceOffsets;
        DynamicList<char> recvData;

        for
        (
            labelRange fromProcs(1, nProcs-1);
            (fromProcs.start() < nProcs && fromProcs.size() > 0);
            (fromProcs.start() += fromProcs.size())
        )
        {
            {
                label endProci = fromProcs.start();

                if (UPstream::master(comm))
                {
                    // Send at least one proc, even if that if it is larger
                    // than maxBufferSize.
                    // Also handles the corner case the first proci has size 0,
                    // but the next one is too large.

                    endProci = fromProcs.start()+1;  // At least one proc

                    for
                    (
                        off_t total = recvSizes[fromProcs.start()];
                        (
                            endProci < nProcs
                         &&
                            (
                                !total
                             || (total + recvSizes[endProci] < maxBufferSize)
                            )
                        );
                        ++endProci
                    )
                    {
                        total += recvSizes[endProci];
                    }
                }

                Pstream::broadcast(endProci, comm);
                fromProcs.size() = (endProci - fromProcs.start());
            }

            const label startOfRequests = UPstream::nRequests();

            // Setup non-blocking send/recv or MPI_Gatherv
            gatherProcData
            (
                comm,
                localData,
                recvSizes,

                fromProcs,

                sliceOffsets,
                recvData,
                commsType  // ie, blocking or non-blocking
            );

            // For sanity checks
            // const label endOfRequests = UPstream::nRequests();

            if (UPstream::master(comm))
            {
                OSstream& os = osPtr();

                // Write received data
                label currRequest = startOfRequests;
                for (const label proci : fromProcs)
                {
                    SubList<char> procSlice
                    (
                        recvData,
                        sliceOffsets[proci+1]-sliceOffsets[proci],
                        sliceOffsets[proci]
                    );

                    if
                    (
                        (UPstream::commsTypes::nonBlocking == commsType)
                     && (proci != UPstream::masterNo())
                     && !procSlice.empty()
                    )
                    {
                        UPstream::waitRequest(currRequest);
                        ++currRequest;
                    }

                    std::streamoff currOffset =
                        decomposedBlockData::writeBlockEntry
                        (
                            os,
                            proci,
                            procSlice
                        );

                    if (proci < blockOffset.size())
                    {
                        blockOffset[proci] = currOffset;
                    }
                }
            }

            if (UPstream::commsTypes::nonBlocking == commsType)
            {
                UPstream::waitRequests(startOfRequests);
            }
        }

        if (UPstream::master(comm))
        {
            ok = osPtr->good();
        }
    }

    if (syncReturnState)
    {
        //- Enable to get synchronised error checking.
        //  Ensures that all procs are as slow as the master
        //  (which does all the writing)
        Pstream::broadcast(ok, comm);
    }

    return ok;
}


bool Foam::decomposedBlockData::read()
{
    autoPtr<ISstream> isPtr;
    fileName objPath(fileHandler().filePath(false, *this, word::null));
    if (UPstream::master(comm_))
    {
        isPtr.reset(new IFstream(objPath));
        IOobject::readHeader(*isPtr);
    }

    return readBlocks(comm_, isPtr, contentData_, commsType_);
}


bool Foam::decomposedBlockData::writeData(Ostream& os) const
{
    IOobject io(*this);
    IOstreamOption streamOpt(os);

    int verValue;
    int fmtValue;
    // Unneeded: word masterName(name());
    fileName masterLocation(instance()/db().dbDir()/local());

    // Re-read my own data to find out the header information
    if (UPstream::master(comm_))
    {
        ISpanStream headerStream(contentData_);
        io.readHeader(headerStream);

        verValue = headerStream.version().canonical();
        fmtValue = static_cast<int>(headerStream.format());
    }

    // Broadcast header information
    Pstream::broadcasts
    (
        comm_,
        verValue,
        fmtValue,
        // Unneeded: masterName
        io.headerClassName(),
        io.note(),
        // Unneeded: io.instance()
        // Unneeded: io.local()
        masterLocation
    );

    streamOpt.version(IOstreamOption::versionNumber::canonical(verValue));
    streamOpt.format(IOstreamOption::streamFormat(fmtValue));

    if (UPstream::is_subrank(comm_))
    {
        decomposedBlockData::writeHeader
        (
            os,
            streamOpt,  // streamOpt for data
            io.headerClassName(),
            io.note(),
            masterLocation,
            name(),
            dictionary()
        );
    }

    // Write the character data
    if (isA<OFstream>(os))
    {
        // Serial file output - can use writeRaw()
        os.writeRaw(contentData_.cdata(), contentData_.size_bytes());
    }
    else
    {
        // Other cases are less fortunate, and no std::string_view
        std::string str(contentData_.cdata(), contentData_.size_bytes());
        os.writeQuoted(str, false);
    }

    if (!Pstream::master(comm_))
    {
        IOobject::writeEndDivider(os);
    }

    return os.good();
}


bool Foam::decomposedBlockData::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    autoPtr<OSstream> osPtr;
    if (UPstream::master(comm_))
    {
        // Note: always write binary. These are strings so readable anyway.
        //       They have already be tokenised on the sending side.

        osPtr.reset(new OFstream(objectPath(), IOstreamOption::BINARY));

        // Update meta-data for current state
        const_cast<regIOobject&>
        (
            static_cast<const regIOobject&>(*this)
        ).updateMetaData();

        decomposedBlockData::writeHeader
        (
            *osPtr,
            streamOpt,  // streamOpt for data
            static_cast<const IOobject&>(*this)
        );
    }

    const labelList recvSizes
    (
        UPstream::listGatherValues<label>(contentData_.size(), comm_)
    );

    List<std::streamoff> blockOffsets;  // Optional
    return writeBlocks
    (
        comm_,
        osPtr,
        blockOffsets,
        contentData_,
        recvSizes,
        UList<stdFoam::span<char>>(),  // dummy proc data (nothing pre-gathered)
        commsType_
    );
}


// ************************************************************************* //
