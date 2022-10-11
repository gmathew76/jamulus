/******************************************************************************\
 * Copyright (c) 2004-2022
 *
 * Author(s):
 *  Volker Fischer
 *
 * Note: We are assuming here that put and get operations are secured by a mutex
 *       and accessing does not occur at the same time.
 *
 ******************************************************************************
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 *
\******************************************************************************/

#include "buffer.h"

/* Network buffer implementation **********************************************/
void CNetBuf::Init ( const int iNewBlockSize, const int iNewNumBlocks, const bool bNUseSequenceNumber, const bool bPreserve )
{
    // store the sequence number activation flag
    bUseSequenceNumber = bNUseSequenceNumber;

    // in simulation mode the size is not changed during operation -> we do
    // not have to implement special code for this case
    // only enter the "preserve" branch, if object was already initialized
    // and the block sizes are the same
    if ( bPreserve && ( !bIsSimulation ) && bIsInitialized && ( iBlockSize == iNewBlockSize ) )
    {
        // extract all data from buffer in temporary storage
        CVector<CVector<uint8_t>> vecvecTempMemory = vecvecMemory; // allocate worst case memory by copying

        if ( !bNUseSequenceNumber )
        {
            int iPreviousDataCnt = 0;

            while ( Get ( vecvecTempMemory[iPreviousDataCnt], iBlockSize ) )
            {
                iPreviousDataCnt++;
            }

            // now resize the buffer to the new size (buffer is empty after this operation)
            Resize ( iNewNumBlocks, iNewBlockSize );

            // copy the previous data back in the buffer (make sure we only copy as much
            // data back as the new buffer size can hold)
            int iDataCnt = 0;

            while ( ( iDataCnt < iPreviousDataCnt ) && Put ( vecvecTempMemory[iDataCnt], iBlockSize ) )
            {
                iDataCnt++;
            }
        }
        else
        {
            // store current complete buffer state in temporary memory
            CVector<int>  veciTempBlockValid ( iNumBlocksMemory );
            const uint8_t iOldSequenceNumberAtGetPos = iSequenceNumberAtGetPos;
            const int     iOldNumBlocksMemory        = iNumBlocksMemory;
            const int     iOldBlockGetPos            = iBlockGetPos;
            const int     iOldBlockPutPos            = iBlockPutPos;
            const EBufState oldEBufState             = eBufState;
            int           iCurBlockPos               = 0;

            int iDeltaPutGetPos = DeltaDistance(iOldBlockPutPos, iOldBlockGetPos);
            int iOldCount = (eBufState == BS_FULL) ? iOldNumBlocksMemory : iDeltaPutGetPos;
            int iNewCount = iOldCount;

            printf("Stats fillGets %d, lowBuffers %d, outOfSeqBuffers %d, outOfSeqRangeBuffers %d\n", fillGets, lowBuffers, outOfSeqBuffers, outOfSeqRangeBuffers );
            printf("Resize oldNumBlocks %d newNumBlocks %d iOldCount = %d, newCount = %d, delta = %d\n", iOldNumBlocksMemory, iNewNumBlocks, iOldCount, iNewCount, iDeltaPutGetPos);

            while ( iBlockGetPos < iNumBlocksMemory )
            {
                veciTempBlockValid[iCurBlockPos] = veciBlockValid[iBlockGetPos];
                vecvecTempMemory[iCurBlockPos++] = vecvecMemory[iBlockGetPos++];
            }

            for ( iBlockGetPos = 0; iBlockGetPos < iOldBlockGetPos; iBlockGetPos++ )
            {
                veciTempBlockValid[iCurBlockPos] = veciBlockValid[iBlockGetPos];
                vecvecTempMemory[iCurBlockPos++] = vecvecMemory[iBlockGetPos];
            }

            printf("Copied %d blocks.\n", iCurBlockPos);

            // now resize the buffer to the new size
            Resize ( iNewNumBlocks, iNewBlockSize );

            //
            // Copy over the data from the temporary vector back into the new buffer.
            //

            iBlockGetPos            = 0; // per definition

            //
            // If the new size is smaller than the number of blocks we have, we discard blocks starting
            // at the get position until we can fit the new buffer.
            //
            if (iOldCount <= iNewNumBlocks) {
                
                iSequenceNumberAtGetPos = iOldSequenceNumberAtGetPos;
                iNewCount = iOldCount;

                // copy over everything
                for ( int iCurPos = 0; iCurPos < iOldCount; iCurPos++ )
                {
                    veciBlockValid[iCurPos] = veciTempBlockValid[iCurPos];
                    vecvecMemory[iCurPos]   = vecvecTempMemory[iCurPos];
                }                
                iNewCount = iOldCount;
                iBlockPutPos = iNewCount;
            }
            else 
            {
                iSequenceNumberAtGetPos = iOldSequenceNumberAtGetPos + (iOldCount - iNewNumBlocks);
                iNewCount = iNewNumBlocks;

                // copy over everything from the tail of the buffer.
                int srcIdx, dstIdx;
                for ( dstIdx = 0, srcIdx = (iOldCount - iNewCount); dstIdx < iNewCount; srcIdx++, dstIdx++)
                {
                    veciBlockValid[dstIdx] = veciTempBlockValid[srcIdx];
                    vecvecMemory[dstIdx]   = vecvecTempMemory[srcIdx];
                }
                iBlockPutPos = 0; // a full buffer.
            }

            //
            // Update the buffer state.
            //

            switch (oldEBufState)
            {
                case BS_EMPTY:
                    break;

                case BS_FILLING:

                    if (iNewCount >= iNewNumBlocks)
                    {
                       eBufState = BS_FULL; 
                    }
                    else if (iNewCount > iNumBlocksMemory / 2)
                    {
                       eBufState = BS_OK; 
                    }
                    else if (iNewCount == 0)
                    {
                        eBufState = BS_EMPTY;
                    }
                    break;

                case BS_OK:

                    if (iNewCount >= iNumBlocksMemory)
                    {
                       eBufState = BS_FULL; 
                    }
                    else if (iNewCount <= iNumBlocksMemory / 2)
                    {
                       eBufState = BS_FILLING; 
                    }
                    else if (iNewCount == 0)
                    {
                        eBufState = BS_EMPTY;
                    }
                    break;

                case BS_FULL:

                    if (iNewCount == 0)
                    {
                        eBufState = BS_EMPTY;
                    }
                    else if (iNewCount <= iNumBlocksMemory / 2)
                    {
                       eBufState = BS_FILLING; 
                    }
                    else if (iNewCount < iNumBlocksMemory)
                    {
                       eBufState = BS_OK;
                    }
                    break;
            }

        }
    }
    else
    {
        Resize ( iNewNumBlocks, iNewBlockSize );
        eBufState = BS_EMPTY;
    }

    // set initialized flag
    bIsInitialized = true;
}

void CNetBuf::Resize ( const int iNewNumBlocks, const int iNewBlockSize )
{
    // allocate memory for actual data buffer
    vecvecMemory.Init ( iNewNumBlocks );
    veciBlockValid.Init ( iNewNumBlocks, 0 ); // initialize with zeros = invalid

    if ( !bIsSimulation )
    {
        for ( int iBlock = 0; iBlock < iNewNumBlocks; iBlock++ )
        {
            vecvecMemory[iBlock].Init ( iNewBlockSize );
        }
    }

    // init buffer pointers and buffer state (empty buffer) and store buffer properties
    iBlockGetPos     = 0;
    iBlockPutPos     = 0;
    eBufState        = BS_EMPTY;
    iBlockSize       = iNewBlockSize;
    iNumBlocksMemory = iNewNumBlocks;
}

bool CNetBuf::Put ( const CVector<uint8_t>& vecbyData, int iInSize )
{
    // if the sequence number is used, we need a complete different way of applying
    // the new network packet
    if ( bUseSequenceNumber )
    {
        // check that the input size is a multiple of the block size
        if ( ( iInSize % ( iBlockSize + iNumBytesSeqNum ) ) != 0 )
        {
            return false;
        }

        // to get the number of input blocks we assume that the number of bytes for
        // the sequence number is much smaller than the number of coded audio bytes
        const int iNumBlocks = /* floor */ ( iInSize / iBlockSize );

        // copy new data in internal buffer
        for ( int iBlock = 0; iBlock < iNumBlocks; iBlock++ )
        {

            int iNewBlockPutPos = iBlockPutPos;
            bool bUpdatePutPos = true;
            const int iBlockOffset = iBlock * ( iBlockSize + iNumBytesSeqNum );

            // extract sequence number of current received block (per definition
            // the sequence number is appended after the coded audio data)
            const int iCurrentSequenceNumber = vecbyData[iBlockOffset + iBlockSize];

            // calculate the sequence number difference and take care of wrap
            int iSeqNumDiff = iCurrentSequenceNumber - static_cast<int> ( iSequenceNumberAtGetPos );

            if ( iSeqNumDiff < -128 )
            {
                iSeqNumDiff += 256;
            }
            else if ( iSeqNumDiff >= 128 )
            {
                iSeqNumDiff -= 256;
            }

            // The 1-byte sequence number wraps around at a count of 256. So, if a packet is delayed
            // further than this we cannot detect it. But it does not matter since such a packet is
            // more than 100 ms delayed so we have a bad network situation anyway. Therefore we
            // assume that the sequence number difference between the received and local counter is
            // correct. 
            //

            //
            // If the received packet is too late (i.e. our get position is already ahead), we don't have
            // any use for the data. Drop it and move to the next block.
            //
            if ( iSeqNumDiff < 0 )
            {
                continue;
            }
            //
            // else if the received packet comes too early so we move the "buffer window" in the
            // future until the received packet is the last packet in the buffer
            //
            else if ( iSeqNumDiff >= iNumBlocksMemory )
            {
                // the received packet comes too early so we move the "buffer window" in the
                // future until the received packet is the last packet in the buffer
                for ( int i = 0; i < iSeqNumDiff - iNumBlocksMemory + 1; i++ )
                {
                    // insert an invalid block at the shifted position
                    veciBlockValid[iBlockGetPos] = 0; // invalidate

                    // we increase the local sequence number and get position and take care of wrap
                    iSequenceNumberAtGetPos++;
                    iBlockGetPos++;

                    if ( iBlockGetPos >= iNumBlocksMemory )
                    {
                        iBlockGetPos -= iNumBlocksMemory;
                    }
                }

                outOfSeqRangeBuffers++;

                // insert the new packet at the end of the buffer since it is too early (since
                // we add an offset to the get position, we have to take care of wrapping)
                iNewBlockPutPos = iBlockGetPos + iNumBlocksMemory - 1;

                if ( iNewBlockPutPos >= iNumBlocksMemory )
                {
                    iNewBlockPutPos -= iNumBlocksMemory;
                }
            }
            else
            {
                // this is the regular case: the received packet fits into the buffer so
                // we will write it at the correct position based on the sequence number
                iNewBlockPutPos = iBlockGetPos + iSeqNumDiff;

                if ( iNewBlockPutPos >= iNumBlocksMemory )
                {
                    iNewBlockPutPos -= iNumBlocksMemory;
                }

                if (DeltaDistance(iNewBlockPutPos, iBlockGetPos) < DeltaDistance(iBlockPutPos, iBlockGetPos))
                {
                    outOfSeqBuffers++;
                    bUpdatePutPos = false;
                }
            }

            // for simulation buffer only update pointer, no data copying
            if ( !bIsSimulation )
            {
                // copy one block of data in buffer
                std::copy ( vecbyData.begin() + iBlockOffset,
                            vecbyData.begin() + iBlockOffset + iBlockSize,
                            vecvecMemory[iNewBlockPutPos].begin() );
            }

            // valid packet added, set flag
            veciBlockValid[iNewBlockPutPos] = 1;

            //
            // Update our put position if we received a packet at the end of the buffer.
            //

            if (bUpdatePutPos)
            {
                iBlockPutPos = iNewBlockPutPos + 1;
                if ( iBlockPutPos >= iNumBlocksMemory )
                {
                    iBlockPutPos -= iNumBlocksMemory;
                }
            }

            //
            // Update the buffer state.
            //

            int iPutGetDelta = DeltaDistance(iBlockPutPos, iBlockGetPos);

            switch (eBufState)
            {
                case BS_EMPTY:

                    if (iPutGetDelta == 0)
                    {
                        eBufState = BS_FULL;
                    }
                    else if (iPutGetDelta <= iNumBlocksMemory/2)
                    {
                        eBufState = BS_FILLING;
                    }
                    else 
                    {
                        eBufState = BS_OK;
                    }
                    break;

                case BS_FILLING:

                    if (iPutGetDelta == 0)
                    {
                        eBufState = BS_FULL;
                    }
                    else if (iPutGetDelta > iNumBlocksMemory/2)
                    {
                        eBufState = BS_OK;
                    }
                    break;

                case BS_OK:
                    if (iPutGetDelta == 0)
                    {
                        eBufState = BS_FULL;
                    }
                    break;

                case BS_FULL:
                    break;
            }

        } // for
    }
    else
    {
        // check if there is not enough space available and that the input size is a
        // multiple of the block size
        if ( ( GetAvailSpace() < iInSize ) || ( ( iInSize % iBlockSize ) != 0 ) )
        {
            return false;
        }

        // copy new data in internal buffer
        const int iNumBlocks = iInSize / iBlockSize;

        for ( int iBlock = 0; iBlock < iNumBlocks; iBlock++ )
        {
            // for simultion buffer only update pointer, no data copying
            if ( !bIsSimulation )
            {
                // calculate the block offset once per loop instead of repeated multiplying
                const int iBlockOffset = iBlock * iBlockSize;

                // copy one block of data in buffer
                std::copy ( vecbyData.begin() + iBlockOffset, vecbyData.begin() + iBlockOffset + iBlockSize, vecvecMemory[iBlockPutPos].begin() );
            }

            // set the put position one block further
            iBlockPutPos++;

            // take care about wrap around of put pointer
            if ( iBlockPutPos == iNumBlocksMemory )
            {
                iBlockPutPos = 0;
            }
        }

        // set buffer state

        const int iDeltaPutGetPos = DeltaDistance(iBlockPutPos, iBlockGetPos);

        switch (eBufState)
        {
            case BS_EMPTY:

                if (iDeltaPutGetPos == 0)
                {
                    eBufState = BS_FULL;
                }
                else if (iDeltaPutGetPos <= iNumBlocksMemory/2)
                {
                    eBufState = BS_FILLING;
                }
                else
                {
                    eBufState = BS_OK;
                }
                break;

            case BS_FILLING:

                if (iDeltaPutGetPos == 0)
                {
                    eBufState = BS_FULL;
                }
                else if (iDeltaPutGetPos > iNumBlocksMemory/2)
                {
                    eBufState = BS_OK;
                }
                break;

            case BS_OK:

                if (iDeltaPutGetPos == 0)
                {
                    eBufState = BS_FULL;
                }

            case BS_FULL:
                break;
        }

    }

    return true;
}

bool CNetBuf::Get ( CVector<uint8_t>& vecbyData, const int iOutSize )
{
    bool bReturn = true;

    // check requested output size and available buffer data
    if ( ( iOutSize == 0 ) || ( iOutSize != iBlockSize ) || ( GetAvailData() < iOutSize ) )
    {
        return false;
    }

    if ((eBufState == BS_FILLING) || (eBufState == BS_EMPTY)) {
        fillGets++;
        return false;
    }

    // if using sequence numbers, we do not use the block put position
    // at all but only determine the state from the "valid block" indicator
    if ( bUseSequenceNumber )
    {
        bReturn = ( veciBlockValid[iBlockGetPos] > 0 );

        // invalidate the block we are now taking from the buffer
        veciBlockValid[iBlockGetPos] = 0; // zero means invalid
    }

    // for simultion buffer or invalid block only update pointer, no data copying
    if ( !bIsSimulation && bReturn )
    {
        // copy data from internal buffer in output buffer
        std::copy ( vecvecMemory[iBlockGetPos].begin(), vecvecMemory[iBlockGetPos].begin() + iBlockSize, vecbyData.begin() );
    }

    // set the get position and sequence number one block further
    iBlockGetPos++;
    iSequenceNumberAtGetPos++; // wraps around automatically

    // take care about wrap around of get pointer
    if ( iBlockGetPos == iNumBlocksMemory )
    {
        iBlockGetPos = 0;
    }

    int delta = DeltaDistance(iBlockPutPos,iBlockGetPos);

    if (delta <= 1) {
        lowBuffers++;
    }

    // set buffer state flag
    switch(eBufState)
    {
        case BS_FULL:

            if (delta == 0)
            {
                eBufState = BS_EMPTY;
            }
            else{
                eBufState = BS_OK;
            }
            break;

        case BS_OK:
        case BS_FILLING:

            if ( delta == 0 )
            {
                eBufState = BS_EMPTY;
            }
            break;

        default:
            break;    
    }

    if ( iBlockPutPos == iBlockGetPos )
    {
        eBufState = BS_EMPTY;
    }
    else {
        eBufState = BS_OK;
    }

    return bReturn;
}

int CNetBuf::GetAvailSpace() const
{
    return iNumBlocksMemory * iBlockSize - GetAvailData();
}

int CNetBuf::GetAvailData() const
{
    int iAvBlocks = DeltaDistance(iBlockPutPos, iBlockGetPos);

    if ( ( iAvBlocks == 0 ) && ( eBufState == BS_FULL ) )
    {
        iAvBlocks = iNumBlocksMemory;
    }

    return iAvBlocks * iBlockSize;
}

int CNetBuf::DeltaDistance(int putPos, int getPos) const
{
    if (putPos >= getPos)
    {
        return putPos - getPos;
    }
    else
    {
        return iNumBlocksMemory - getPos + putPos;
    }
}    

/* Network buffer with statistic calculations implementation ******************/
CNetBufWithStats::CNetBufWithStats() :
    CNetBuf ( false ), // base class init: no simulation mode
    iMaxStatisticCount ( MAX_STATISTIC_COUNT ),
    bUseDoubleSystemFrameSize ( false ),
    dAutoFilt_WightUpNormal ( IIR_WEIGTH_UP_NORMAL ),
    dAutoFilt_WightDownNormal ( IIR_WEIGTH_DOWN_NORMAL ),
    dAutoFilt_WightUpFast ( IIR_WEIGTH_UP_FAST ),
    dAutoFilt_WightDownFast ( IIR_WEIGTH_DOWN_FAST ),
    dErrorRateBound ( ERROR_RATE_BOUND ),
    dUpMaxErrorBound ( UP_MAX_ERROR_BOUND )
{
    // Define the sizes of the simulation buffers,
    // must be NUM_STAT_SIMULATION_BUFFERS elements!
    // Avoid the buffer length 1 because we do not have a solution for a
    // sample rate offset correction. Caused by the jitter we usually get bad
    // performance with just one buffer.
    viBufSizesForSim[0] = 2;
    viBufSizesForSim[1] = 3;
    viBufSizesForSim[2] = 4;
    viBufSizesForSim[3] = 5;
    viBufSizesForSim[4] = 6;
    viBufSizesForSim[5] = 7;
    viBufSizesForSim[6] = 8;
    viBufSizesForSim[7] = 9;
    viBufSizesForSim[8] = 10;
    viBufSizesForSim[9] = 11;

    // set all simulation buffers in simulation mode
    for ( int i = 0; i < NUM_STAT_SIMULATION_BUFFERS; i++ )
    {
        SimulationBuffer[i].SetIsSimulation ( true );
    }
}

void CNetBufWithStats::GetErrorRates ( CVector<double>& vecErrRates, double& dLimit, double& dMaxUpLimit )
{
    // get all the averages of the error statistic
    vecErrRates.Init ( NUM_STAT_SIMULATION_BUFFERS );

    for ( int i = 0; i < NUM_STAT_SIMULATION_BUFFERS; i++ )
    {
        vecErrRates[i] = ErrorRateStatistic[i].GetAverage();
    }

    // get the limits for the decisions
    dLimit      = dErrorRateBound;
    dMaxUpLimit = dUpMaxErrorBound;
}

void CNetBufWithStats::Init ( const int iNewBlockSize, const int iNewNumBlocks, const bool bNUseSequenceNumber, const bool bPreserve )
{
    // call base class Init
    CNetBuf::Init ( iNewBlockSize, iNewNumBlocks, bNUseSequenceNumber, bPreserve );

    // inits for statistics calculation
    if ( !bPreserve )
    {
        // set the auto filter weights and max statistic count
        if ( bUseDoubleSystemFrameSize )
        {
            dAutoFilt_WightUpNormal   = IIR_WEIGTH_UP_NORMAL_DOUBLE_FRAME_SIZE;
            dAutoFilt_WightDownNormal = IIR_WEIGTH_DOWN_NORMAL_DOUBLE_FRAME_SIZE;
            dAutoFilt_WightUpFast     = IIR_WEIGTH_UP_FAST_DOUBLE_FRAME_SIZE;
            dAutoFilt_WightDownFast   = IIR_WEIGTH_DOWN_FAST_DOUBLE_FRAME_SIZE;
            iMaxStatisticCount        = MAX_STATISTIC_COUNT_DOUBLE_FRAME_SIZE;
            dErrorRateBound           = ERROR_RATE_BOUND_DOUBLE_FRAME_SIZE;
            dUpMaxErrorBound          = UP_MAX_ERROR_BOUND_DOUBLE_FRAME_SIZE;
        }
        else
        {
            dAutoFilt_WightUpNormal   = IIR_WEIGTH_UP_NORMAL;
            dAutoFilt_WightDownNormal = IIR_WEIGTH_DOWN_NORMAL;
            dAutoFilt_WightUpFast     = IIR_WEIGTH_UP_FAST;
            dAutoFilt_WightDownFast   = IIR_WEIGTH_DOWN_FAST;
            iMaxStatisticCount        = MAX_STATISTIC_COUNT;
            dErrorRateBound           = ERROR_RATE_BOUND;
            dUpMaxErrorBound          = UP_MAX_ERROR_BOUND;
        }

        for ( int i = 0; i < NUM_STAT_SIMULATION_BUFFERS; i++ )
        {
            // init simulation buffers with the correct size
            SimulationBuffer[i].Init ( iNewBlockSize, viBufSizesForSim[i], bNUseSequenceNumber );

            // init statistics
            ErrorRateStatistic[i].Init ( iMaxStatisticCount, true );
        }

        // reset the initialization counter which controls the initialization
        // phase length
        ResetInitCounter();

        // init auto buffer setting with a meaningful value, also init the
        // IIR parameter with this value
        iCurAutoBufferSizeSetting = 6;
        dCurIIRFilterResult       = iCurAutoBufferSizeSetting;
        iCurDecidedResult         = iCurAutoBufferSizeSetting;
    }
}

void CNetBufWithStats::ResetInitCounter()
{
    // start initialization phase of IIR filtering, use a quarter the size
    // of the error rate statistic buffers which should be ok for a good
    // initialization value (initialization phase should be as short as
    // possible)
    iInitCounter = iMaxStatisticCount / 4;
}

bool CNetBufWithStats::Put ( const CVector<uint8_t>& vecbyData, const int iInSize )
{
    // call base class Put
    const bool bPutOK = CNetBuf::Put ( vecbyData, iInSize );

    // update statistics calculations
    for ( int i = 0; i < NUM_STAT_SIMULATION_BUFFERS; i++ )
    {
        ErrorRateStatistic[i].Update ( !SimulationBuffer[i].Put ( vecbyData, iInSize ) );
    }

    return bPutOK;
}

bool CNetBufWithStats::Get ( CVector<uint8_t>& vecbyData, const int iOutSize )
{
    // call base class Get
    const bool bGetOK = CNetBuf::Get ( vecbyData, iOutSize );

    // update statistics calculations
    for ( int i = 0; i < NUM_STAT_SIMULATION_BUFFERS; i++ )
    {
        ErrorRateStatistic[i].Update ( !SimulationBuffer[i].Get ( vecbyData, iOutSize ) );
    }

    // update auto setting
    UpdateAutoSetting();

    return bGetOK;
}

void CNetBufWithStats::UpdateAutoSetting()
{
    int  iCurDecision      = 0; // dummy initialization
    int  iCurMaxUpDecision = 0; // dummy initialization
    bool bDecisionFound;

    // Get regular error rate decision -----------------------------------------
    // Use a specified error bound to identify the best buffer size for the
    // current network situation. Start with the smallest buffer and
    // test for the error rate until the rate is below the bound.
    bDecisionFound = false;

    for ( int i = 0; i < NUM_STAT_SIMULATION_BUFFERS - 1; i++ )
    {
        if ( ( !bDecisionFound ) && ( ErrorRateStatistic[i].GetAverage() <= dErrorRateBound ) )
        {
            iCurDecision   = viBufSizesForSim[i];
            bDecisionFound = true;
        }
    }

    if ( !bDecisionFound )
    {
        // in case no buffer is below bound, use largest buffer size
        iCurDecision = viBufSizesForSim[NUM_STAT_SIMULATION_BUFFERS - 1];
    }

    // Get maximum upper error rate decision -----------------------------------
    // Use a specified error bound to identify the maximum upper error rate
    // to identify if we have a too low buffer setting which gives a very
    // bad performance constantly. Start with the smallest buffer and
    // test for the error rate until the rate is below the bound.
    bDecisionFound = false;

    for ( int i = 0; i < NUM_STAT_SIMULATION_BUFFERS - 1; i++ )
    {
        if ( ( !bDecisionFound ) && ( ErrorRateStatistic[i].GetAverage() <= dUpMaxErrorBound ) )
        {
            iCurMaxUpDecision = viBufSizesForSim[i];
            bDecisionFound    = true;
        }
    }

    if ( !bDecisionFound )
    {
        // in case no buffer is below bound, use largest buffer size
        iCurMaxUpDecision = viBufSizesForSim[NUM_STAT_SIMULATION_BUFFERS - 1];

        // This is a worst case, something very bad had happened. Hopefully
        // this was just temporary so that we initiate a new initialization
        // phase to get quickly back to normal buffer sizes (hopefully).
        ResetInitCounter();
    }

    // Post calculation (filtering) --------------------------------------------
    // Define different weights for up and down direction. Up direction
    // filtering shall be slower than for down direction since we assume
    // that the lower value is the actual value which can be used for
    // the current network condition. If the current error rate estimation
    // is higher, it may be a temporary problem which should not change
    // the current jitter buffer size significantly.
    // For the initialization phase, use lower weight values to get faster
    // adaptation.
    double       dWeightUp, dWeightDown;
    const double dHysteresisValue   = FILTER_DECISION_HYSTERESIS;
    bool         bUseFastAdaptation = false;

    // check for initialization phase
    if ( iInitCounter > 0 )
    {
        // decrease init counter
        iInitCounter--;

        // use the fast adaptation
        bUseFastAdaptation = true;
    }

    // if the current detected buffer setting is below the maximum upper bound
    // decision, then we enable a booster to go up to the minimum required
    // number of buffer blocks (i.e. we use weights for fast adaptation)
    if ( iCurAutoBufferSizeSetting < iCurMaxUpDecision )
    {
        bUseFastAdaptation = true;
    }

    if ( bUseFastAdaptation )
    {
        dWeightUp   = dAutoFilt_WightUpFast;
        dWeightDown = dAutoFilt_WightDownFast;
    }
    else
    {
        dWeightUp   = dAutoFilt_WightUpNormal;
        dWeightDown = dAutoFilt_WightDownNormal;
    }

    // apply non-linear IIR filter
    MathUtils().UpDownIIR1 ( dCurIIRFilterResult, static_cast<double> ( iCurDecision ), dWeightUp, dWeightDown );

    //### TEST: BEGIN ###//
    // TEST store important detection parameters in file for debugging
    /*
    static FILE* pFile = fopen ( "test.dat", "w" );
    static int icnt = 0;
    if ( icnt == 50 )
    {
        fprintf ( pFile, "%d %e\n", iCurDecision, dCurIIRFilterResult );
        fflush ( pFile );
        icnt = 0;
    }
    else
    {
        icnt++;
    }
    */
    //### TEST: END ###//

    // apply a hysteresis
    iCurAutoBufferSizeSetting = MathUtils().DecideWithHysteresis ( dCurIIRFilterResult, iCurDecidedResult, dHysteresisValue );

    // Initialization phase check and correction -------------------------------
    // sometimes in the very first period after a connection we get a bad error
    // rate result -> delete this from the initialization phase
    if ( iInitCounter == iMaxStatisticCount / 8 )
    {
        // check error rate of the largest buffer as the indicator
        if ( ErrorRateStatistic[NUM_STAT_SIMULATION_BUFFERS - 1].GetAverage() > dErrorRateBound )
        {
            for ( int i = 0; i < NUM_STAT_SIMULATION_BUFFERS; i++ )
            {
                ErrorRateStatistic[i].Reset();
            }
        }
    }
}


bool CNetBuf::TestImpl(void)
{
    CVector<uint8_t> testData;
    CVector<uint8_t>  outData;

    testData.Init(65);
    outData.Init(64);

    Init ( 64, 5, true, false );

    if (eBufState != BS_EMPTY) {
        printf("x. Unexpected State %d\n", eBufState);
        return false;
    }

    if (iBlockPutPos != 0) {
        printf("x. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("x. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 0) {
        printf("x. Unexpected deltaDistance %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }        

    testData[0] = 0;
    testData[64] = 0;
    if (!Put ( testData, 65)) {
        printf("Error put[0]\n");
        return false;
    }

    if (eBufState != BS_FILLING) {
        printf("0. Unexpected State %d", eBufState);
        return false;
    }

    if (iBlockPutPos != 1) {
        printf("0. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("0. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 1) {
        printf("0. Unexpected deltaDistance %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }        

   if (Get ( outData, 64)) {
       printf("Unexpected get[0_0]\n");
       return false;
   }

    testData[0] = 1;
    testData[64] = 1;
    if (!Put ( testData, 65)) {
        printf("Error put[1]\n");
        return false;
    }

    if (eBufState != BS_FILLING) {
        printf("1. Unexpected State %d\n", eBufState);
        return false;
    }

    if (iBlockPutPos != 2) {
        printf("1. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("1. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 2) {
        printf("1. Unexpected deltaDistance %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }    

   if (Get ( outData, 64)) {
       printf("Unexpected success get[0_1]\n");
       return false;
   }

    testData[0] = 2;
    testData[64] = 2;
    if (!Put ( testData, 65)) {
        printf("Error put[2]\n");
        return false;
    }

    if (eBufState != BS_OK) {
        printf("2. Unexpected State %d\n", eBufState);
        return false;
    }

    if (iBlockPutPos != 3) {
        printf("2. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("2. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 3) {
        printf("2. Unexpected deltaDistance %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }

   if (!Get ( outData, 64)) {
        printf("Unexpected get error get[0]\n");
        return false;
   }
   if (outData[0] != 0)
   {
        printf("Unexpected data returned. get[0]\n");
        return false;
   }

    if (eBufState != BS_OK) {
        printf("2_0. Unexpected State %d\n", eBufState);
        return false;
    }

    if (iBlockPutPos != 3) {
       printf("2_0. Unexpected iBlockPutPos %d\n", iBlockPutPos);
       return false;
    }
    if (iBlockGetPos != 1) {
        printf("2_0. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 2) {
        printf("2_0. Unexpected deltaDistance %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }

    testData[0] = 3;
    testData[64] = 3;
    if (!Put ( testData, 65)) {
        printf("Error put[3]\n");
        return false;
    }
    testData[0] = 4;
    testData[64] = 4;
    if (!Put ( testData, 65)) {
        printf("Error put[4]\n");
        return false;
    }

    if (eBufState != BS_OK) {
        printf("4. Unexpected State %d\n", eBufState);
        return false;
    }

    if (iBlockPutPos != 0) {
        printf("4. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 1) {
        printf("4. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 4) {
        printf("5. Unexpected deltaDistance %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }

    testData[0] = 5;
    testData[64] = 5;

    if (!Put ( testData, 65)) {
        printf("5. Unexpected error put[5]\n");
        return false;
    }

    if (eBufState != BS_FULL) {
        printf("5. Unexpected State %d\n", eBufState);
        return false;
    }

    if (iBlockPutPos != 1) {
        printf("5. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 1) {
        printf("5. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 0) {
        printf("5. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }

   if (!Get ( outData, 64)) {
        printf("Unexpected get error get[1]\n");
        return false;
   }

    if (iBlockPutPos != 1) {
        printf("5_1. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 2) {
        printf("5_1. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_OK) {
        printf("5_1. Unexpected State %d\n", eBufState);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 4) {
        printf("5_1. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }    

   if (outData[0] != 1)
   {
        printf("Unexpected data returned. get[1]\n");
        return false;
   }

   if (!Get ( outData, 64)) {
        printf("Unexpected get error get[2]\n");
        return false;
   }

    if (iBlockPutPos != 1) {
        printf("5_2. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 3) {
        printf("5_2. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_OK) {
        printf("5_2. Unexpected State %d\n", eBufState);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 3) {
        printf("5_2. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }    

   if (outData[0] != 2)
   {
        printf("Unexpected data returned. get[2]\n");
        return false;
   }

    // resize the buffer to 10.
    Init ( 64, 10, true, true );

   if (iBlockPutPos != 3) {
        printf("5R. Unexpected iBlockPutPos after resize %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("5R. Unexpected iBlockGetPos after resize  %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_FILLING) {
        printf("5R. Unexpected State after resize  %d\n", eBufState);
        return false;
    }

    testData[0] = 8;
    testData[64] = 8;

    if (!Put ( testData, 65)) {
        printf("8. Unexpected error put[8]\n");
        return false;
    }

    if (iBlockPutPos != 6) {
        printf("8. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("8. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 6) {
        printf("8. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }

    if (eBufState != BS_OK) {
        printf("8. Unexpected State %d\n", eBufState);
        return false;
    }

    testData[0] = 9;
    testData[64] = 9;

    if (!Put ( testData, 65)) {
        printf("9. Unexpected error put[9]\n");
        return false;
    }

    if (iBlockPutPos != 7) {
        printf("9. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("9. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 7) {
        printf("9. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }

    if (eBufState != BS_OK) {
        printf("9. Unexpected State %d\n", eBufState);
        return false;
    }

    testData[0] = 7;
    testData[64] = 7;

    if (!Put ( testData, 65)) {
        printf("7. Unexpected error put[7]\n");
        return false;
    }

   // 0  1  2  3  4  5  6  7  8  9 
   // 3, 4, 5, ?, 7, 8, 9, x, x, x
   // ^                    ^

    if (iBlockPutPos != 7) {
        printf("7. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("7. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 7) {
        printf("7. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }

    if (eBufState != BS_OK) {
        printf("7. Unexpected State %d\n", eBufState);
        return false;
    }

   if (!Get ( outData, 64)) {
        printf("Unexpected get error get[3]\n");
        return false;
   }

    if (iBlockPutPos != 7) {
        printf("8_1. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 1) {
        printf("8_1. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_OK) {
        printf("8_1. Unexpected State %d\n", eBufState);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 6) {
        printf("8_1. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }

   if (outData[0] != 3)
   {
        printf("Unexpected data returned. get[3]\n");
        return false;
   }        

   // 0  1  2  3  4  5  6  7  8  9 
   // x, 4, 5, ?, 7, 8, 9, x, x, x
   //    ^                 
   //                      ^

    // resize the buffer to 4.
    Init ( 64, 4, true, true );

   // 0  1  2  3
   // ?, 7, 8, 9
   // ^        
   // ^

    if (iBlockPutPos != 0) {
        printf("9R. Unexpected iBlockPutPos after resize %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("9R. Unexpected iBlockGetPos after resize  %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_FULL) {
        printf("9R. Unexpected State after resize  %d\n", eBufState);
        return false;
    }

   if (Get ( outData, 64)) {
        printf("Unexpected get success get[6]. Got %d\n", outData[0]);
        return false;
   }

   // 0  1  2  3
   // x, 7, 8, 9
   //    ^        
   // ^

    if (iBlockPutPos != 0) {
        printf("9_x. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 1) {
        printf("9_x. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_OK) {
        printf("9_x. Unexpected State %d\n", eBufState);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 3) {
        printf("9_x. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }    

   if (!Get ( outData, 64)) {
        printf("Unexpected get error get[7].\n");
        return false;
   }

   // 0  1  2  3
   // x, x, 8, 9
   //       ^        
   // ^

    if (iBlockPutPos != 0) {
        printf("9_7. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 2) {
        printf("9_7. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_OK) {
        printf("9_7. Unexpected State %d\n", eBufState);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 2) {
        printf("9_7. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    }    

   if (outData[0] != 7)
   {
        printf("Unexpected data returned. get[7]\n");
        return false;
   }        

   testData[0] = 17;
   testData[64] = 17;

    if (!Put ( testData, 65)) {
        printf("17. Unexpected error put[7]\n");
        return false;
    }

   // 0  1  2  3
   // ?, ?, ?, 17
   // ^        
   // ^

    if (iBlockPutPos != 0) {
        printf("17_x. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("17_x. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_FULL) {
        printf("17_x. Unexpected State %d\n", eBufState);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 0) {
        printf("17_x. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    } 

   testData[0] = 18;
   testData[64] = 18;

    if (!Put ( testData, 65)) {
        printf("17. Unexpected error put[7]\n");
        return false;
    }

   // 0   1   2   3
   // 18, ?,  ?,  17
   //     ^        
   //     ^

    if (iBlockPutPos != 1) {
        printf("18_x. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 1) {
        printf("18_x. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_FULL) {
        printf("18_x. Unexpected State %d\n", eBufState);
        return false;
    }

    if (DeltaDistance(iBlockPutPos, iBlockGetPos) != 0) {
        printf("18_x. Unexpected deltaDistances %d\n", DeltaDistance(iBlockPutPos, iBlockGetPos));
    } 

  if (Get ( outData, 64)) {
        printf("Unexpected get success get[15]. Got %d\n", outData[0]);
        return false;
   }

   // 0   1   2   3
   // 18, x,  ?,  17
   //         ^        
   //     ^

    if (iBlockPutPos != 1) {
        printf("18_15. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 2) {
        printf("18_15. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_OK) {
        printf("18_15. Unexpected State %d\n", eBufState);
        return false;
    }

  if (Get ( outData, 64)) {
        printf("Unexpected get success get[16]. Got %d\n", outData[0]);
        return false;
   }

   // 0   1   2   3
   // 18, x,  x,  17
   //             ^        
   //     ^

    if (iBlockPutPos != 1) {
        printf("18_16. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 3) {
        printf("18_16. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_OK) {
        printf("18_16. Unexpected State %d\n", eBufState);
        return false;
    }

    if (!Get ( outData, 64)) {
        printf("Unexpected get error get[17].\n");
        return false;
    }

   // 0   1   2   3
   // 18, x,  x,  x
   // ^        
   //     ^

    if (iBlockPutPos != 1) {
        printf("18_17. Unexpected iBlockPutPos %d\n", iBlockPutPos);
        return false;
    }
    if (iBlockGetPos != 0) {
        printf("18_17. Unexpected iBlockGetPos %d\n", iBlockGetPos);
        return false;
    }

    if (eBufState != BS_OK) {
        printf("18_17. Unexpected State %d\n", eBufState);
        return false;
    }

   if (outData[0] != 17)
   {
        printf("Unexpected data returned. get[17]\n");
        return false;
   } 

    return true;
}