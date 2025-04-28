function [txBlocksUnique, txBitsUnique, rxInfoBits, rxBlocks, rvsTransmitted, ...
    rvDistribution, RLC_entity, blockErrorsPostHARQ, blockErrorsPostRLCARQ, blockErrors] = metricsUpdate(blockError, txBlocksUnique, txBitsUnique, rxInfoBits, rxBlocks, ...
             rvsTransmitted, transportBlockSizes, harqEntity, env, ...
            rvDistribution, blockErrorsPostHARQ, blockErrorsPostRLCARQ, rvSequenceLength, idxEsNo, RLC_entity, idxCodeword, blockErrors)

if (blockError == 0)

    % Successfully received a unique block transmitted - HARQ
    % Process Finished
    txBlocksUnique(idxEsNo) = txBlocksUnique(idxEsNo) + 1;
    txBitsUnique(idxEsNo) = txBitsUnique(idxEsNo) + transportBlockSizes;

    % Count only successfully received blocks and bits
    rxInfoBits(idxEsNo) = rxInfoBits(idxEsNo) + transportBlockSizes;
    rxBlocks(idxEsNo) = rxBlocks(idxEsNo) + 1;

    % Sucessfully received block, after a total number of
    % retransmissions
    rvsTransmitted(idxEsNo) = rvsTransmitted(idxEsNo) + (harqEntity.TransmissionNumber + 1);

    % Update rv distribution, unrolled for parfor
    if env.PSCH.HARQmode == "Active"
        rvDistribution(1:harqEntity.TransmissionNumber+1,idxEsNo) = rvDistribution(1:harqEntity.TransmissionNumber+1,idxEsNo) + 1;
    else
        RLC_entity(harqEntity.HARQProcessID + 1) = 0;
    end

elseif (blockError == 1)
    blockErrors = blockErrors + 1;
    if env.PSCH.HARQmode == "Active"
        % Block error after the last r.v.
        if (harqEntity.TransmissionNumber == (rvSequenceLength-1))
            % When the last rv has been sent and received in error,
            % increase total number of block errors.

            % Failed to receive a unique block that was transmitted.
            % HARQ Processes that have not finished are not counted
            % below.
            txBlocksUnique(idxEsNo) = txBlocksUnique(idxEsNo) + 1;
            if env.transmissionType == "DL"
                txBitsUnique(idxEsNo) = txBitsUnique(idxEsNo) + transportBlockSizes(idxCodeword);
            elseif env.transmissionType == "UL"
                txBitsUnique(idxEsNo) = txBitsUnique(idxEsNo) + transportBlockSizes;
            end

            blockErrorsPostHARQ(idxEsNo) = blockErrorsPostHARQ(idxEsNo) + 1;
            rvsTransmitted(idxEsNo) = rvsTransmitted(idxEsNo) + rvSequenceLength; % rvSequenceLength: 4 rvs by default

            % Update rv distribution, unrolled for parfor
            rvDistribution(1:rvSequenceLength,idxEsNo) = rvDistribution(1:rvSequenceLength,idxEsNo) + 1;
        end
    else
        if RLC_entity(harqEntity.HARQProcessID + 1) == env.RLCARQ_MaxTx
            
            blockErrorsPostRLCARQ(idxEsNo) = blockErrorsPostRLCARQ(idxEsNo) +  1;
            
            txBlocksUnique(idxEsNo) = txBlocksUnique(idxEsNo) + 1;
            if env.transmissionType == "DL"
                txBitsUnique(idxEsNo) = txBitsUnique(idxEsNo) + transportBlockSizes(idxCodeword);
            elseif env.transmissionType == "UL"
                txBitsUnique(idxEsNo) = txBitsUnique(idxEsNo) + transportBlockSizes;
            end
        end

        blockErrorsPostHARQ(idxEsNo) = blockErrorsPostHARQ(idxEsNo) + 1;
        rvsTransmitted(idxEsNo) = rvsTransmitted(idxEsNo) + rvSequenceLength;
    end
end

end