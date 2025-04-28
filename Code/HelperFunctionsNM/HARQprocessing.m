% HARQ and RLC ARQ processing and Transport Block generation

for idxCodeword = 1:pdschNR.NumCodewords
    % Check if the HARQ Entity Requires New Data (Sequence Timeout or Successful Reception)
    % When harqEntity.NewData(idxCodeword) is set to 1, the rv and TransmissionNumber properties are reset.
    if (harqEntity.NewData(idxCodeword))
        % Create the transport block
        if env.PSCH.HARQmode == "Active"
            transportBlock = randi([0 1], transportBlockSizes(idxCodeword), 1);
            % Set block (needs to be in for loop because it must be set before encoding)
            setTransportBlock(encoderDLSCH, transportBlock, idxCodeword-1, harqEntity.HARQProcessID);
        else
            if (RLC_entity(harqEntity.HARQProcessID + 1) < env.RLCARQ_MaxTx) && (RLC_entity(harqEntity.HARQProcessID + 1) ~= 0)
                RLC_entity(harqEntity.HARQProcessID + 1) = RLC_entity(harqEntity.HARQProcessID + 1) + 1;
            else
                transportBlock = randi([0 1], transportBlockSizes(idxCodeword), 1);
                % Set block
                setTransportBlock(encoderDLSCH, transportBlock, idxCodeword-1, harqEntity.HARQProcessID);
                RLC_entity(harqEntity.HARQProcessID + 1) = 1;
            end

            RLC_txNum(RLC_entity(harqEntity.HARQProcessID + 1)) = RLC_txNum(RLC_entity(harqEntity.HARQProcessID + 1)) + 1;
        end
        % If reception failed due to sequence timeout, flush the buffer
        if harqEntity.SequenceTimeout(idxCodeword)
            resetSoftBuffer(decoderDLSCH, idxCodeword-1, harqEntity.HARQProcessID);
        end
    end
end
