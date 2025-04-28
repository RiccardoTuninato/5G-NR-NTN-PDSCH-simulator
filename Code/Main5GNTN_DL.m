%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% NR NTN PDSCH Throughput %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 5G NTN Link-level PDSCH simulator (only DL) - Ver 1.0
%
% - ITU LMS or NTN TDL channel models supported
% - PDSCH, DMRS and PTRS symbols
% - BLER and SE evaluation
% - Parameters can be set in loadParameters.m
% 
% Mainly based on MATLAB 5G toolbox and Satellite communications toolbox
% MATLAB example: https://nl.mathworks.com/help/satcom/ug/nr-ntn-pdsch-throughput.html 
%
% Software developed by Politecnico di Torino for TEC-ESC during the
% activity "Physical Layer Simulator for 5G New Radio and OTFS".
% ESA Contract No. 4000141791/23/NL/KK/nh
%
% For questions contact: riccardo.tuninato@ext.esa.int
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

addpath("HelperFunctionsNM")
addpath("ChannelGeneration")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     SETUP     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Transmission parameters
env = loadParameters();

%% Initialise Objects
[decoderDLSCH, encoderDLSCH, carrier, pdschNR, harqEntity, pdschNRExtension] = initialiseObjects_DL(env);

%% System and simulation metrics initialization

nPoints = length(env.EsNodBRange);
rxInfoBits = zeros(1, nPoints);             % Total number of received (information) bits
rxBlocks = zeros(1, nPoints);               % Total number of received blocks
txBitsTransmitted = zeros(1, nPoints);      % Total number of transmitted bits
txBitsUnique = zeros(1, nPoints);           % Transmitted information bits, no repetitions
txBlocksTransmitted = zeros(1, nPoints);    % Total number of transmitted blocks, including repetitions
txBlocksUnique = zeros(1, nPoints);         % Total number of unique blocks transmitted, no repetitions
blockErrors = zeros(1, nPoints);            % Total number of block errors
blockErrorsPostHARQ = zeros(1, nPoints);    % Total number of block errors after HARQ retxs
rvsTransmitted = zeros(1, nPoints);         % Counter for redundancy version transmitted

rvSequenceLength = length(env.PSCH.HARQrvSequence); 
rvDistribution = zeros(rvSequenceLength, nPoints);
% Example 1: if rvSequence = [0 2 3 1], rv0 = rv(1); rv2 = rv(2); rv3 = rv(3); rv1 = rv(4);
% Example 2: if rvSequence = [0 3 0 3], rv0_1 = rv(1); rv3_1 = rv(2); rv0_2 = rv(3); rv3_2 = rv(4);

RLC_entity = zeros(env.PSCH.NHARQProcesses, nPoints);
RLC_txNum = zeros(env.RLCARQ_MaxTx, nPoints);
blockErrorsPostRLCARQ = zeros(1, nPoints);

% Total number of slots in simulation period
nSlots = env.numFrames * carrier.SlotsPerFrame;
nMaxSlots = env.numMaxFrames * carrier.SlotsPerFrame;

if env.channelModel == "Fading"
    [chan_init, chDelay] = chanInit(env);
end

%%  Disp output that summarize main parameters 
disp("Simulation of "+env.transmissionType+" at "+env.fc/1e9+" GHz. "+env.PSCH.SCS/1e3+" kHz SCS, "+env.PSCH.NPRB+" PRB ")
disp("Number of simulated 10 ms frames: "+env.numFrames+", resulting in "+10*env.numFrames/(15e3/env.PSCH.SCS)+" slots.")
disp("LDPC algo "+env.LDPCDecodingAlgorithm+" with "+env.MaximumLDPCIterationCount+" iteration(s).")
if env.channelModel == "AWGN"
    disp("AWGN channel and "+env.eqType+" equalization.")
else
    disp("Fading channel, "+env.FadingChanType+" model, K factor "+env.KFactor+" and user speed "+env.terminalSpeed+" [km/h] and "+env.eqType+" equalization.")
end

%% Start simulation cycle with different EsNo values

for idxEsNo = 1:length(env.EsNodBRange)

    blockError = 0;

    EsNodB = env.EsNodBRange(idxEsNo);
    disp("Iteration number "+idxEsNo+" for EsNo "+EsNodB+" dB")

    decoderDLSCH.reset();
    encoderDLSCH.reset();
    harqEntity = HARQEntity(0:(env.PSCH.NHARQProcesses-1), env.PSCH.HARQrvSequence, pdschNR.NumCodewords);
    
    waveformInfo = nrOFDMInfo(carrier);

    RLC_entity = zeros(env.PSCH.NHARQProcesses, 1);
    RLC_txNum = zeros(env.RLCARQ_MaxTx, 1);

    if env.channelModel == "Fading"
        chan = chan_init;
    end

    % While min. block error threshold or simulation period not over
    stopCondition = 0; currentSlot = 0;
    while ~stopCondition

        % Update carrier slot
        carrier.NSlot = currentSlot;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%  TRANSMITTER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Generate PDSCH channel information
        [pdschIndices, pdschIndicesInfo] = nrPDSCHIndices(carrier, pdschNR);

        % Get transport block size -  TS 38.214
        transportBlockSizes = nrTBS(pdschNR.Modulation, pdschNR.NumLayers, numel(pdschNR.PRBSet), pdschIndicesInfo.NREPerPRB, env.PSCH.TargetCodeRate, env.PSCH.Overhead);

        % HARQ and RLC ARQ processing and Transport Block generation
        HARQprocessing

        % Encode DLSCH blocks
        codedTransportBlock = encoderDLSCH(pdschNR.Modulation, pdschNR.NumLayers, ...
            pdschIndicesInfo.G, harqEntity.RedundancyVersion, harqEntity.HARQProcessID);

        % Generate PDSCH Symbols (Complex)
        pdschSymbols = nrPDSCH(carrier, pdschNR, codedTransportBlock);
        dmrsIndices = nrPDSCHDMRSIndices(carrier, pdschNR);
        dmrsSymbols = nrPDSCHDMRS(carrier, pdschNR);

        % Perform implementation-specific PDSCH PT-RS mapping
        if env.PSCH.EnablePTRS == true
            ptrsIndices = nrPDSCHPTRSIndices(carrier, pdschNR);
            ptrsSymbols = nrPDSCHPTRS(carrier, pdschNR);     
        end

        % Generate resource grid
        pdschGrid = nrResourceGrid(carrier);

        % Set symbols in resource grid
        pdschGrid(pdschIndices) = pdschSymbols;
        pdschGrid(dmrsIndices) = dmrsSymbols;
        
        if env.PSCH.EnablePTRS == true
            pdschGrid(ptrsIndices) = ptrsSymbols;
        end
        pschGrid = pdschGrid;

        % Time domain
        [txWaveform, infoWf] = nrOFDMModulate(carrier, pschGrid);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%  CHANNEL  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Fading channel
        if env.channelModel == "Fading"
            [chOutWaveform, ChanCoeff, sampleTimes, chan] = fadingChannel(env, txWaveform, currentSlot, chan);
        else
            chOutWaveform = txWaveform;
        end

        offset = 0;
        if env.perfectChanEstim
            pathGains = 1; pathFilters = 1;
            if env.channelModel == "Fading"
                pathGains = ChanCoeff; % pathFilters = 1; offset = 0; % LMS
                if env.FadingChanType == "NTN_TDL"
                    chanInfo = info(chan.ChannelFilter);
                    pathFilters = chanInfo.ChannelFilterCoefficients.';
                    [offset, mag] = nrPerfectTimingEstimate(pathGains, pathFilters);
                end
            end
        else
            % Assume perfect timing synchronization with TDL channel
            % (time sync is not relevant in this version)
            if env.channelModel == "Fading" && env.FadingChanType == "NTN_TDL"
                pathGains = ChanCoeff; % pathFilters = 1; offset = 0; % LMS
                chanInfo = info(chan.ChannelFilter);
                pathFilters = chanInfo.ChannelFilterCoefficients.';
                [offset, mag] = nrPerfectTimingEstimate(pathGains, pathFilters);
            end
        end
        
        % Complex AWGN channel
        EsNo = 10^(EsNodB/10);
        N0 = 1/sqrt(2.0*double(waveformInfo.Nfft)*EsNo);
        noise = N0 * complex(randn(size(chOutWaveform)),randn(size(chOutWaveform)));
        rxWaveform = chOutWaveform + noise;
        if env.channelModel == "Fading"
            % Time-sample alignment
            rxWaveform = circshift(rxWaveform,-1*offset,1); % rxWaveform(1+offset:end-offset+1,:);
            rxWaveform = rxWaveform(1:size(txWaveform,1),:);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%  RECEIVER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Demodulate noisy grid
        rxGrid = nrOFDMDemodulate(carrier, rxWaveform);

        if env.perfectChanEstim
            % Perfect estimation
            noiseGrid = nrOFDMDemodulate(carrier, noise);
            noiseEst = var(noiseGrid(:));
            %estChGrid = (rxGrid - noiseGrid)./pdschGrid;
            estChGrid = nrPerfectChannelEstimate(carrier, pathGains, pathFilters);
        else
            % Practical channel estimation between the received grid and
            % each transmission layer, using the PSCH DM-RS for each layer
            [estChGrid,noiseEst] = hSubbandChannelEstimate(carrier, rxGrid, ...
                dmrsIndices, dmrsSymbols, pdschNRExtension.PRGBundleSize, 'CDMLengths', pdschNR.DMRS.CDMLengths);

            % Average noise estimate across PRGs and layers
            noiseEst = mean(noiseEst,'all');
        end

        [pdschRx,pdschEst] = nrExtractResources(pdschIndices, rxGrid, estChGrid);
        
        % Equalization
        if env.eqType == "MMSEIRC"
            [pdschEq,csi] = nrEqualizeMMSEIRC(pdschRx,pdschEst,noiseEst,1,0);
        else
            pdschEq = pdschRx;
        end

        % Demodulate PDSCH Symbols
        [rxSoftBits, rxSymbols] = nrPDSCHDecode(pdschEq, pdschNR.Modulation, pdschNR.NID, pdschNR.RNTI, noiseEst); % LLR scaled by noise variance/power

        % Scale LLRs by CSI
        if env.eqType == "MMSEIRC"
            csi = nrLayerDemap(csi); % CSI layer demapping
            Qm= length(rxSoftBits{1})/length(rxSymbols{1}); % bits per symbol
            csi{1} = repmat(csi{1}.',Qm,1);   % expand by each bit per symbol
            rxSoftBits{1} = rxSoftBits{1} .* csi{1}(:);   % scale
        end

        % Decode DLSCH blocks
        decoderDLSCH.TransportBlockLength = transportBlockSizes;
        [decodedBits, blockError] = decoderDLSCH(rxSoftBits, pdschNR.Modulation, ...
            pdschNR.NumLayers, harqEntity.RedundancyVersion, harqEntity.HARQProcessID);

        % Update and advance current HARQ process
        updateAndAdvance(harqEntity,blockError,transportBlockSizes,pdschIndicesInfo.G);

        % Store metrics
        % Total blocks/bits transmitted overall
        txBitsTransmitted(idxEsNo) = txBitsTransmitted(idxEsNo) + transportBlockSizes; % Total number of transmitted bits
        txBlocksTransmitted(idxEsNo) = txBlocksTransmitted(idxEsNo) + 1;
        
        % Store and update performance metrics
        [txBlocksUnique, txBitsUnique, rxInfoBits, rxBlocks, rvsTransmitted, rvDistribution, RLC_entity, blockErrorsPostHARQ, blockErrorsPostRLCARQ, blockErrors(idxEsNo)] = metricsUpdate(blockError, ...
            txBlocksUnique, txBitsUnique, rxInfoBits, rxBlocks, rvsTransmitted, transportBlockSizes, harqEntity, env, rvDistribution, blockErrorsPostHARQ, blockErrorsPostRLCARQ, ...
            rvSequenceLength, idxEsNo, RLC_entity, idxCodeword, blockErrors(idxEsNo));


        % Check stopping condition
        condition1 = currentSlot < nMaxSlots;
        condition2 = (blockErrorsPostHARQ(idxEsNo) < env.errorThreshold) || (currentSlot < nSlots);
        if ~(condition1 && condition2)
            stopCondition = 1;
        end

        % Change slot
        currentSlot = currentSlot + 1;
    end

    RLC_entity(:, idxEsNo) = RLC_entity;
    RLC_txNum(:, idxEsNo) = RLC_txNum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RESULTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute final metrics and clean variables
nSlots = txBlocksTransmitted; 
bitsEfficiency = rxInfoBits./(carrier.SymbolsPerSlot * nSlots);
spectralEfficiency = bitsEfficiency/env.PSCH.chBandwidth/(1/env.PSCH.SCS);
effectiveChannelRate = rxInfoBits./txBitsTransmitted;

%% Plotting

figure;
subplot(2,2,1);
if env.PSCH.HARQmode == "Disable"
    semilogy(env.EsNodBRange, blockErrorsPostRLCARQ./txBlocksUnique, 'o-','linewidth',2,'Markersize',8);
else
    semilogy(env.EsNodBRange, blockErrorsPostHARQ./txBlocksUnique, 'o-','linewidth',2,'Markersize',6);
end
xlabel('E_{s}/N_{0} [dB]');
xlim([min(env.EsNodBRange) max(env.EsNodBRange)])
ylabel('BLER');
grid on;

subplot(2,2,2);
if env.PSCH.HARQmode == "Disable"
    plot(env.EsNodBRange, sum(RLC_txNum,1)./(txBlocksUnique),'o-','linewidth',2,'Markersize',6);
    ylabel('# Transmitted TB');
else
    plot(env.EsNodBRange, rvsTransmitted./(txBlocksUnique),'o-','linewidth',2,'Markersize',6);
    ylabel('# Transmitted RVs');
end
xlabel('E_{s}/N_{0} [dB]');
grid on;

subplot(2,2,3);
plot(env.EsNodBRange, spectralEfficiency, 'o-','linewidth',2,'Markersize',6)
xlabel('E_{s}/N_{0} [dB]');
ylabel('Spectral Efficiency [info bits/Hz]');
grid on;

subplot(2,2,4);
plot(env.EsNodBRange, effectiveChannelRate, 'o-','linewidth',2,'Markersize',6)
ylabel('Effective Code Rate [including retransmissions]');
xlabel('E_{s}/N_{0} [dB]');
grid on;

titleStr = {'SCS =' + string(num2str(env.PSCH.SCS/1e3)) + 'kHz, ' + env.channelModel + ' (' + env.PSCH.modulation + '-' + env.PSCH.LabelCodeRate + '), BW = ' + num2str(env.PSCH.chBandwidth/1e6) + ' [MHz]', ...
    'HARQ:' + env.PSCH.HARQmode +  ',[' + num2str(env.PSCH.HARQrvSequence) + ']'};
sgtitle(titleStr);

filename = "Results/";

if env.channelModel == "Fading"
    filename =  filename + "Perf_%s_%.2f_SCS%d_PRB_%d_chan_%s_Kfact_%d_Tspeed_%d_DS_%d_RLCtx_%d.mat";
    filename = sprintf(filename, env.PSCH.modulation, env.PSCH.TargetCodeRate, ...
        env.PSCH.SCS, env.PSCH.NPRB, env.FadingChanType, env.KFactor, env.terminalSpeed,  env.delaySpread, env.RLCARQ_MaxTx);
else
    filename = filename + "Perf_%s_%.2f_SCS%d_PRB_%d_eq_%s_chan_%s_RLCtx_%d.mat";
    filename = sprintf(filename, env.PSCH.modulation, env.PSCH.TargetCodeRate, ...
        env.PSCH.SCS/1e3, env.PSCH.NPRB, env.eqType, env.channelModel, env.RLCARQ_MaxTx);
end

save(filename, "env", "waveformInfo", "carrier", "pdschNR", "txBitsTransmitted","txBitsUnique", "txBlocksTransmitted","txBlocksUnique", ...
    "blockErrorsPostHARQ", "spectralEfficiency", "effectiveChannelRate", "rvsTransmitted", ...
    "rxInfoBits", "nSlots", "bitsEfficiency", "rvsTransmitted", "rxBlocks", "blockErrorsPostRLCARQ", ...
    "RLC_entity", "RLC_txNum", "rvDistribution")
