function [decoderDLSCH, encoderDLSCH, carrier, pdschNR, harqEntity, pdschNRExtension] = initialiseObjects_DL(env)

%% Create encoder/decoder
% Create encoder object for DL-PDSCH
% TS 38.212 Section 7.2
% CRC -> Code Block Segmentation + CRC -> LDPC ->
% Rate Matching -> Code Block Concatenation

    
% Create UL-SCH encoder System object to perform transport channel encoding
encoderDLSCH = nrDLSCH;

if (env.PSCH.NHARQProcesses > 1)
    encoderDLSCH.MultipleHARQProcesses = true;
else
    encoderDLSCH.MultipleHARQProcesses = false;
end

encoderDLSCH.TargetCodeRate = env.PSCH.TargetCodeRate;

% Create decoder object for DL-PDSCH
% TS 38.212 Section 7.2
% Rate Recovery -> LDPC Decoding (soft, NMS) ->
% Desegmentation -> CRC Decoding

decoderDLSCH = nrDLSCHDecoder;

if (env.PSCH.NHARQProcesses > 1)
    decoderDLSCH.MultipleHARQProcesses = true;
else
    decoderDLSCH.MultipleHARQProcesses = false;
end

decoderDLSCH.TargetCodeRate = env.PSCH.TargetCodeRate;
decoderDLSCH.LDPCDecodingAlgorithm = env.LDPCDecodingAlgorithm;
decoderDLSCH.MaximumLDPCIterationCount = env.MaximumLDPCIterationCount;

%% Configure carrier
% TS 38.212 Sections 4.2, 4.3, 4.4
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = env.PSCH.SCS/1e3;
carrier.CyclicPrefix = env.PSCH.CPtype;
carrier.NSizeGrid = env.PSCH.NPRB;

% SymbolsPerSlot, SlotsPerSubframe, SlotsPerFrame are automatically
% computed with regards to SCS and CPtype

%% Configure PDSCH Channel
% TS 38.211

pdschNR = nrPDSCHConfig();
pdschNR.RNTI = 1;
pdschNR.NID = 0; % Scrambling sequence 0
pdschNR.Modulation = env.PSCH.modulation;
pdschNR.PRBSet = (0:env.PSCH.NPRB-1); % All resource blocks are used (full band)
pdschNR.NumLayers = env.PSCH.NumLayers;
pdschNR.VRBToPRBInterleaving = env.PSCH.VRBToPRBInterleaving;
pdschNR.EnablePTRS = env.PSCH.EnablePTRS;
pdschNR.MappingType = 'A';
pdschNR.PRBSetType = 'PRB'; % To be compatible with Nicolo's assumptions. After 2023a, default is VRB to PRB

% CSI-RS Reserved 
% Two reserved PRB patterns
reservedConfig1 = nrPDSCHReservedConfig();
reservedConfig1.SymbolSet = [];
reservedConfig1.PRBSet = [];
reservedConfig1.Period = [];

reservedConfig2 = nrPDSCHReservedConfig();
reservedConfig2.PRBSet = [23 24 25 26 27 28];
reservedConfig2.SymbolSet = [];
reservedConfig2.Period = 40;

pdschNR.ReservedPRB = {reservedConfig1, reservedConfig2};

% PDSCH DM-RS configuration
pdschNR.DMRS.DMRSTypeAPosition = 2;      % First DM-RS symbol position for type A
pdschNR.DMRS.DMRSLength = 1;             % Specify double front-loaded DM-RS
pdschNR.DMRS.DMRSAdditionalPosition = 0; % Specify an additional DM-RS pair
pdschNR.DMRS.DMRSConfigurationType = 1;  % DM-RS configuration type (1,2)
pdschNR.DMRS.NumCDMGroupsWithoutData = 1;% CDM groups without data [1, 2, 3]; 1 corresponds to group zero
pdschNR.DMRS.DMRSPortSet = 0;
pdschNR.DMRS.NIDNSCID = 1;               % DM-RS scrambling identity (0...65535)
pdschNR.DMRS.NSCID = 0;                  % DM-RS scrambling initialization (0,1)

% PT-RS configuration (TS 38.211 Section 7.4.1.2)
pdschNR.EnablePTRS = env.PSCH.EnablePTRS;
pdschNR.PTRS.TimeDensity = 1;
pdschNR.PTRS.FrequencyDensity = 2;
pdschNR.PTRS.REOffset = '00';
% PT-RS antenna port, subset of DM-RS port set. Empty corresponds to lowest
% DM-RS port number
pdschNR.PTRS.PTRSPortSet = [];

% No VRB>PRB interleaving, no reference signals by default
% Scrambling using default cell identity (1)

pdschNRExtension = struct();
pdschNRExtension.PRGBundleSize = [];     % 2, 4, or [] to signify "wideband"


%% Create HARQ Entity (HARQ Management)
harqSequence = 0:(env.PSCH.NHARQProcesses-1);
harqEntity = HARQEntity(harqSequence, env.PSCH.HARQrvSequence, pdschNR.NumCodewords);

% Contains:
% - HARQ ID
% - RedundancyVersion (Current, per codeword)
% - TransmissionNumber (# of retransmissions for that block)
% - Flag that indicates new data required (successful reception or sequence
% timeout)
% - Flag that indicates sequence timeout

% Starting points of the circular buffer are such that RV0 and RV3 are
% self-decodable, soft-combining approach (store & combine rather than
% discard failed packet)
% [0,2,3,1] > every second retransmission is self decodable
% New data indicator for code-block group

% After exhausting RVs > CBG is dropped (sequence timeout)

% Stop & Wait
% Multiple processes in parallel to sustain a high throughput (some
% processes are hung waiting for new transmission or timeout)
% Transmitter informs which process the data relates to

end