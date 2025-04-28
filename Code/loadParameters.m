%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% System parameters %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function env = loadParameters()

% env is a structure variable containing all the parameters

% Constants
env.c = 3e8; % [m/s] speed of light
env.RE = 6371; % [km] distance between Geocenter and Satellite

% Simulation parameters
env.numFrames = 6; % Minimum number of 10 ms frames to simulate
env.numMaxFrames = 6; % Maximum number of 10 ms frames to simulate
env.errorThreshold = 50; % Min of X errors/SNR

env.EsNodBRange = [17.8 18 18.2 18.3]; 

% 5G NTN System parameters
env.freqBand = "Ka-Band"; % S-Band, L-Band, Ka-Band
env.transmissionType = "DL"; %[UL, DL]
if env.transmissionType == "UL"
    % Note: that if transform precoding is enabled, the number of layers should be set to 1.
    env.PSCH.ULPrecodingFlag = true; 
end
switch env.freqBand  
    case "L-Band"
        if env.transmissionType == "UL"
            env.fc = 1643.5e6;
        elseif env.transmissionType == "DL"
            env.fc = 1543.5e6;
        end
    case "S-Band"
        if env.transmissionType == "UL"
            env.fc = 2000e6;
        elseif env.transmissionType == "DL"
            env.fc = 2185e6;
        end
     case "Ka-Band"
        if env.transmissionType == "UL"
            env.fc = 30e9;
        elseif env.transmissionType == "DL"
            env.fc = 20e9;
        end
end
env.PSCH.SCS = 120e3; %Subcarrieri spacing in [Hz], from [15, 30, 60, 120, 240, 480, 960] kHz - Table 4.2-1 of 38.211
env.PSCH.PRBSize = 140;
env.PSCH.NPRB = 140;
env.PSCH.chBandwidth = env.PSCH.PRBSize * 12 * env.PSCH.SCS;
env.PSCH.NFFT = max(128, 2^nextpow2(ceil(env.PSCH.PRBSize*12 / 0.85))); % Nfft results in a maximum occupancy of 85% and Nfft minimum is 128
env.PSCH.sampRate = env.PSCH.NFFT*env.PSCH.SCS; 
env.PSCH.CPtype = "Normal"; % Cyclic Prefix type, "Normal" or "Extended"
if env.PSCH.CPtype == "Normal"
    env.PSCH.N_OFDM_symb = 14;
else
    env.PSCH.N_OFDM_symb = 12;
end

% Sat parameters
env.elevationAngle = 90; % [degree]
env.orbit = "LEO";
env.satDistance = 1200; % [km]
env.slantRange = sqrt((env.RE*sin(deg2rad(env.elevationAngle)))^2 + env.satDistance^2 + 2*env.satDistance*env.RE) - env.RE*sin(deg2rad(env.elevationAngle)); % [km]
env.RTT = 2*env.slantRange*1e3/env.c; % RTT in [s] = twice the propagation time
env.beamConf = "fixed";

% UE and scenario parameters
env.deviceType = "Handheld";
env.terminalSpeed = 50; % [km/h]
env.environment = 'custom'; % [vehicular urban suburban rural rural wooded]
env.delaySpread = 20e-9; % [s]
env.FreeSpacePathLoss = 10*log10((4*pi*env.slantRange*1e3*env.fc/env.c)^2);

% Transmitter parameters
env.PSCH.NHARQProcesses = 16; % max(16, ceil(env.RTT/(env.PSCH.N_OFDM_symb*1/env.PSCH.SCS))); % Max. 16 HARQ processes active, maximum from standard
% 5G Toolbox accept a maximum of 32 HARQ processes. To change this, modify the following files:
% toolbox\5g\5g\+nr5g\+internal\validateParameters -> line 32: {'scalar','integer','>=',0,'<=',N},fcnName,'HARQID');
% toolbox\5g\5g\nrDLSCHDecoder -> line 220: MaxNumHARQProcesses = N;
env.PSCH.HARQmode = "Disable"; %"Disable"; 
env.PSCH.HARQrvSequence = [0 2 3 1]; % Redundancy version sequence for HARQ. At most 4 transmissions, i.e., 3 retransmissions
env.PSCH.TargetCodeRate = 948/1024; % Look at Table 5.1.3.1-(1,2,3,4) TS 38 214
env.PSCH.LabelCodeRate = '948/1024'; % Rate label
env.PSCH.modulation = "64QAM"; %""64QAM", "16QAM"; %"QPSK";

% RLC ARQ
if env.PSCH.HARQmode == "Disable"
    env.RLCARQ = "Disable"; %"Disable"; %"Active";
    env.PSCH.HARQrvSequence = [0]; % Only RV 0 is transmitted if HARQ is disabled (self-decodable)
    env.RLCARQ_MaxTx = 1; % Maximum number of RLC ARQ transmissions, i.e. retransmissions + 1
else
    env.PSCH.RLCARQ = "Disable"; %"Disable"; %"Active"
    env.RLCARQ_MaxTx = 1; % Maximum number of RLC ARQ transmissions, i.e. retransmissions + 1
end

% Channel parameters
env.channelModel = "AWGN"; %"Fading"; %"AWGN";
if strcmp(env.channelModel, "Fading")
    env.KFactor = 15; % Rician K factor of the channel
    env.MaximumDopplerShift = env.fc * (env.terminalSpeed * (1000/3600) / env.c);    
    env.FadingChanType = "LMS"; % (LMS, Lutz, NTN_TDL, Rician)
    env.ResetChannel = 1; % 1 reset (release)| 0 time evolution with fixed seed
    env.chanSeed = 73;
    if env.FadingChanType == "NTN_TDL" 
        env.SetKFactor_manually = true; %% flag to set manually the K-factor
    end
else
    env.FadingChanType = "";
end

env.PSCH.NumLayers = 1; % Codeword transmitted in a single layer (single antenna)
env.PSCH.VRBToPRBInterleaving = 0;
env.PSCH.EnablePTRS = 0;
env.timeDensityPTRS = 2;
env.PSCH.Overhead = 0; % Extra overhead on PDPSCH

% Receiver parameters
env.LDPCDecodingAlgorithm = "Layered belief propagation"; %"Layered belief propagation"; "Normalized min-sum"
env.MaximumLDPCIterationCount = 25;
env.eqType = "MMSEIRC"; % "MMSEIRC" "none"
env.perfectChanEstim = true;

end

