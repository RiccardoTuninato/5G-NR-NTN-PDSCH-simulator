

% LMS Channel
% chan = p681LMSChannel;
clc
clear

addpath("HelperFunctionsNM")

env = loadParameters();

[decoderDLSCH, encoderDLSCH, carrier, pdschNR, harqEntity] = initialiseObjects_DL(env);

Ts = 1/env.PSCH.SCS;
nSlots = env.numFrames * carrier.SlotsPerFrame;
Tot_time = nSlots * Ts * env.PSCH.N_OFDM_symb;

if env.LMSChanType == "p681LMSChannel"
    % For ITU-R P.681 LMS channel
    chan = p681LMSChannel;
    % Environment type
    chan.Environment = env.environment; %env.environment; %"Urban";
    % Carrier frequency (in Hz)
    chan.CarrierFrequency = env.fc; %3.8e9;
    % Elevation angle with respect to ground plane (in degrees)
    chan.ElevationAngle = 90; %env.elevationAngle; % 45;
    % Speed of movement of ground terminal (in m/s)
    chan.MobileSpeed = env.terminalSpeed * (1000/3600); %2;
    % Direction of movement of ground terminal (in degrees)
    chan.AzimuthOrientation = 0;
    
    % Parameters for the states properties
    % Parameters of state duration distribution in dB, specified as a 2-by-2 matrix
    chan.StateDistribution = [1e2, 1e-5; 1e2, 1e0]; % [muG, muB; sigmaG, sigmaB]
    % Minimum duration of each state in meters, specified as a two-element row vector. The first element corresponds to good state and the second element corresponds to bad state.
    chan.MinStateDuration = [1e4, 1e-1]; % [G, B]
    
    % Parameters of direct path amplitude distribution (dB)
    chan.DirectPathDistribution = [0, 0; 0, 0]; % [muMaG, muMaB; sigmaMaG, sigmaMaB] Good = Bad

    % Coefficients to compute the multipath power, specified as a 2-by-2 matrix
    chan.MultipathPowerCoefficients = [0, 0; -env.KFactor, -env.KFactor]; %[0, 0; -env.KFactor, -env.KFactor];  %  [h1G, h1B; h2G, h2B] (default [-0.0481 0.9434; -14.7450 -1.7555])
    % Coefficients to compute standard deviation of direct path amplitude in all states, specified as a 2-by-2 matrix
    chan.StandardDeviationCoefficients = [0, 0; 0, 0];  %  [g1G, g1B; g2G, g2B]
    
    % Direct path amplitude correlation distance (m)
    L_corr = 7.9; % default [1.7910 1.7910] 7.9
    chan.TransitionLengthCoefficients = [0; 0];
    chan.DirectPathCorrelationDistance = [L_corr, L_corr];
    chan.StateProbabilityRange = [1e-5, 1e-5; 1-1e-5, 1-1e-5];
    
    if env.environment == "custom"
        filename = "LMS/ITUChan_freqBand_%s_%s_elAngle_%d_Kfact_%d_Tspeed_%d.mat"; % _Nframes_%d
        filename = sprintf(filename, env.freqBand, env.transmissionType, env.elevationAngle, env.KFactor, env.terminalSpeed);
    else
        filename = "LMS/ITUChan_freqBand_%s_%s_elAngle_%d_env_%s_Tspeed_%d.mat"; % _Nframes_%d
        filename = sprintf(filename, env.freqBand, env.transmissionType, env.elevationAngle, env.environment, env.terminalSpeed); %, env.numFrames
    end
    % 
elseif env.LMSChanType == "lutzLMSChannel"
    % For Lutz LMS channel
    chan = lutzLMSChannel;
    % Rician K-factor (in dB)
    chan.KFactor = env.KFactor;
    % Lognormal fading parameters (in dB)
    chan.LogNormalFading = [0 0]; %[-13.6 3.8];
    % State duration distribution
    chan.StateDurationDistribution = "Exponential";
    % Mean state duration (in seconds)
    chan.MeanStateDuration = [1e3 0];
    % Maximum Doppler shift (in Hz)
    chan.MaximumDopplerShift = env.MaximumDopplerShift;

    filename = "LMS/LutzChan_freqBand_%s_%s_Kfact_%d_Tspeed_%d.mat"; %_Nframes_%d
    filename = sprintf(filename, env.freqBand, env.transmissionType, env.KFactor, env.terminalSpeed); %, env.numFrames
else
    ricianChan = comm.RicianChannel( ...
        SampleRate=env.fc, ...
        PathDelays=pathDelays, ...
        AveragePathGains=avgPathGains, ...
        KFactor=env.KFactor, ...
        MaximumDopplerShift=env.MaximumDopplerShift, ...
        ChannelFiltering=false, ...
        Visualization='Impulse and frequency responses');

    % For Rician channel
    chan = RicianChannel;
    % Rician K-factor (in dB)
    chan.KFactor = env.KFactor;
    % Lognormal fading parameters (in dB)
    chan.LogNormalFading = [0 0]; %[-13.6 3.8];
    % State duration distribution
    chan.StateDurationDistribution = "Exponential";
    % Mean state duration (in seconds)
    chan.MeanStateDuration = [1e3 0];
    % Maximum Doppler shift (in Hz)
    chan.MaximumDopplerShift = env.MaximumDopplerShift;

    filename = "LMS/LutzChan_freqBand_%s_%s_Kfact_%d_Tspeed_%d.mat"; %_Nframes_%d
    filename = sprintf(filename, env.freqBand, env.transmissionType, env.KFactor, env.terminalSpeed); %, env.numFrames
end



chan.RandomStream = "mt19937ar with seed";

carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = env.PSCH.SCS/1e3;
carrier.CyclicPrefix = env.PSCH.CPtype;
carrier.NSizeGrid = env.PSCH.NPRB;
waveformInfoLocal = nrOFDMInfo(carrier);

% Sampling rate (in Hz)
chan.SampleRate = env.PSCH.sampRate/env.chanSampScale; %400;

chan.InitialState = "Good";
chan.FadingTechnique = "Filtered Gaussian noise"; 

% Channel duration (in seconds)
chanDur = 1/env.PSCH.SCS*env.PSCH.N_OFDM_symb;
% Random input waveform
for sub = 0:carrier.SlotsPerSubframe-1
    numSamples_temp(sub+1) = sum(waveformInfoLocal.SymbolLengths(sub*env.PSCH.N_OFDM_symb+1:(sub+1)*env.PSCH.N_OFDM_symb)); %floor(chan.SampleRate*chanDur)+1;
end
numSamples = ceil(max(numSamples_temp)/env.chanSampScale);

% numSubframes = env.numFrames*10;
chan.NumSamples = numSamples; %sum(numSamples_temp)*numSubframes/1e2; 

chCoeff_real = int32(zeros(nSlots, numSamples));
chCoeff_imag = int32(zeros(nSlots, numSamples));

chan.ChannelFiltering = false;

%%{
for nslot = 0:nSlots-1
    
    % Set random number generator with seed
    if env.ResetChannel == 1 || nslot == 1 || chan.MobileSpeed == 0 
        seed = nslot; chan.Seed = seed; rng(seed);
    end

    numSamples_slot = ceil(numSamples_temp(mod(nslot,carrier.SlotsPerSubframe)+1)/env.chanSampScale);
    chan.NumSamples = numSamples_slot; %sum(numSamples_temp)*numSubframes/1e2;
    
    % Pass the input signal through channel
    [channelCoefficients,sampleTimes,stateSeries] = step(chan);
    chCoeff_real(nslot+1, 1:numSamples_slot) = int32(real(channelCoefficients)*1e6); %round(real(channelCoefficients),4);
    chCoeff_imag(nslot+1, 1:numSamples_slot) = int32(imag(channelCoefficients)*1e6);
    % Ch_pwr(nslot) = sum(abs(channelCoefficients).^2)/length(channelCoefficients);
    % Ch_pwr_real(nslot) = sum(abs(real(channelCoefficients)).^2)/length(channelCoefficients);
    % Ch_pwr_imag(nslot) = sum(abs(imag(channelCoefficients)).^2)/length(channelCoefficients);
    
    % if sum(stateSeries) ~= length(stateSeries)
    %     disp("change state");
    % end
    %figure; plot(abs(channelCoefficients));
    if env.ResetChannel == 1 || chan.MobileSpeed == 0
        release(chan);
    end
end
%}

%{
% Set random number generator with seed
seed = 73; chan.Seed = seed; rng(seed);


[channelCoefficients,sampleTimes,stateSeries] = step(chan);
%figure; plot(abs(channelCoefficients))
chCoeff_real = int32(real(channelCoefficients)*1e6);
chCoeff_imag = int32(imag(channelCoefficients)*1e6);
%}

%save(filename, "env", "chCoeff_real","chCoeff_imag", "chCoeff_real_sign","chCoeff_imag_sign", "sampleTimes","stateSeries","-v7.3");
save(filename, "env", "chCoeff_real","chCoeff_imag", "sampleTimes","stateSeries","-v7.3");
% figure(1)
% plot(sampleTimes,20*log10(abs(in)),sampleTimes,20*log10(abs(fadedWave)))
% title(['Power Profile of Waveform for Duration ' num2str(chanDur) ' seconds'])
% legend('Input Waveform', 'Faded Waveform')
% xlabel('Time (in s)')
% ylabel('Power (in dB)')