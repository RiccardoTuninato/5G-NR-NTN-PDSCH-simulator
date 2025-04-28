function [chan, chDelay, ntnTDLParams] = chanInit(env)

% Fading Channel
   
addpath("HelperFunctionsNM")
ntnTDLParams = NaN;
chDelay = 0;

%env = loadParameters();
%[decoderDLSCH, encoderDLSCH, carrier, pdschNR, harqEntity] = initialiseObjects(env);

if env.FadingChanType == "LMS"
    % For ITU-R P.681 LMS channel
    chan = p681LMSChannel;
    % Environment type
    chan.Environment = env.environment; %env.environment; %"Urban";
    % Carrier frequency (in Hz)
    chan.CarrierFrequency = env.fc; %3.8e9;
    % Elevation angle with respect to ground plane (in degrees)
    chan.ElevationAngle = env.elevationAngle; % 45;
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
    L_corr = 7.9;
    chan.TransitionLengthCoefficients = [0; 0];
    chan.DirectPathCorrelationDistance = [L_corr, L_corr];
    chan.StateProbabilityRange = [1e-5, 1e-5; 1-1e-5, 1-1e-5];
    chan.InitialState = "Good";
     
elseif env.FadingChanType == "Lutz"
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

    chan.InitialState = "Good";

elseif env.FadingChanType == "NTN_TDL"

    chanParams = struct;
    chanParams.CarrierFrequency = env.fc;             % In Hz
    chanParams.ElevationAngle = env.elevationAngle;   % In degrees
    chanParams.SatelliteAltitude = env.satDistance;   % In m
    chanParams.SatelliteSpeed = 0;               % In m/s
    chanParams.MobileSpeed = env.terminalSpeed * (1000/3600);           % In m/s
    chanParams.SampleRate = env.PSCH.sampRate;        % In Hz

    chanParams.NumTransmitAntennas = 1;
    chanParams.NumReceiveAntennas = 1;

    % Initialize the NTN TDL channel parameters in a structure
    chanParams.NTNChannelType = "TDL";
    chanParams.DelayProfile = "NTN-TDL-D";
    chanParams.DelaySpread = env.delaySpread;                          % In s
    chanParams.TransmissionDirection = "Downlink";
    %ntnTDLParams.MIMOCorrelation = "Low";
    %ntnTDLParams.Polarization = "Co-Polar";
    % Modify the below parameters, when DelayProfile is set to Custom
    %chanParams.PathDelays = 0;                               % In s
    %chanParams.AveragePathGains = 0;                         % In dB
    %chanParams.FadingDistribution = "Rayleigh";    

    chan = nrTDLChannel;
    chan = hSetupNTNChannel(chanParams);

    if (env.SetKFactor_manually)
        K_factor_model = chan.BaseChannel.AveragePathGains(1) - pow2db(sum(db2pow(chan.BaseChannel.AveragePathGains(2:end))));
        chan.BaseChannel.AveragePathGains(2:end) = chan.BaseChannel.AveragePathGains(2:end)-env.KFactor + K_factor_model;
        chan.BaseChannel.KFactorFirstTap = chan.BaseChannel.AveragePathGains(1) - chan.BaseChannel.AveragePathGains(2);

        %%% normalize the cumulative power of the channel paths
        PathGains_lin = db2pow(chan.BaseChannel.AveragePathGains);
        totPower = sum(PathGains_lin);
        chan.BaseChannel.AveragePathGains = chan.BaseChannel.AveragePathGains - pow2db(totPower);
    else
        chan.BaseChannel.KFactorFirstTap = env.KFactor;
    end

    if env.delaySpread > 0
        %%% renormalize the delays
        mean_delay = sum(chan.BaseChannel.PathDelays.*db2pow(chan.BaseChannel.AveragePathGains))/sum(db2pow(chan.BaseChannel.AveragePathGains));
        DS_rms = sqrt(sum((abs(chan.BaseChannel.PathDelays-mean_delay).^2).*db2pow(chan.BaseChannel.AveragePathGains))/sum(db2pow(chan.BaseChannel.AveragePathGains)));
        DS_norm = DS_rms/env.delaySpread;
        chan.BaseChannel.PathDelays = chan.BaseChannel.PathDelays/DS_norm;
    end

    %%%% denormalized the delays (in sec)
    %chan.BaseChannel.PathDelays = chan.BaseChannel.PathDelays*env.delaySpread;

    %chan.BaseChannel.AveragePathGains = [-3.2355e-04 -40 -60]
    %chan.InitialTime = 100;

    tdlChanInfo = info(chan.BaseChannel);
    maxChDelay = tdlChanInfo.MaximumChannelDelay;
    chDelay = tdlChanInfo.ChannelFilterDelay;

    %chanInfo = info(chan.ChannelFilter);
    %pathFilters = chanInfo.ChannelFilterCoefficients.';
    %pathFilters = pathFilters.*10.^(chan.BaseChannel.AveragePathGains/10);

else

    RicianChannel = comm.RicianChannel( ...
    SampleRate=env.fc, ...
    KFactor=env.KFactor, ...
    MaximumDopplerShift=env.MaximumDopplerShift);
        
    chan = RicianChannel;
end

chan.RandomStream = "mt19937ar with seed";
chan.FadingTechnique = "Filtered Gaussian noise"; 

% Sampling rate (in Hz)
chan.SampleRate = env.PSCH.sampRate; %400;
chan.ChannelFiltering = true;
    
% Set random number generator with seed
chan.Seed = env.chanSeed; rng(env.chanSeed); %+73; chan.Seed = seed; rng(seed);

end