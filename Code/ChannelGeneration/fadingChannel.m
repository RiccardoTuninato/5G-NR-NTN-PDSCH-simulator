function [outWave, channelCoefficients, sampleTimes, chan] = fadingChannel(env, in, nslot, chan)

% Fading Channel

if env.ResetChannel == 1
     
    if env.FadingChanType == "NTN_TDL"
        chan.BaseChannel.reset;
        chan.ChannelFilter.reset;
    else
        release(chan);  
    end
    % Set random number generator with seed
    seed = nslot; chan.Seed = seed; rng(seed); %+73; chan.Seed = seed; rng(seed);
end

if env.FadingChanType == "NTN_TDL"
    tdlChanInfo = info(chan.BaseChannel);
    chDelay = tdlChanInfo.ChannelFilterDelay;
    in = [in; zeros(chDelay,size(chDelay,2))];
    % Generate the faded waveform for NTN TDL channel
    [outWave, channelCoefficients, sampleTimes] = hGenerateNTNChannel(chan,in);
else
    % Pass the input signal through channel
    [outWave, channelCoefficients, sampleTimes,stateSeries] = step(chan, in);
end


end