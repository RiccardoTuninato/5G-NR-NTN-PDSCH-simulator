function [fadedWave] = LMS_channel_fix(env, in, ch_coeff)

% LMS Channel

% filename = "LMS/LMSChan_freqBand_%s_%s_elAngle_%d_env_%s.mat";
% filename = sprintf(filename, env.freqBand, env.transmissionType, env.elevationAngle, env.environment);
% load(filename);

if size(in) == size(ch_coeff(1:length(in)))
    fadedWave = in.*ch_coeff(1:length(in));
else
    fadedWave = in.*ch_coeff(1:length(in)).';
end

%save(filename, "ch_coeff")

end
