function [T_proc] = evaluateUEprocTime(env)

mu = log2(env.PSCH.SCS/15e3); 
%Ts = 1/env.PSCH.SCS;

% Section 5.3 TS 38.214
% N1 from Table 5.3-1
if env.UEcapability == 1
    switch mu
        case 0
            N1 = 8;
        case 1
            N1 = 10;
        case 2
            N1 = 17;
        case 3
            N1 = 20;
        otherwise
            N1 = 20; % Not defined in the table
    end
else % N1 from Table 5.3-2
    switch mu
        case 0
            N1 = 3;
        case 1
            N1 = 4.5;
        case 2
            N1 = 9;
        otherwise
            N1 = 9; % Not defined in the table
    end
end

%For operation with shared spectrum channel access, ext T is calculated 
% according to [4, TS 38.211], otherwise ext T = 0
T_ext = 0;

Deltaf_max =  480*1e3; % 480 kHz
N_f = 4096;
Tc = 1/(Deltaf_max*N_f); % Basic time unit for NR; see clause 4.1: Tc = 1/(Deltaf_max * N_f)
k = 64; % k = Ts/Tc Ts = 1/(Deltaf_ref * N_f,ref); Deltaf_ref = 15*1e3 Hz   N_f,ref = 2048

% For the PDSCH mapping type A as given in clause 7.4.1.1 of [4, TS 38.211]: 
% if the last symbol of PDSCH is on the i-th symbol of the slot where i < 7, 
% then d1,1 = 7 - i, otherwise d1,1 = 0 
d_11 = 0;
% NOTE: other possibilities for d_11, but not interesting to us now

% If a PUCCH of a larger priority index would overlap with PUCCH/PUSCH of a 
% smaller priority index, d2 for the PUCCH of a larger priority is set as 
% reported by the UE; otherwise d2 = 0
d_2 = 0;

% next uplink symbol
T_proc = (N1 + d_11 + d_2)*(2048 + 144)*k*(2^(-1*mu))*Tc + T_ext;

end