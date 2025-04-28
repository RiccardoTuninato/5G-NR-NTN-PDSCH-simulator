function [out,csi] = nrEqualizeMMSEIRC(rxSym,Hest,nVar,mode,irc)
%nrEqualizeMMSE MMSE Equalization
%   [OUT,CSI] = nrEqualizeMMSE(RXSYM,HEST,NVAR,MODE) returns equalized symbols
%   OUT by performing MMSE or MMSE-IRC equalization on extracted resource elements of a
%   physical channel RXSYM using the estimated channel information HEST and
%   estimated noise variance NVAR, whereas MODE selects MMSE(1) or MMSE-IRC(2). 
%   The function also returns soft channel
%   state information CSI.
% 
%   Both OUT and CSI are of same size NRE-by-P. RXSYM is of size NRE-by-R,
%   HEST is of size NRE-by-R-by-P and NVAR is a real nonnegative
%   scalar value.
%   Where,
%   NRE - Number of resource elements of a physical channel whose values
%         are extracted from each K-by-L plane of received grid, where K
%         represents the number of subcarriers and L represents the number
%         of OFDM symbols
%   P   - Number of transmission planes.
%   R   - Number of receive antennas
%
%   Example:
%   % Perform MMSE equalization on extracted resource elements of PBCH
%   % using nrEqualizeMMSE. 
%
%   % Create symbols and indices for a PBCH transmission
%   ncellid = 146;
%   v = 0;
%   E = 864;
%   cw = randi([0 1],E,1);
%   pbchTxSym = nrPBCH(cw,ncellid,v);
%   pbchInd = nrPBCHIndices(ncellid);
%
%   % Generate an empty resource array for one transmitting antenna and
%   % populate it with PBCH symbols using generated PBCH indices
%   P = 1;
%   txGrid = complex(zeros([240 4 P]),zeros([240 4 P]));
%   txGrid(pbchInd) = pbchTxSym;
%
%   % Perform OFDM modulation
%   txWaveform = ofdmmod(txGrid,256,[22 18 18 18],[1:8 249:256].');
%
%   % Create channel matrix and apply channel to transmitted waveform
%   R = 4;
%   H = fft(eye(max([P R])));
%   H = H(1:P,1:R);
%   H = H / norm(H);
%   rxWaveform = txWaveform * H;
%
%   % Permute the channel matrix to 1-by-1-by-R-by-P and use it to create
%   % the channel estimation grid of size 240-by-4-by-R-by-P
%   hEstGrid = repmat(permute(H.',[3 4 1 2]),[240 4]);
%   nEst = 0.1;
%
%   % Perform OFDM demodulation
%   rxGrid = ofdmdemod(rxWaveform,256,[22 18 18 18],0,[1:8 249:256].');
%
%   % Use nrExtractResources to extract symbols from received and channel
%   % estimate grids in preparation for PBCH decoding
%   [pbchRxSym,pbchHestSym] = nrExtractResources(pbchInd,rxGrid,hEstGrid);
%   figure;
%   plot(pbchRxSym,'o:');
%   title('Received PBCH constellation');
%
%   % Decode PBCH with extracted resource elements
%   [pbchEqSym,csi] = nrEqualizeMMSE(pbchRxSym,pbchHestSym,nEst,1);
%   pbchBits = nrPBCHDecode(pbchEqSym,ncellid,v);
%   figure;
%   plot(pbchEqSym,'o:');
%   title('Equalized PBCH constellation');
%
%   See also nrExtractResources, nrPerfectChannelEstimate,
%   nrPerfectTimingEstimate.

%   Copyright 2018-2019 The MathWorks, Inc.

%#codegen

    narginchk(5,5);

    % Validate received symbols
    fcnName = 'nrEqualizeMMSEIRC';
    validateattributes(rxSym,{'double','single'},{'2d'},fcnName,'RXSYM');

    % Validate estimated channel information
    validateattributes(Hest,{'double','single'},{'3d'},fcnName,'HEST');

    % Validate noise variance
    validateattributes(nVar,{'double','single'},{'scalar','real','nonnegative','finite'},fcnName,'NVAR');

    % Validate the dimensions of rxSym and Hest
    coder.internal.errorIf(size(rxSym,1) ~= size(Hest,1),'nr5g:nrEqualizeMMSE:UnequalNumOfREs',size(rxSym,1),size(Hest,1));
    coder.internal.errorIf(size(rxSym,2) ~= size(Hest,2),'nr5g:nrEqualizeMMSE:UnequalNumOfRxAnts',size(rxSym,2),size(Hest,2));

    % Extract input dimensions
    nRE = size(Hest,1); % Number of resource elements
    R = size(Hest,2);   % Number of receive antennas
    P = size(Hest,3);   % Number of transmission planes

    % Initialize output based on the input dimensions
    out = zeros(nRE,P,'like',rxSym);
    csi = zeros(nRE,P,class(rxSym));

    if (R == 1 && P == 1 && mode == 1)
        % For SISO case
        Hdash = conj(Hest);
        csi = Hdash.*Hest + nVar;
        G = Hdash./csi;
        out = G.*rxSym;
    elseif (R == 1 && P == 1 && mode == 2)
        % For SISO case        
        Hdash = conj(Hest);
        csi = Hdash.*Hest + irc;
        G = Hdash./csi;
        out = G.*rxSym;
    else
        % For non-SISO case

        % Permute to R-by-P-NRE to allow quicker access
        Hest = permute(Hest,[2 3 1]);

        Itx = eye(P);
        n0eye = nVar*Itx;
        for REIdx=1:nRE
            H = Hest(:,:,REIdx);
            temp = (H'*H + n0eye)\Itx;
            csi(REIdx,:) = 1./real(diag(temp));
            G = temp*H';
            out(REIdx,:) = G*(rxSym(REIdx,:).');
        end
    end

end
