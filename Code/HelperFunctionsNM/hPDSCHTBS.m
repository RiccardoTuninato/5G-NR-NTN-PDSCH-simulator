%hPDSCHTBS 5G PDSCH transport block size determination
%   TBS = hPDSCHTBS(CHS,NRED) returns the transport block size associated
%   with a PDSCH, as defined in TS 38.214 5.1.3.2. 
%
%   The PDSCH specific configuration input, CHS, must be a structure
%   including the fields:
%   TargetCoderate      - Code rate used to calculate transport block size
%   PRBSet              - PRBs allocated to the PDSCH (0-based indices)
%   NLayers             - Number of layers
%   Modulation          - Modulation ('QPSK','16QAM','64QAM','256QAM')
%
%   The second input, NRED, is the number of RE allocated per PRB to
%   the PDSCH, accounting for DM-RS, CDM groups and any additional
%   overhead (Xoh-PDSCH).
% 
%   See also hDLSCH, hDLSCHDecode, hPDSCHResources.

%   Copyright 2018 The MathWorks, Inc.

function tbs = hPDSCHTBS(pdsch,nred)

    % Set up the relevant input parameter values
    nprb = length(pdsch.PRBSet);
    rate = pdsch.TargetCodeRate;
    qm = 2*find(strcmpi(pdsch.Modulation,{'QPSK','16QAM','64QAM','256QAM'}));
    nlayers = pdsch.NLayers;

    % Get total number of REs allocated for PDSCH
    nre = min(156,nred)*nprb;
    
    % Obtain intermediate number of information bits 
    ninfo = nre*rate*qm*nlayers;
    
    % TBS determination
    if ninfo <= 3824
        % Get quantized intermediate number of information bits
        n = max(3,floor(log2(ninfo))-6);
        ninfod = max(24,2^n*floor(ninfo/2^n));
        % Search the TBS table
        tbsTab = tbstable();
        for t = 1:length(tbsTab)
            tbs = tbsTab(t);
            if tbs >= ninfod
                return;
            end
        end
    else % ninfo > 3824
        
        % Get quantized intermediate number of information bits
        n = floor(log2(ninfo-24))-5;
        ninfod = max(3840,2^n*round((ninfo-24)/2^n));
        if rate <= 0.25
            C = ceil((ninfod+24)/3816); 
            tbs = 8*C*ceil((ninfod+24)/(8*C))-24;
        else
            if ninfod > 8424
                C = ceil((ninfod+24)/8424); 
                tbs = 8*C*ceil((ninfod+24)/(8*C))-24;
            else
                tbs = 8*ceil((ninfod+24)/8)-24;
            end
        end
    end
    
end

% TS 38.214 table 5.1.3.2-2
function tbs = tbstable()
    % Construct the static table 
    persistent tbstable;
    if isempty(tbstable)
        tbs = [24		336		1288	
              32		352		1320	
              40		368		1352	
              48		384		1416		
              56		408		1480		
              64		432		1544		
              72		456		1608		
              80		480		1672		
              88		504		1736		
              96		528		1800		
             104		552		1864		
             112		576		1928		
             120		608		2024		
             128		640		2088		
             136		672		2152		
             144		704		2216		
             152		736		2280		
             160		768		2408		
             168		808		2472		
             176		848		2536		
             184		888		2600		
             192		928		2664		
             208		984		2728		
             224		1032 	2792		
             240		1064	2856		
             256		1128	2976		
             272		1160	3104		
             288		1192	3240		
             304		1224	3368		
             320		1256	3496];		
        tbstable = tbs(:);
        tbstable = [tbstable; 3624; 3752; 3824];
    end
    tbs = tbstable;
end
