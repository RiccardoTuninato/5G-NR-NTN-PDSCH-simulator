function [cw,symbols, detransformedPTRS, demod, xInd,yInd] = nrPUSCHDecodeModified(varargin)
%nrPUSCHDecode Physical uplink shared channel decoding
%   [CW,SYMBOLS] = nrPUSCHDecode(...) returns a vector of soft bits CW and
%   vector of received constellation symbols SYMBOLS, resulting from the
%   inverse operation of physical uplink shared channel processing as
%   defined in TS 38.211 Section 6.3.1. The processing consists of MIMO
%   deprecoding, transform deprecoding, layer demapping, symbol
%   demodulation, and descrambling. Note that this function performs
%   descrambling with uplink control information (UCI) placeholder only
%   when the function signature has input arguments related to control
%   information. For all the other function signatures, the function
%   performs only data descrambling. Also, note that, for the function
%   signatures with configuration objects having transform precoding
%   enabled, the function assumes input symbols contain both the data and
%   PT-RS symbols (if applicable) and uses only data symbols for further
%   processing. For all the other function signatures, the function assumes
%   input symbols are only data symbols.
%
%   [CW,SYMBOLS] = nrPUSCHDecode(SYM,MODULATION,NID,RNTI) performs physical
%   uplink shared channel (PUSCH) demodulation given received PUSCH symbols
%   SYM, modulation scheme MODULATION ('pi/2-BPSK', 'QPSK', '16QAM',
%   '64QAM', '256QAM'), scrambling identity NID (0...1023) and radio
%   network temporary identifier RNTI (0...65535). Note that MIMO
%   deprecoding and transform deprecoding both are disabled.
%
%   [CW,SYMBOLS] = nrPUSCHDecode(SYM,MODULATION,NID,RNTI,NVAR) allows the
%   variance of additive white Gaussian noise on the received PUSCH symbols
%   to be specified via the parameter NVAR, a nonnegative real scalar. The
%   default value is 1e-10.
%
%   [CW,SYMBOLS] = nrPUSCHDecode(SYM,MODULATION,NID,RNTI,NVAR,TPRECODE,MRB)
%   allows transform deprecoding to be enabled or disabled via the
%   parameter TPRECODE (false, true). MRB is the allocated PUSCH bandwidth
%   in resource blocks. Note that MIMO deprecoding is disabled.
%
%   [CW,SYMBOLS] = nrPUSCHDecode(SYM,MODULATION,NID,RNTI,NVAR,TPRECODE,MRB,TXSCHEME,NLAYERS,TPMI)
%   allows the transmission scheme to be specified via the parameter
%   TXSCHEME ('nonCodebook', 'codebook'). For TXSCHEME = 'codebook', MIMO
%   deprecoding is performed. NLAYERS is the number of transmission layers
%   (1...4) and TPMI is the transmitted precoding matrix indicator
%   (0...27).
%
%   [CW,SYMBOLS] = nrPUSCHDecode(CARRIER,PUSCH,SYM,NVAR) returns a vector
%   of soft bits CW and vector of received constellation symbols SYMBOLS
%   resulting from the inverse operation of the physical uplink shared
%   channel processing as defined in TS 38.211 Section 6.3.1, for the
%   specified carrier configuration CARRIER and the uplink shared channel
%   configuration PUSCH. The SYM is the received symbols for each layer and
%   NVAR is the optional noise variance scaling factor of the soft bits.
%   The default value of NVAR is 1e-10.
%
%   [CW,SYMBOLS] = nrPUSCHDecode(CARRIER,PUSCH,TCR,TBS,OACK,OCSI1,OCSI2,SYM,NVAR)
%   returns a vector of soft bits CW and vector of received constellation
%   symbols SYMBOLS resulting from the inverse operation of the physical
%   uplink shared channel processing as defined in TS 38.211 Section 6.3.1.
%   TCR is a scalar with value in between 0 and 1, specifying target code
%   rate. TBS is a scalar nonnegative value, specifying the transport block
%   size. OACK is a scalar nonnegative value, specifying the HARQ-ACK
%   payload length. OCSI1 is a scalar nonnegative value, specifying the CSI
%   part 1 payload length. OCSI2 is a scalar nonnegative value, specifying
%   the CSI part 2 payload length. SYM is the received PUSCH symbols for
%   each layer and NVAR is optional noise variance. The default value of
%   NVAR is 1e-10. Note that this syntax handles the UCI placeholders in
%   the descrambling.
%
%   CARRIER is a carrier configuration object, as described in
%   <a href="matlab:help('nrCarrierConfig')">nrCarrierConfig</a> with properties:
%
%   NCellID           - Physical layer cell identity (0...1007) (default 1)
%   SubcarrierSpacing - Subcarrier spacing in kHz
%                       (15 (default), 30, 60, 120, 240, 480, 960)
%   CyclicPrefix      - Cyclic prefix ('normal' (default), 'extended')
%   NSizeGrid         - Number of resource blocks in carrier resource
%                       grid (1...275) (default 52)
%   NStartGrid        - Start of carrier resource grid relative to common
%                       resource block 0 (CRB 0) (0...2199) (default 0)
%   NSlot             - Slot number (default 0)
%
%   PUSCH is the physical uplink shared channel configuration object, as
%   described in <a href="matlab:help('nrPUSCHConfig')">nrPUSCHConfig</a> with properties:
%
%   Modulation         - Modulation scheme
%                        ('QPSK' (default), 'pi/2-BPSK', '16QAM', '64QAM', '256QAM')
%   NumLayers          - Number of transmission layers (1...4) (default 1)
%   MappingType        - Mapping type of physical uplink shared channel
%                        ('A' (default), 'B')
%   SymbolAllocation   - Symbol allocation of physical uplink shared
%                        channel (default [0 14]). This property is a
%                        two-element vector. First element represents the
%                        start of OFDM symbol in a slot. Second element
%                        represents the number of contiguous OFDM symbols
%   PRBSet             - PRBs allocated for physical uplink shared channel
%                        within a BWP (0-based) (default 0:51)
%   TransformPrecoding - Flag to enable transform precoding
%                        (0 (default), 1). 0 indicates that transform
%                        precoding is disabled, and the waveform type is
%                        CP-OFDM. 1 indicates that transform precoding is
%                        enabled, and waveform type is DFT-s-OFDM
%   TransmissionScheme - Transmission scheme of physical uplink shared
%                        channel ('nonCodebook' (default), 'codebook')
%   TPMI               - Transmitted precoding matrix indicator (0...27)
%                        (default 0)
%   FrequencyHopping   - Frequency hopping configuration
%                        ('neither' (default), 'intraSlot', 'interSlot')
%   Interlacing        - Enable interlacing (default false)
%   RBSetIndex         - Resource block set index (default 0)
%   InterlaceIndex     - Interlace indices (0...9) (default 0)
%   BetaOffsetACK      - Beta offset of HARQ-ACK (default 20)
%   BetaOffsetCSI1     - Beta offset of CSI part 1 (default 6.25)
%   BetaOffsetCSI2     - Beta offset of CSI part 2 (default 6.25)
%   UCIScaling         - UCI scaling factor. The nominal value is one of
%                        {0.5, 0.65, 0.8, 1 (default)}
%   NID                - Scrambling identity (0...1023) (default []). Use
%                        empty ([]) to set the value to NCellID
%   RNTI               - Radio network temporary identifier (0...65535)
%                        (default 1)
%   NRAPID             - Random access preamble index to initialize the
%                        scrambling sequence for msgA on PUSCH (0...63)
%                        (default [])
%   DMRS               - PUSCH-specific DM-RS configuration object, as
%                        described in <a href="matlab:help('nrPUSCHDMRSConfig')">nrPUSCHDMRSConfig</a> with properties:
%       DMRSConfigurationType   - DM-RS configuration type (1 (default), 2).
%                                 When transform precoding is enabled, the
%                                 value must be 1
%       DMRSTypeAPosition       - Position of first DM-RS OFDM symbol in a
%                                 slot (2 (default), 3)
%       DMRSLength              - Number of consecutive DM-RS OFDM symbols
%                                 (1 (default), 2). When intra-slot
%                                 frequency hopping is enabled, the value
%                                 must be 1. Value of 1 indicates
%                                 single-symbol DM-RS. Value of 2 indicates
%                                 double-symbol DM-RS
%       DMRSAdditionalPosition  - Maximum number of DM-RS additional
%                                 positions (0...3) (default 0). When
%                                 intra-slot frequency hopping is enabled,
%                                 the value must be either 0 or 1
%       DMRSPortSet             - DM-RS antenna port set (0...11)
%                                 (default []). The default value implies
%                                 that the values are in the range from 0
%                                 to NumLayers-1
%       CustomSymbolSet         - Custom DM-RS symbol locations (0-based)
%                                 (default []). This property is used to
%                                 override the standard defined DM-RS
%                                 symbol locations. Each entry corresponds
%                                 to a single-symbol DM-RS
%       NumCDMGroupsWithoutData - Number of CDM groups without data
%                                 (1...3) (default 2). When transform
%                                 precoding is enabled, the value must be 2
%   EnablePTRS         - Enable or disable the PT-RS configuration
%                        (0 (default), 1)
%   PTRS               - PUSCH-specific PT-RS configuration object, as
%                        described in <a href="matlab:help('nrPUSCHPTRSConfig')">nrPUSCHPTRSConfig</a> with properties:
%       TimeDensity            - PT-RS time density (1 (default), 2, 4)
%    These properties are applicable, when transform precoding is set to 0:
%       FrequencyDensity       - PT-RS frequency density (2 (default), 4)
%       REOffset               - Resource element offset
%                                ('00' (default), '01', '10', '11')
%       PTRSPortSet            - PT-RS antenna port set (default []). The
%                                default value of empty ([]) implies the
%                                value is equal to the lowest DM-RS antenna
%                                port configured
%    These properties are applicable, when transform precoding is set to 1:
%       NumPTRSSamples         - Number of PT-RS samples (2 (default), 4)
%       NumPTRSGroups          - Number of PT-RS groups (2 (default), 4, 8)
%
%   For operation with shared spectrum channel access for FR1, set
%   Interlacing = true and specify the allocated frequency resources using
%   the RBSetIndex and InterlaceIndex properties of the PUSCH
%   configuration. The PRBSet, FrequencyHopping, and SecondHopStartPRB
%   properties are ignored.
%
%   Example 1:
%   % Generate PUSCH symbols for a codeword of 8064 bits, using 256QAM
%   % modulation and 2 layers (defaulting to transform precoding
%   % disabled and non-codebook based transmission), and perform PUSCH 
%   % demodulation.
%
%   modulation = '256QAM';
%   nlayers = 2;
%   ncellid = 17;
%   rnti = 111;
%
%   cw = randi([0 1],8064,1);
%   sym = nrPUSCH(cw,modulation,nlayers,ncellid,rnti);
%   rxcw = nrPUSCHDecode(sym,modulation,ncellid,rnti);
%   isequal(cw,double(rxcw<0))
%
%   Example 2:
%   % Generate PUSCH symbols for a codeword of 8064 bits, using QPSK
%   % modulation, 1 layer, transform precoding, 4 antenna ports and
%   % codebook-based transmission, and perform PUSCH demodulation.
%
%   modulation = 'QPSK';
%   nlayers = 1;
%   ncellid = 17;
%   rnti = 111;
%   transformPrecode = true;
%   MRB = 6;
%   txScheme = 'codebook';
%   nports = 4;
%   TPMI = 1;
%   nVar = 0; % noise variance assumed to be zero
%
%   cw = randi([0 1],8064,1);
%   sym = nrPUSCH(cw,modulation,nlayers,ncellid,rnti,transformPrecode,MRB,txScheme,nports,TPMI);
%   rxcw = nrPUSCHDecode(sym,modulation,ncellid,rnti,nVar,transformPrecode,MRB,txScheme,nlayers,TPMI);
%   isequal(cw,double(rxcw<0))
%
%   Example 3:
%   % Generate the PUSCH symbols and decode the data bits for the same
%   % configuration specified in Example 1 with the usage of objects.
%
%   carrier = nrCarrierConfig;
%   carrier.NCellID = 17;
%
%   pusch = nrPUSCHConfig;
%   pusch.Modulation = '256QAM';
%   pusch.NumLayers = 2;
%   pusch.NID = [];
%   pusch.RNTI = 111;
%   pusch.TransformPrecoding = 0;
%   pusch.TransmissionScheme = 'nonCodebook';
%
%   cw = randi([0 1],8064,1);
%   sym = nrPUSCH(carrier,pusch,cw);
%   rxcw = nrPUSCHDecode(carrier,pusch,sym);
%   isequal(cw,double(rxcw<0))
%
%   See also nrPUSCH, nrPUSCHCodebook, nrPUSCHDescramble, nrULSCHDecoder,
%   nrPUSCHConfig, nrCarrierConfig.

%   Copyright 2018-2023 The MathWorks, Inc.

%#codegen
    
    narginchk(3,10);
    
    % Parse and validate inputs
    objSyntax = any(nargin == [3 8 9]) || ...
            (nargin == 4 && isa(varargin{1},'nrCarrierConfig'));
    if objSyntax
        carrier = varargin{1};           % Carrier configuration object
        pusch = varargin{2};             % PUSCH configuration object
        coder.internal.errorIf(~(isa(carrier,'nrCarrierConfig') && isscalar(carrier)),...
            'nr5g:nrPXSCH:InvalidCarrierInput');
        coder.internal.errorIf(~(isa(pusch,'nrPUSCHConfig') && isscalar(pusch)),...
            'nr5g:nrPUSCH:InvalidPUSCHInput');
        validateConfig(pusch);
        if nargin >= 8
            sym = varargin{8};           % Received PUSCH symbols
            if nargin == 9
                % Noise variance
                nVar = varargin{9};
            else
                nVar = 1e-10;
            end
        else
            sym = varargin{3};           % Received PUSCH symbols 
            if nargin == 4
                % Noise variance
                nVar = varargin{4};
            else
                nVar = 1e-10;
            end
        end
        modulation = pusch.Modulation;   % Modulation scheme
        rnti = pusch.RNTI;               % Radio network temporary identifier
        nrapid = pusch.NRAPID;           % Random access preamble index for msgA on PUSCH
        if isempty(pusch.NID)
            % If PUSCH scrambling identity is empty, use physical layer
            % cell identity
            nid = carrier.NCellID;
        else
            nid = pusch.NID(1);
        end
        transformPrecode = pusch.TransformPrecoding;
        if transformPrecode
            %MRB = numel(nr5g.internal.allocatedPhysicalResourceBlocks(carrier,pusch));
            MRB = length(pusch.PRBSet);
        else
            MRB = 1;
        end
        txScheme = pusch.TransmissionScheme;
        nlayers = pusch.NumLayers;
        TPMI = pusch.TPMI;
    else
        sym = varargin{1};
        modulation = varargin{2};
        nid = varargin{3};
        rnti = varargin{4};
        nrapid = []; % msgA on PUSCH support is only for configuration object syntax

        if nargin>4
            nVar = varargin{5};
        else
            nVar = 1e-10;
        end
        if nargin>5
            narginchk(7,10);
            transformPrecode = varargin{6};
            MRB = varargin{7};
        else
            transformPrecode = false;
            MRB = 1;
        end
        if (nargin>7)
            narginchk(10,10);
            txScheme = varargin{8};
            nlayers = varargin{9};
            TPMI = varargin{10};
        else
            txScheme = 'nonCodebook';
            nlayers = max([size(sym,2) 1]);
            TPMI = 0;
        end
    end
    fcnName = 'nrPUSCHDecode';
    validateattributes(sym,{'double','single'},{'finite'},fcnName,'SYM');
    modlist = {'pi/2-BPSK','QPSK','16QAM','64QAM','256QAM'};
    validatestring(modulation,modlist,fcnName,'MODULATION');
    validateattributes(nlayers,{'numeric'}, ...
        {'scalar','integer','>=',1,'<=',4},fcnName,'NLAYERS');
    validateattributes(transformPrecode,{'numeric','logical'}, ...
        {'scalar'},fcnName,'TPRECODE');
    schemelist = {'nonCodebook' 'codebook'};
    txScheme = validatestring(txScheme,schemelist,fcnName,'TXSCHEME');
    
    % MIMO deprecoding, TS 38.211 Section 6.3.1.5
    if (isempty(sym))
        deprecoded = zeros(0,1);
    else
        if (strcmpi(txScheme,'codebook'))
            nports = size(sym,2);
            W = nrPUSCHCodebook(nlayers,nports,TPMI,transformPrecode);
        else % 'nonCodebook'
            W = eye(nlayers);
        end
        deprecoded = sym * pinv(W);
    end
    
    detransformedPTRS = [];
    % Transform deprecoding, TS 38.211 Section 6.3.1.4
    if (transformPrecode) && ~isempty(deprecoded)
        detransformed = nrTransformDeprecode(deprecoded,MRB);
        if objSyntax && varargin{2}.EnablePTRS
            % Get the data symbols ignoring PT-RS symbols, if PT-RS is
            % enabled
            pusch = varargin{2};
            % Generate PT-RS indices
            ptrsInd = nrPUSCHPTRSIndices(carrier,pusch);
            % Get PUSCH resource information
            rmInfo = nr5g.internal.pusch.resourcesInfo(carrier,pusch);
            % Initialize a temporary variable with PUSCH allocation,
            % ignoring DM-RS OFDM symbols
            numDataSym = numel(rmInfo.PRBSet)*12*numel(rmInfo.PUSCHSymbolSet);
            temp = zeros(numDataSym,1);
            temp(ptrsInd(:,1)) = 1;
            % Ensure the number of rows in the input SYM is equal to the
            % number of data symbols containing PT-RS
            validateattributes(sym,{'double','single'},{'nrows',numDataSym},fcnName,'SYM');
            
            % Extract data symbols, excluding PT-RS
            detransformedData = detransformed(temp == 0,:);

            % Extract data symbols, excluding PT-RS
            detransformedPTRS = detransformed(temp == 1,:);

        else
            % When function signature is without carrier input or when
            % PT-RS is disabled with carrier input, all the symbols are
            % assumed to be data symbols
            detransformedData = detransformed;
        end
    else
        % When transform precoding is disabled, all the symbols are data
        % symbols
        detransformedData = deprecoded;
    end

    % Layer demapping, TS 38.211 Section 6.3.1.3
    if (isempty(detransformedData))
        symbols = detransformedData;
    else
        symbols = nrLayerDemap(detransformedData);
        symbols = symbols{1};
    end

    % Demodulation, TS 38.211 Section 6.3.1.2
    demod = nrSymbolDemodulate(symbols(:),modulation,nVar);
    
    % Descrambling, TS 38.211 Section 6.3.1.1
    if (nargin == 8 || nargin == 9) && ~isempty(demod)
        % Parse and validate inputs
        tcr   = varargin{3};
        tbs   = varargin{4};
        oack  = varargin{5};
        ocsi1 = varargin{6};
        ocsi2 = varargin{7};
        validateattributes(tcr,{'numeric'},{'real','scalar','>',0,'<',1},fcnName,'TCR');
        validateattributes(tbs,{'numeric'},{'scalar','integer','nonnegative'},fcnName,'TBS');
        validateattributes(oack,{'numeric'},{'scalar','integer','nonnegative'},fcnName,'OACK');
        validateattributes(ocsi1,{'numeric'},{'scalar','integer','nonnegative'},fcnName,'OCSI1');
        validateattributes(ocsi2,{'numeric'},{'scalar','integer','nonnegative'},fcnName,'OCSI2');
        % Get UCI placeholder locations
        [xInd,yInd] = getPlaceHolderLocations(pusch,tcr,tbs,oack,ocsi1,ocsi2);
        % Pass only the UCI placeholder locations that are within the
        % length of demodulated bits, to the nrPUSCHDescramble function
        demodLen = length(demod);
        xInd = xInd(xInd <= demodLen);
        yInd = yInd(yInd <= demodLen);
    else
        xInd = zeros(0,1);
        yInd = zeros(0,1);
    end
    cw = nrPUSCHDescramble(demod,nid,rnti,nrapid,xInd,yInd);

end

function [xInd,yInd] = getPlaceHolderLocations(pusch,tcr,tbs,oack,ocsi1,ocsi2)
%getPlaceHolderLocations UCI placeholder locations
%
%   [XIND,YIND] = getPlaceHolderLocations(PUSCH,TCR,TBS,OACK,OCSI1,OCSI2)
%   returns the UCI placeholder 'x' locations, XIND, and UCI placeholder
%   'y' locations, YIND, by encoding the HARQ-ACK, CSI part 1, and CSI part 2,
%   and then multiplexing the codeword.

    % Get the rate-match and resource information
    [info,resInfo] = nr5g.internal.pusch.getULSCHInfo(pusch,tcr,tbs,oack,ocsi1,ocsi2);
    resInfo.Modulation = pusch.Modulation;
    resInfo.NumLayers = pusch.NumLayers;

    % Check if any of the UCI payload length is 1 or 2, and get the UCI
    % placeholder indices
    ackFlag = (oack <= 2) && (oack > 0) && (info.GACK > 0);
    csi1Flag = (ocsi1 <= 2) && (ocsi1 > 0) && (info.GCSI1 > 0);
    csi2Flag = (ocsi2 <= 2) && (ocsi2 > 0) && (info.GCSI2 > 0);
    if any([ackFlag csi1Flag csi2Flag])
        % Perform data and control mapping
        [cwTemp,indInfo] = nr5g.internal.pusch.dataAndControlMapping(...
            resInfo,info.GULSCH,info.GACK,info.GCSI1,info.GCSI2,info.GACKRvd);
        % Encode the UCI (HARQ-ACK, CSI part 1, CSI part 2) with data
        % having all ones, if payload length is 1 or 2. And, map the coded
        % UCI bits into the codeword
        cw = zeros(size(cwTemp));
        if ackFlag
            cack = uciEncode(ones(oack,1),info.GACK,pusch.Modulation);
            cw(indInfo.ACKIndices) = cack(1:length(indInfo.ACKIndices));
        end
        if csi1Flag
            ccsi1 = uciEncode(ones(ocsi1,1),info.GCSI1,pusch.Modulation);
            cw(indInfo.CSI1Indices) = ccsi1(1:length(indInfo.CSI1Indices));
        end
        if csi2Flag
            if ~isempty(indInfo.CSI2Indices)
                if ~isempty(indInfo.CSI2ACKIndices)
                    % Find the logical indices to access coded CSI part 2
                    % that are part of actual codeword
                    sortCSI2Ind = sort([indInfo.CSI2Indices;indInfo.CSI2ACKIndices]);
                    csi2LogicalInd = ismember(sortCSI2Ind,indInfo.CSI2Indices);
                else
                    csi2LogicalInd = true(length(indInfo.CSI2Indices),1);
                end
                ccsi2 = uciEncode(ones(ocsi2,1),info.GCSI2,pusch.Modulation);
                cw(indInfo.CSI2Indices) = ccsi2(csi2LogicalInd);
            end
        end
        % Get the UCI placeholder indices
        xInd = find(cw == -1);
        yInd = find(cw == -2);
    else
        xInd = zeros(0,1);
        yInd = zeros(0,1);
    end

end

function out = uciEncode(in,E,modulation)
%uciEncode Performs UCI encoding for small block lengths of 1 or 2 bits
%
%   OUT = uciEncode(IN,E,MODULATION) returns the output OUT by encoding the
%   input bits IN of length 1 or 2 and repeats the encoded bits up to
%   length E. The encoding is according to TS 38.212 Sections 5.3.3.1 and
%   5.3.3.2.

    coded = nr5g.internal.smallEncode12(in,modulation);
    out = coded(mod(0:E-1,length(coded))+1,1);

end
