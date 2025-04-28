function [sym,ptrs,modulated] = nrPUSCHModified(varargin)
%nrPUSCH Physical uplink shared channel
%   [SYM,PTRSSYM] = nrPUSCH(...) returns a complex matrix SYM containing
%   the physical uplink shared channel (PUSCH) modulation symbols as
%   defined in TS 38.211 Sections 6.3.1.1 - 6.3.1.5. The processing
%   consists of scrambling, symbol modulation, layer mapping, transform
%   precoding, and MIMO precoding. The function also returns the precoded
%   phase tracking reference signal (PT-RS) symbols, PTRSSYM, which are
%   mapped to the resource grid. PT-RS is handled only in the function
%   signatures having configuration objects. For all the other function
%   signatures having no configuration objects, the output PTRSSYM is
%   empty.
%
%   In case of transform precoding enabled, the function maps the data
%   modulated symbols and PT-RS symbols at appropriate locations prior to
%   transform precoding process. The function then performs the transform
%   precoding and MIMO precoding (if applicable), to generate the effective
%   modulated symbols SYM. Also, the effective modulated symbols at the
%   PT-RS locations prior to transform precoding are returned in the output
%   PTRSSYM.
%
%   SYM = nrPUSCH(CW,MODULATION,NLAYERS,NID,RNTI) performs PUSCH modulation
%   given codeword bit vector CW, modulation scheme MODULATION
%   ('pi/2-BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'), number of layers
%   NLAYERS (1...4), scrambling identity NID (0...1023) and radio network
%   temporary identifier RNTI (0...65535). Note that transform precoding
%   and MIMO precoding both are disabled.
%
%   SYM = nrPUSCH(CW,MODULATION,NLAYERS,NID,RNTI,TPRECODE,MRB) allows
%   transform precoding to be enabled or disabled via the parameter
%   TPRECODE (false, true). MRB is the allocated PUSCH bandwidth in
%   resource blocks. Note that MIMO precoding is disabled.
%
%   SYM = nrPUSCH(CW,MODULATION,NLAYERS,NID,RNTI,TPRECODE,MRB,TXSCHEME,NPORTS,TPMI)
%   allows the transmission scheme to be specified via the parameter
%   TXSCHEME ('nonCodebook', 'codebook'). For TXSCHEME = 'codebook', MIMO
%   precoding is performed. NPORTS is the number of antenna ports (1, 2, 4)
%   and TPMI is the transmitted precoding matrix indicator (0...27).
%
%   [SYM,PTRSSYM] = nrPUSCH(CARRIER,PUSCH,CW) returns a complex matrix SYM
%   containing the physical uplink shared channel modulation symbols as
%   defined in TS 38.211 Sections 6.3.1.1 - 6.3.1.5, given the carrier
%   configuration CARRIER, uplink shared channel configuration PUSCH, and
%   codeword bit vector CW. The function also returns the precoded PT-RS
%   symbols, PTRSSYM, which are mapped to the resource grid.
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
%   PUSCH is a physical uplink shared channel configuration object, as
%   described in <a href="matlab:help('nrPUSCHConfig')">nrPUSCHConfig</a> with properties:
%
%   NSizeBWP           - Size of the bandwidth part (BWP) in terms
%                        of number of physical resource blocks (PRBs)
%                        (1...275) (default []). The default value implies
%                        the value is equal to the size of carrier resource
%                        grid
%   NStartBWP          - Starting PRB index of BWP relative to CRB 0
%                        (0...2473) (default []). The default value implies
%                        the value is equal to the start of carrier
%                        resource grid
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
%                        precoding is disabled and the waveform type is
%                        CP-OFDM. 1 indicates that transform precoding is
%                        enabled and waveform type is DFT-s-OFDM
%   TransmissionScheme - Transmission scheme of physical uplink shared
%                        channel ('nonCodebook' (default), 'codebook')
%   NumAntennaPorts    - Number of antenna ports (1 (default), 2, 4)
%   TPMI               - Transmitted precoding matrix indicator (0...27)
%                        (default 0)
%   FrequencyHopping   - Frequency hopping configuration
%                        ('neither' (default), 'intraSlot', 'interSlot')
%   SecondHopStartPRB  - PRB start for second hop relative to the BWP
%                        (0-based) (0...274) (default 1)
%   Interlacing        - Enable interlacing (default false)
%   RBSetIndex         - Resource block set index (default 0)
%   InterlaceIndex     - Interlace indices (0...9) (default 0)
%   NID                - PUSCH scrambling identity (0...1023) (default [])
%   RNTI               - Radio network temporary identifier (0...65535)
%                        (default 1)
%   NRAPID             - Random access preamble index to initialize the
%                        scrambling sequence for msgA on PUSCH (0...63)
%                        (default [])
%   DMRS               - PUSCH-specific DM-RS configuration object, as
%                        described in <a href="matlab:help('nrPUSCHDMRSConfig')">nrPUSCHDMRSConfig</a> with properties:
%       DMRSConfigurationType  - DM-RS configuration type (1 (default), 2).
%                                When transform precoding is enabled, the
%                                value must be 1
%       DMRSTypeAPosition      - Position of first DM-RS OFDM symbol in a
%                                slot (2 (default), 3)
%       DMRSLength             - Number of consecutive DM-RS OFDM symbols
%                                (1 (default), 2). When intra-slot
%                                frequency hopping is enabled, the value
%                                must be 1. Value of 1 indicates
%                                single-symbol DM-RS. Value of 2 indicates
%                                double-symbol DM-RS
%       DMRSAdditionalPosition - Maximum number of DM-RS additional
%                                positions (0...3) (default 0). When
%                                intra-slot frequency hopping is enabled,
%                                the value must be either 0 or 1
%       DMRSPortSet            - DM-RS antenna port set (0...11)
%                                (default []). The default value implies
%                                that the values are in the range from 0 to
%                                NumLayers-1
%       CustomSymbolSet        - Custom DM-RS symbol locations (0-based)
%                                (default []). This property is used to
%                                override the standard defined DM-RS symbol
%                                locations. Each entry corresponds to a
%                                single-symbol DM-RS
%    These properties are applicable, when transform precoding is set to 0:
%       NIDNSCID               - DM-RS scrambling identity (0...65535)
%                                (default []). Use empty ([]) to set the
%                                value to NCellID
%       NSCID                  - DM-RS scrambling initialization
%                                (0 (default), 1)
%    These properties are applicable, when transform precoding is set to 1:
%       NRSID                  - DM-RS scrambling identity (0...1007)
%                                (default []). Use empty ([]) to
%                                set the value to NCellID
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
%       NID                    - PT-RS scrambling identity (0...1007)
%                                (default []). Use empty ([]) to set the
%                                value to DM-RS scrambling identity NRSID
%
%   [SYM,PTRSSYM] = nrPUSCH(...,NAME,VALUE) specifies an additional option
%   as a NAME,VALUE pair to allow control over the datatype of the output
%   symbols SYM and PTRSSYM:
%
%   'OutputDataType' - 'double' for double precision (default)
%                      'single' for single precision
%
%   For operation with shared spectrum channel access for FR1, set
%   Interlacing = true and specify the allocated frequency resources using
%   the RBSetIndex and InterlaceIndex properties of the PUSCH
%   configuration. The PRBSet, FrequencyHopping, and SecondHopStartPRB
%   properties are ignored.
%
%   Example 1:
%   % Generate PUSCH symbols for a codeword of 8064 bits, using 16QAM
%   % modulation and 2 layers, and defaulting to transform precoding
%   % disabled and non-codebook based transmission.
%
%   modulation = '16QAM';
%   nlayers = 2;
%   ncellid = 17;
%   rnti = 111;
%
%   cw = randi([0 1],8064,1);
%   sym = nrPUSCH(cw,modulation,nlayers,ncellid,rnti);
%   size(sym)
%
%   Example 2:
%   % Generate PUSCH symbols for a codeword of 8064 bits, using 256QAM
%   % modulation, 1 layer, transform precoding, 4 antenna ports and
%   % codebook-based transmission.
%
%   modulation = '256QAM';
%   nlayers = 1;
%   ncellid = 17;
%   rnti = 111;
%   transformPrecode = true;
%   MRB = 6;
%   txScheme = 'codebook';
%   nports = 4;
%   TPMI = 1;
%
%   cw = randi([0 1],8064,1);
%   sym = nrPUSCH(cw,modulation,nlayers,ncellid,rnti,transformPrecode,MRB,txScheme,nports,TPMI);
%   size(sym)
%
%   Example 3:
%   % Generate the PUSCH symbols for the configuration specified in Example
%   % 1 with the usage of objects.
%
%   carrier = nrCarrierConfig;
%   carrier.NCellID = 17;
%
%   pusch = nrPUSCHConfig;
%   pusch.Modulation = '16QAM';
%   pusch.NumLayers = 2;
%   pusch.NID = [];
%   pusch.RNTI = 111;
%   pusch.TransformPrecoding = 0;
%   pusch.TransmissionScheme = 'nonCodebook';
%
%   cw = randi([0 1],8064,1);
%   sym = nrPUSCH(carrier,pusch,cw);
%   size(sym)
%
%   See also nrPUSCHDecode, nrPUSCHScramble, nrPUSCHCodebook, nrULSCH,
%   nrPUSCHConfig, nrCarrierConfig, nrIntraCellGuardBandsConfig.
 
%   Copyright 2018-2023 The MathWorks, Inc.

%#codegen
    
    narginchk(3,12);
    
    % Parse and validate inputs
    objSyntax = nargin == 3 || (nargin == 5 && isa(varargin{1},'nrCarrierConfig'));
    if objSyntax
        carrier = varargin{1};          % Carrier configuration object
        pusch = varargin{2};            % PUSCH configuration object
        coder.internal.errorIf(~(isa(carrier,'nrCarrierConfig') && isscalar(carrier)),...
            'nr5g:nrPXSCH:InvalidCarrierInput');
        coder.internal.errorIf(~(isa(pusch,'nrPUSCHConfig') && isscalar(pusch)),...
            'nr5g:nrPUSCH:InvalidPUSCHInput');
        validateConfig(pusch);
        cw = varargin{3};               % Codeword bit vector
        modulation = pusch.Modulation;  % Modulation scheme
        nlayers = pusch.NumLayers;      % Number of layers
        if isempty(pusch.NID)
            % If PUSCH scrambling identity is empty, use physical layer
            % cell identity
            nid = carrier.NCellID;
        else
            nid = pusch.NID(1);
        end
        rnti = pusch.RNTI;              % Radio network temporary identifier
        nrapid = pusch.NRAPID;          % Random access preamble index for msgA on PUSCH
        transformPrecode = pusch.TransformPrecoding;
        if transformPrecode
            %MRB = numel(nr5g.internal.allocatedPhysicalResourceBlocks(carrier,pusch));
            MRB = length(pusch.PRBSet);
        else
            MRB = 1;
        end
        txScheme = pusch.TransmissionScheme;
        nports = pusch.NumAntennaPorts;
        TPMI = pusch.TPMI;
        nvarg = 4;
    else
        narginchk(5,12);
        cw = varargin{1};
        modulation = varargin{2};
        nlayers = varargin{3};
        nid = varargin{4};
        rnti = varargin{5};
        nrapid = []; % msgA on PUSCH support is only for configuration object syntax

        if (nargin>5)
            if (isstring(varargin{6}) || ischar(varargin{6}))
                % ...,NAME,VALUE)
                transformPrecode = false;
                MRB = 1;
                nvarg = 6;
            else
                % ...,TPRECODE,MRB,...
                narginchk(7,12);
                transformPrecode = varargin{6};
                MRB = varargin{7};
                nvarg = 8;
            end
        else
            transformPrecode = false;
            MRB = 1;
            nvarg = 6;
        end
        if (nargin>7)
            if (nargin<10)
                % ...,NAME,VALUE)
                txScheme = 'nonCodebook';
            else
                % ...,TXSCHEME,NPORTS,TPMI,...
                txScheme = varargin{8};
                nports = varargin{9};
                TPMI = varargin{10};
                nvarg = 11;
            end
        else
            txScheme = 'nonCodebook';
        end
    end
    fcnName = 'nrPUSCH';
    modlist = {'pi/2-BPSK','QPSK','16QAM','64QAM','256QAM'};
    validatestring(modulation,modlist,fcnName,'MODULATION');
    validateattributes(nlayers,{'numeric'}, ...
        {'scalar','integer','>=',1,'<=',4},fcnName,'NLAYERS');
    validateattributes(transformPrecode,{'numeric','logical'}, ...
        {'scalar'},fcnName,'TPRECODE');
    schemelist = {'nonCodebook' 'codebook'};
    txScheme = validatestring(txScheme,schemelist,fcnName,'TXSCHEME');
    
    % Scrambling, TS 38.211 Section 6.3.1.1
    scrambled = nrPUSCHScramble(cw,nid,rnti,nrapid);
    
    % Modulation, TS 38.211 Section 6.3.1.2
    modulated = nrSymbolModulate(scrambled,modulation,varargin{nvarg:end});
    
    % Layer mapping, TS 38.211 Section 6.3.1.3
    layered = nrLayerMap(modulated,nlayers);
    
    % Transform precoding, TS 38.211 Section 6.3.1.4
    if objSyntax
        % Generate PT-RS symbols
        ptrsSym = nrPUSCHPTRS(varargin{1},varargin{2},varargin{nvarg:end});
    else
        ptrsSym = zeros(0,1,class(layered));
    end
    if (transformPrecode)
        if ~isempty(ptrsSym) && ~isempty(layered)
            % Generate PT-RS indices
            ptrsInd = nrPUSCHPTRSIndices(carrier,pusch);
            % Get PUSCH resource information
            rmInfo = nr5g.internal.pusch.resourcesInfo(carrier,pusch);
            % Create a temporary variable to map the layered symbols and
            % PT-RS symbols at appropriate locations
            if ~isempty(rmInfo.PRBSet) && ~isempty(rmInfo.PUSCHSymbolSet)
                nRBSC = 12; % Number of subcarriers in a resource block
                tmp = complex(zeros(numel(rmInfo.PRBSet)*nRBSC*numel(rmInfo.PUSCHSymbolSet),...
                    nlayers,class(layered)));
                % Apply scaling factor to PT-RS symbols and map them at the
                % PT-RS locations
                tmp(ptrsInd) = nr5g.internal.pusch.ptrsPowerFactorDFTsOFDM(pusch.Modulation)*ptrsSym;
                dataLogicalInd = (tmp == 0);
                % Validate the input codeword length
                qm = numel(scrambled)/numel(modulated);
                validateattributes(cw,{'double','int8'},...
                    {'numel',nnz(dataLogicalInd)*qm},fcnName,'CW');
                % Map the layered symbols in the locations, other than PT-RS
                tmp(dataLogicalInd) = layered;
            else
                tmp = layered;
            end
        else
            ptrsInd = uint32(zeros(0,1));
            tmp = layered;
        end
        % Perform transform precoding
        transformed = nrTransformPrecode(tmp,MRB);
    else
        ptrsInd = uint32(zeros(0,1));
        transformed = layered;
    end
    
    % MIMO precoding, TS 38.211 Section 6.3.1.5
    if (strcmpi(txScheme,'codebook'))
        W = nrPUSCHCodebook(nlayers,nports,TPMI,transformPrecode);
    else % 'nonCodebook'
        W = eye(nlayers);
    end
    sym = transformed * W;

    % Get the additional PT-RS output
    if transformPrecode
        % In case of transform precoding enabled, directly access PT-RS
        % from the precoded symbols, depending on the ptrsInd variable
        ptrs = sym(ptrsInd(:,1),:);
    else
        % In case of transform precoding disabled, assign ptrsSym
        ptrs = ptrsSym;
    end

end
