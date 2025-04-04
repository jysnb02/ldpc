classdef (StrictDefaults)TurboDecoder < matlab.System
%TurboDecoder Decode input using a turbo decoder.
%   TURBODEC = comm.TurboDecoder creates a turbo decoder System object,
%   TURBODEC. This object uses the a-posteriori probability (APP)
%   constituent decoder to iteratively decode the parallel-concatenated
%   convolutionally encoded input data.
%
%   TURBODEC = comm.TurboDecoder(Name, Value) creates a turbo decoder
%   object, TURBODEC, with the specified property Name set to the specified
%   Value. You can specify additional name-value pair arguments in any
%   order as (Name1, Value1, ... , NameN, ValueN).
%
%   TURBODEC = comm.TurboDecoder(TRELLIS, INTERLVRINDICES, NUMITER) creates
%   a turbo decoder object, TURBODEC, with the TrellisStructure property
%   set to TRELLIS, the InterleaverIndices property set to INTERLVRINDICES
%   and the NumIterations property set to NUMITER.
%
%   Step method syntax:
%
%   Y = step(TURBODEC, X) decodes the input data, X, using the parallel
%   concatenated convolutional coding scheme that you specify using the
%   TrellisStructure and InterleaverIndices properties. It returns the
%   binary decoded data, Y. Both X and Y are column vectors of double or
%   single precision data type. When the constituent convolutional code
%   represents a rate 1/N code, the step method sets the length of the
%   output vector, Y, to (M-2*numTails)/(2*N-1), where M is the input
%   vector length and numTails is given by
%   log2(TrellisStructure.numStates)*N. The output length, L, is the same
%   as the length of the interleaver indices.
%
%   Y = step(TURBODEC, X, INTERLVRINDICES) uses the INTERLVRINDICES
%   specified as an input. INTERLVRINDICES is a column vector containing
%   integer values from 1 to L with no repeated values. The lengths of the
%   INTERLVRINDICES input and the Y output are the same.
%
%   Y = step(TURBODEC, X, INTERLVRINDICES, ININDICES) uses the ININDICES
%   specified as an input to account for input bit reordering and
%   puncturing. ININDICES is a column vector containing integer values of
%   length equal to the input X length. ININDICES vector values must be
%   relative to the fully encoded data for the coding scheme for all
%   streams, including the tail bits.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(obj, x) and y = obj(x) are
%   equivalent.
%
%   TurboDecoder methods:
%
%   step     - Perform turbo decoding (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create turbo decoder object with same property values
%   isLocked - Locked status (logical)
%
%   TurboDecoder properties:
%
%   TrellisStructure         - Trellis structure of constituent 
%                              convolutional code
%   InterleaverIndicesSource - Source of interleaving indices
%   InterleaverIndices       - Interleaving indices
%   InputIndicesSource       - Source of input indices
%   InputIndices             - Input indices
%   Algorithm                - Decoding algorithm
%   NumScalingBits           - Number of scaling bits
%   NumIterations            - Number of decoding iterations
%
%   % Example:
%   %   Transmit turbo-encoded blocks of data over a BPSK-modulated AWGN
%   %   channel, decode using an iterative turbo decoder and display errors
% 
%   noiseVar = 4; frmLen = 256;
%   s = RandStream('mt19937ar', 'Seed', 11);
%   intrlvrIndices = randperm(s, frmLen);
% 
%   turboEnc  = comm.TurboEncoder('TrellisStructure', poly2trellis(4, ...
%            [13 15 17], 13), 'InterleaverIndices', intrlvrIndices);
%   chan  = comm.AWGNChannel('NoiseMethod', 'Variance', 'Variance', noiseVar);
%   turboDec  = comm.TurboDecoder('TrellisStructure', poly2trellis(4, ...
%            [13 15 17], 13), 'InterleaverIndices', intrlvrIndices, ...
%            'NumIterations', 4);
%   ber = comm.ErrorRate;
% 
%   for frmIdx = 1:8
%       data = randi(s, [0 1], frmLen, 1);
%       encodedData = turboEnc(data);
%       modSignal = pskmod(encodedData,2);
%       receivedSignal = chan(modSignal);
% 
%       % Convert received signal to log-likelihood ratios for decoding
%       receivedBits  = turboDec((-2/(noiseVar/2))*real(receivedSignal));
%     
%       errorStats = ber(data, receivedBits);
%   end
%   fprintf('Error rate = %f\nNumber of errors = %d\nTotal bits = %d\n', ...
%   errorStats(1), errorStats(2), errorStats(3))     
%
%   See also comm.TurboEncoder, comm.APPDecoder.

%   Copyright 2011-2022 The MathWorks, Inc.

%#codegen

properties (Nontunable)
  %TrellisStructure Trellis structure of constituent convolutional code
  %   Specify the trellis as a MATLAB structure that contains the trellis
  %   description of the constituent convolutional code. Use the istrellis
  %   function to check if a structure is a valid trellis structure. The
  %   default is the result of poly2trellis(4, [13 15], 13).
  %
  %   See also istrellis, poly2trellis.
  TrellisStructure = poly2trellis(4, [13 15], 13);
  %InterleaverIndicesSource Source of interleaver indices
  %   Specify the source of the interleaver indices as one of 'Property' |
  %   'Input port'. The default is 'Property'. When you set this property
  %   to 'Input port' the object uses the interleaver indices specified as
  %   an input to the object. When you set this property to
  %   'Property', the object uses the interleaver indices that you specify
  %   in the InterleaverIndices property.
  InterleaverIndicesSource = 'Property';
  %InterleaverIndices Interleaver indices
  %   Specify the mapping used to permute the input bits at the encoder as
  %   a column vector of integers. The default is (64:-1:1).'. This mapping
  %   is a vector with the number of elements equal to length, L,
  %   of the output of the object. Each element must be an integer
  %   between 1 and L, with no repeated values.
  InterleaverIndices = (64:-1:1).';
  %InputIndicesSource Source of input indices
  %   Specify the source of the input indices as one of 'Auto' |
  %   'Property' | 'Input port'. The default is 'Auto'. When you set
  %   this property to 'Input port' the object uses the input indices
  %   specified as an input to the object. When you set this property to
  %   'Property', the object uses the input indices that you specify in
  %   the InputIndices property. The 'Auto' setting evaluates indices that
  %   puncture the second systematic stream and include all tail bits in
  %   the input.
  InputIndicesSource (1,:) char {matlab.system.mustBeMember(InputIndicesSource,{'Auto','Property','Input port'})} = 'Auto';
  %InputIndices Input indices
  %   Specify the bit ordering and puncturing used on the fully encoded
  %   data as a column vector of integers. The length of InputIndices
  %   must be equal to the length of the decoder input. The default is
  %   getTurboIOIndices(64,2,3). Each element specified must be integer
  %   valued.
  InputIndices = getTurboIOIndices(64,2,3);
  %Algorithm Decoding algorithm
  %   Specify the decoding algorithm that the object uses for decoding as
  %   one of 'True APP' | 'Max*' | 'Max'. The default is 'TrueAPP'. When
  %   you set this property to 'True APP', the object implements true a
  %   posteriori probability decoding. When you set this property to any
  %   other value, the object uses approximations to increase the speed of
  %   the computations.
  Algorithm = 'True APP';
  %NumScalingBits Number of scaling bits
  %   Specify the number of bits the constituent decoders use to scale the
  %   input data to avoid losing precision during computations. The
  %   constituent decoders multiply the input by 2^NumScalingBits and
  %   divide the pre-output by the same factor. NumScalingBits must be a
  %   scalar integer between 0 and 8. This property applies when you set
  %   the Algorithm property to 'Max*'. The default is 3.
  NumScalingBits = 3;
  %NumIterations Number of decoding iterations
  %   Specify the number of decoding iterations used for each call to the
  %   object. The default is 6. The object will iterate and provide
  %   updates to the log-likelihood ratios (LLR) of the uncoded output
  %   bits. The output of the object is the hard-decision output of
  %   the final LLR update.
  NumIterations = 6;
end

properties (Constant, Hidden)
    InterleaverIndicesSourceSet = comm.CommonSets.getSet('SpecifyInputs');
    AlgorithmSet = comm.CommonSets.getSet('Algorithm');
end

properties(Access = private, Nontunable)
    % Constituent components
    cAPPDec1;       % decoder1
    cAPPDec2;       % decoder2
    % Commonly used attributes - set for default TrellisStructure
    pK = 1;             % number of uncoded bits
    pN = 2;             % number of coded bits
    pMLen = 3;          % memory length of a constituent encoder
    pNumTails = 6;      % number of tail bits per constituent encoder
end

methods
    % CONSTRUCTOR
    function obj = TurboDecoder(varargin)
        setProperties(obj, nargin, varargin{:}, 'TrellisStructure', ...
                      'InterleaverIndices', 'NumIterations');
    end
    
    function set.InterleaverIndices(obj, value)
        validateattributes(value,...
            {'numeric'}, {'real', 'finite', 'positive', 'integer', ...
            'vector'}, '','InterleaverIndices');     
        obj.InterleaverIndices = value;
    end
    
    function set.InputIndices(obj, value)
            validateattributes(value, {'numeric'}, ...
            {'real', 'finite', 'positive', 'integer', 'vector'}, ...
            '', 'InputIndices'); 

        obj.InputIndices = value;
    end
    
    function set.NumScalingBits(obj, value)
        validateattributes(value, ...
            {'numeric'}, {'real', 'finite', 'nonnegative', 'integer', ...
            'scalar', '>=', 0, '<=', 8}, '', 'NumScalingBits');    
        obj.NumScalingBits = value;
    end
    
    function set.NumIterations(obj, value)
        validateattributes(value,...
            {'numeric'}, {'real', 'finite', 'positive', 'integer', ...
            'scalar'}, '', 'NumIterations'); 
        obj.NumIterations = value;
    end    
end
    
methods(Access = protected)   % System object APIs

    %% Validate inputs
    function validateInputsImpl(~, x, varargin)
        if ~isfloat(x)
            matlab.system.internal.error(...
                'MATLAB:system:invalidInputDataType','X','floating-point');
        else % is a float
            coder.internal.errorIf(isa(x, 'embedded.fi'),...
                                   'comm:TurboDecoder:invalidFiInput');
        end
    end
    
    %% Number of Inputs
    function num = getNumInputsImpl(obj)
        num = 1 + strcmp(obj.InterleaverIndicesSource, 'Input port') ...
            + strcmp(obj.InputIndicesSource, 'Input port');
    end
    
    %% Size propagators
    function flag = isInputSizeMutableImpl(obj,~)
        if strcmp(obj.InterleaverIndicesSource, 'Input port')
            flag = true; % vars only via input port
        else
            flag = false;
        end
    end

    function varargout = isOutputFixedSizeImpl(obj)
        if strcmp(obj.InterleaverIndicesSource, 'Input port')
            varargout = {false};   
        else
            varargout = {true};    
        end
    end
    
    function varargout = getOutputSizeImpl(obj)
        if strcmp(obj.InterleaverIndicesSource, 'Input port')
            inSize = propagatedInputSize(obj, 2);
            varargout = {inSize};  
        else
            varargout = {[length(obj.InterleaverIndices) 1]}; 
        end
    end
    
    %% Type propagators
    function varargout = getOutputDataTypeImpl(obj)
        varargout = {propagatedInputDataType(obj, 1)};   
    end
    
    %% Complexity propagators

    function varargout = isOutputComplexImpl(~)
        varargout = {false}; 
    end

    %% Setup
    function setupImpl(obj, varargin)
        [isOk, status] = istrellis(obj.TrellisStructure);
        coder.internal.errorIf(~isOk,'comm:TurboDecoder:InvalidTrellis', status );

        % Check for mixed modes
        coder.internal.errorIf( ...
           ( strcmp(obj.InterleaverIndicesSource,'Input port') && ...
             strcmp(obj.InputIndicesSource,'Property') ) || ...
           ( strcmp(obj.InterleaverIndicesSource,'Property') && ...
             strcmp(obj.InputIndicesSource,'Input port') ), ...
            'comm:TurboDecoder:InvalidMode');

        if strcmp(obj.Algorithm, 'Max*')
            obj.cAPPDec1 = comm.APPDecoder('TrellisStructure',...
                obj.TrellisStructure,...
                'TerminationMethod', 'Terminated',...
                'Algorithm', obj.Algorithm,...
                'NumScalingBits', obj.NumScalingBits,...
                'CodedBitLLROutputPort', false);
            obj.cAPPDec2 = comm.APPDecoder('TrellisStructure',...
                obj.TrellisStructure,...
                'TerminationMethod', 'Terminated',...
                'Algorithm', obj.Algorithm,...
                'NumScalingBits', obj.NumScalingBits,...
                'CodedBitLLROutputPort', false);
        else
            obj.cAPPDec1 = comm.APPDecoder('TrellisStructure',...
                obj.TrellisStructure,...
                'TerminationMethod', 'Terminated',...
                'Algorithm', obj.Algorithm, ...
                'CodedBitLLROutputPort', false);
            obj.cAPPDec2 = comm.APPDecoder('TrellisStructure',...
                obj.TrellisStructure,...
                'TerminationMethod', 'Terminated',...
                'Algorithm', obj.Algorithm, ...
                'CodedBitLLROutputPort', false);
        end

        obj.pK = log2(obj.TrellisStructure.numInputSymbols);
        coder.internal.errorIf(obj.pK~=1,'comm:TurboEncoder:InvalidTrellis1ByN');
        obj.pN = log2(obj.TrellisStructure.numOutputSymbols);
        obj.pMLen = log2(obj.TrellisStructure.numStates);
        obj.pNumTails = obj.pMLen*(obj.pN);
    end
    
    %% Step
    function y = stepImpl(obj, x, varargin)

        if strcmp(obj.InterleaverIndicesSource, 'Input port')
            interlvrIndices = varargin{1};
            blkLen = length(interlvrIndices); 
        else
            interlvrIndices = get(obj, 'InterleaverIndices');
            blkLen = length(interlvrIndices); 
        end     
        interlvrIndicesCol = interlvrIndices(:);
        coder.internal.assert(isequal((1:blkLen).', sort(interlvrIndicesCol)),...
            'comm:TurboDecoder:InvalidIndices'); 

        if strcmp(obj.InputIndicesSource,'Input port')
            if strcmp(obj.InterleaverIndicesSource,'Input port')
                inIndices = varargin{2};
            else
                inIndices = varargin{1};
            end
        elseif strcmp(obj.InputIndicesSource,'Property')
            inIndices = get(obj,'InputIndices');
        else % Auto
            inIndices = getTurboIOIndices(blkLen,obj.pN,obj.pMLen,'LTE');
        end
        
        if blkLen > 0
            coder.internal.assert(length(x) == length(inIndices),...
                'comm:TurboDecoder:InputLengthMismatch',length(inIndices)); 

            fullInLen = (blkLen+obj.pMLen)*obj.pN*2;
            coder.internal.assert(max(inIndices(:))<=fullInLen,...
                'comm:TurboDecoder:InvalidInIndices',fullInLen);

            typex = class(x);
            % Bit reordering, de-puncturing support
            fullIn = zeros(fullInLen,1,typex);
            fullIn(inIndices,1) = x;
            fullInr = reshape(fullIn,2*obj.pN,blkLen+obj.pMLen);

            Lc1_in = reshape(fullInr(1:obj.pN,:),[],1);
            Lc2_in = reshape(fullInr(obj.pN+1:2*obj.pN,:),[],1);

            Lu1_in = zeros(blkLen+obj.pMLen,1,typex);
                        
            % Turbo Decode
            out1 = zeros(blkLen,1,typex);
            for iterIdx = 1:obj.NumIterations
                Lu1_out = step(obj.cAPPDec1,Lu1_in,Lc1_in);
                tmp = Lu1_out((1:blkLen).', 1);
                tmp2 = tmp(:);
                Lu2_out = step(obj.cAPPDec2, ...
                    [tmp2(interlvrIndicesCol); zeros(obj.pMLen,1,typex)], Lc2_in);
                
                out1(interlvrIndicesCol, 1) = Lu2_out((1:blkLen).',1);
                Lu1_in = [out1; zeros(obj.pMLen,1,typex)];
            end
            
            % Calculate llr and decoded bits - for the final iteration
            llr = out1 + tmp2;
            y = cast((llr>=0),typex);
        else
            % allow empty input
            y = [];
        end
    end

    %%
    function releaseImpl(obj)
        release(obj.cAPPDec1);
        release(obj.cAPPDec2);
    end
    
    function flag = isInactivePropertyImpl(obj, prop)
        flag = false;
        switch prop
            case 'InterleaverIndices'
                flag = strcmp(obj.InterleaverIndicesSource,'Input port');
            case 'InputIndices'
                if ( strcmp(obj.InputIndicesSource,'Input port') || ...
                     strcmp(obj.InputIndicesSource,'Auto') )    
                    flag = true;
                end
            case 'NumScalingBits'
                if ~strcmp(obj.Algorithm,'Max*')
                    flag = true;
                end
        end
    end
  
    function s = saveObjectImpl(obj)
        s = saveObjectImpl@matlab.System(obj);
        if isLocked(obj)
            s.pK        = obj.pK;
            s.pN        = obj.pN;
            s.pMLen     = obj.pMLen;
            s.pNumTails = obj.pNumTails;
            s.cAPPDec1  = matlab.System.saveObject(obj.cAPPDec1);
            s.cAPPDec2  = matlab.System.saveObject(obj.cAPPDec2);
        end
    end
  
    function loadObjectImpl(obj, s, wasLocked)
        if wasLocked
            obj.pK        = s.pK;
            obj.pN        = s.pN;
            obj.pMLen     = s.pMLen;
            obj.pNumTails = s.pNumTails;
            obj.cAPPDec1  = matlab.System.loadObject(s.cAPPDec1);
            obj.cAPPDec2  = matlab.System.loadObject(s.cAPPDec2);
        end
        % Call the base class method
        loadObjectImpl@matlab.System(obj, s);
    end
    
    function icon = getIconImpl(~)        
        icon = sprintf('Turbo\nDecoder');
    end
  
    function varargout = getInputNamesImpl(obj)
        varargout = cell(1, getNumInputs(obj));        
        if strcmp(obj.InterleaverIndicesSource, 'Input port')
            varargout{1} = 'In';
            varargout{2} = 'IntrInd';
            if strcmp(obj.InputIndicesSource,'Input port')
                varargout{3} = 'InInd';
            end            
        else
            if strcmp(obj.InputIndicesSource,'Input port')
                varargout{1} = 'In';
                varargout{2} = 'InInd';
            else            
                varargout{1} = '';
            end
        end
    end
    
    function varargout = getOutputNamesImpl(~)
        % Always have no label for the output port
        varargout = {''};
    end
    
end

methods(Static, Access = protected)
% System Block customization functions
  function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'comm.TurboDecoder', ...
            'Title','comm:TurboDecoder:TurboDecoderTitle', ...
            'Text','comm:TurboDecoder:TurboDecoderDescription' );
  end    
    
    % Specify aliases to existing block parameters
    function groups = getPropertyGroupsImpl
       

        propsList = {'TrellisStructure', 'InterleaverIndicesSource', 'InterleaverIndices', ...
                              'InputIndicesSource','InputIndices','Algorithm','NumScalingBits','NumIterations'};
        paramsList = cell(1,numel(propsList));
        paramsList{1} = matlab.system.display.internal.Property(propsList{1},'Description',['comm:TurboDecoder:',propsList{1}],'Alias', 'trellis');
        paramsList{2} = matlab.system.display.internal.Property(propsList{2},'Description',['comm:TurboDecoder:',propsList{2}]);
        paramsList{3} = matlab.system.display.internal.Property(propsList{3},'Description',['comm:TurboDecoder:',propsList{3}],'Alias', 'intrIndices');
        paramsList{4} = matlab.system.display.internal.Property(propsList{4},'Description',['comm:TurboDecoder:',propsList{4}]);
        paramsList{5} = matlab.system.display.internal.Property(propsList{5},'Description',['comm:TurboDecoder:',propsList{5}]);
        paramsList{6} = matlab.system.display.internal.Property(propsList{6},'Description',['comm:TurboDecoder:',propsList{6}],'Alias', 'algorithm'); 
        paramsList{7} = matlab.system.display.internal.Property(propsList{7},'Description',['comm:TurboDecoder:',propsList{7}],'Alias', 'numScaleBits'); 
        paramsList{8} = matlab.system.display.internal.Property(propsList{8},'Description',['comm:TurboDecoder:',propsList{8}],'Alias', 'numIter'); 

        params = matlab.system.display.Section(...
                'Title', 'comm:TurboDecoder:Parameters', ...
                'PropertyList',paramsList);


        groups = params;
    end
end

end
