classdef (StrictDefaults)TurboEncoder < matlab.System
%TurboEncoder Encode binary data using a turbo encoder.
%   TURBOENC = comm.TurboEncoder creates a System object, TURBOENC, that
%   encodes binary data using a turbo encoder.
%
%   TURBOENC = comm.TurboEncoder(Name, Value) creates a turbo encoder
%   object, TURBOENC, with the specified property Name set to the specified
%   Value. You can specify additional name-value pair arguments in any
%   order as (Name1, Value1, ... , NameN, ValueN).
%
%   TURBOENC = comm.TurboEncoder(TRELLIS, INTERLVRINDICES) creates a turbo
%   encoder object, TURBOENC, with the TrellisStructure property set to
%   TRELLIS, and the InterleaverIndices property set to INTERLVRINDICES.
%
%   Step method syntax:
%
%   Y = step(TURBOENC, X) encodes the binary data, X, using the parallel
%   concatenated convolutional encoding scheme that you specify by the
%   TrellisStructure and InterleaverIndices properties. It returns the
%   encoded data, Y. Both X and Y are column vectors of data type numeric,
%   logical, or unsigned fixed point with word length 1 (fi object).
%   When the constituent convolutional encoder represents a rate 1/N code,
%   the step method sets the length of the output vector, Y, to L*(2*N-1)+
%   2*numTails, where L is the input vector length and numTails is given by
%   log2(TrellisStructure.numStates)*N. The tail bits, due to the
%   termination, are appended at the end after the input bits are encoded.
%
%   Y = step(TURBOENC, X, INTERLVRINDICES) uses the INTERLVRINDICES
%   specified as an input. INTERLVRINDICES is a column vector containing
%   integer values from 1 to L with no repeated values. The length of the
%   data input X and INTERLVRINDICES input must be the same and equal to L.
%
%   Y = step(TURBOENC, X, INTERLVRINDICES, OUTINDICES) uses the OUTINDICES
%   specified as an input to account for bit reordering and puncturing.
%   OUTINDICES is a column vector containing integer values of length equal
%   to the desired output length. OUTINDICES vector values must be relative
%   to the fully encoded data for the coding scheme for all streams,
%   including the tail bits.
%
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(obj, x) and y = obj(x) are
%   equivalent.
%
%   TurboEncoder methods:
%
%   step     - Perform turbo encoding (see above)
%   release  - Allow property value and input characteristics changes
%   clone    - Create turbo encoder object with same property values
%   isLocked - Locked status (logical)
%
%   TurboEncoder properties:
%
%   TrellisStructure         - Trellis structure of constituent 
%                              convolutional code
%   InterleaverIndicesSource - Source of interleaving indices
%   InterleaverIndices       - Interleaving indices
%   OutputIndicesSource      - Source of output indices
%   OutputIndices            - Output indices
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
%   See also comm.TurboDecoder, comm.ConvolutionalEncoder.

%   Copyright 2011-2022 The MathWorks, Inc.

%#codegen

properties (Nontunable)
  %TrellisStructure Trellis structure of constituent convolutional code
  %   Specify the trellis as a MATLAB structure that contains the trellis
  %   description of the constituent convolutional code. Use the istrellis
  %   function to check if a structure is a valid trellis structure. The
  %   default is the result of poly2trellis(4,[13 15],13).
  %
  %   See also istrellis, poly2trellis.
  TrellisStructure = poly2trellis(4,[13 15],13);
  
  %InterleaverIndicesSource Source of interleaver indices
  %   Specify the source of the interleaver indices as one of 'Property' |
  %   'Input port'. The default is 'Property'. When you set this property
  %   to 'Input port' the object uses the interleaver indices specified as
  %   an input to the object. When you set this property to
  %   'Property', the object uses the interleaver indices that you specify
  %   in the InterleaverIndices property.
  InterleaverIndicesSource = 'Property';

  %InterleaverIndices Interleaver indices
  %   Specify the mapping used to permute the input bits as a column vector
  %   of integers. The default is (64:-1:1).'. This mapping is a vector
  %   with the number of elements equal to length, L, of the input to the
  %   object. Each element must be an integer between 1 and L, with no
  %   repeated values.
  InterleaverIndices = (64:-1:1).';

  %OutputIndicesSource Source of output indices
  %   Specify the source of the output indices as one of 'Auto' |
  %   'Property' | 'Input port'. The default is 'Auto'. When you set
  %   this property to 'Input port' the object uses the output indices
  %   specified as an input to the object. When you set this property to
  %   'Property', the object uses the output indices that you specify in
  %   the OutputIndices property. The 'Auto' setting evaluates indices that
  %   puncture the second systematic stream and include all tail bits.
  OutputIndicesSource (1,:) char {matlab.system.mustBeMember(OutputIndicesSource,{'Auto','Property','Input port'})} = 'Auto';

  %OutputIndices Output indices
  %   Specify the bit ordering and puncturing used on the fully encoded
  %   data as a column vector of integers. The length of OutputIndices
  %   specifies the length of the encoder output. The default is
  %   getTurboIOIndices(64,2,3) to correspond to the default trellis and
  %   block length. Each element specified must be an integer value.
  OutputIndices = getTurboIOIndices(64,2,3);

end

properties (Constant, Hidden)
    InterleaverIndicesSourceSet = comm.CommonSets.getSet('SpecifyInputs');
end

properties(Access = private, Nontunable)
    % Constituent components
    cConvEnc1;      % encoder1
    cConvEnc2;      % encoder2
    % Commonly used attributes - set for default TrellisStructure
    pK = 1;             % number of uncoded bits
    pN = 2;             % number of coded bits
    pMLen = 3;          % memory length of a constituent encoder
    pNumTails = 6;      % number of tail bits per constituent encoder
end

methods
    % CONSTRUCTOR
    function obj = TurboEncoder(varargin)
        setProperties(obj, nargin, varargin{:}, 'TrellisStructure', ...
                      'InterleaverIndices');
    end
    
    function set.InterleaverIndices(obj, value)
            validateattributes(value, {'numeric'}, ...
            {'real', 'finite', 'positive', 'integer', 'vector'}, ...
            '', 'InterleaverIndices'); 

        obj.InterleaverIndices = value;
    end
    
    function set.OutputIndices(obj, value)
            validateattributes(value, {'numeric'}, ...
            {'real', 'finite', 'positive', 'integer', 'vector'}, ...
            '', 'OutputIndices'); 

        obj.OutputIndices = value;
    end
end

methods(Access = protected)   % System object APIs

    %% Number of Inputs
    function num = getNumInputsImpl(obj)
        num = 1 + strcmp(obj.InterleaverIndicesSource, 'Input port') ...
            + strcmp(obj.OutputIndicesSource, 'Input port');
    end

    %% Size propagators
    function flag = isInputSizeMutableImpl(obj,~)
        if strcmp(obj.InterleaverIndicesSource, 'Input port')
            flag = true;   % vars only via input port
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
        if strcmp(obj.OutputIndicesSource, 'Auto')
            inSize = propagatedInputSize(obj,1);
            pn = log2(obj.TrellisStructure.numOutputSymbols);
            pmlen = log2(obj.TrellisStructure.numStates);
            pnumtails = pmlen*pn;

            outLen = inSize(1) * (2*pn - 1) + 2*pnumtails;
        elseif strcmp(obj.OutputIndicesSource, 'Property')
            outLen = length(obj.OutputIndices);
        else
            if strcmp(obj.InterleaverIndicesSource, 'Input port')
                inSize = propagatedInputSize(obj,3);
            else
                inSize = propagatedInputSize(obj,2);
            end
            outLen = inSize(1);
        end
        varargout = {[outLen, 1]}; 
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
        coder.internal.errorIf(~isOk,'comm:TurboEncoder:InvalidTrellis', status);

        % Check for mixed modes
        coder.internal.errorIf( ...
           ( strcmp(obj.InterleaverIndicesSource,'Input port') && ...
             strcmp(obj.OutputIndicesSource,'Property') ) || ...
           ( strcmp(obj.InterleaverIndicesSource,'Property') && ...
             strcmp(obj.OutputIndicesSource,'Input port') ), ...
            'comm:TurboEncoder:InvalidMode');

        obj.cConvEnc1 = comm.ConvolutionalEncoder('TrellisStructure',...
            obj.TrellisStructure, 'TerminationMethod', 'Terminated');
        obj.cConvEnc2 = comm.ConvolutionalEncoder('TrellisStructure',...
            obj.TrellisStructure, 'TerminationMethod', 'Terminated');

        obj.pK = log2(obj.TrellisStructure.numInputSymbols);
        coder.internal.errorIf(obj.pK~=1,'comm:TurboEncoder:InvalidTrellis1ByN');

        obj.pN = log2(obj.TrellisStructure.numOutputSymbols);
        obj.pMLen = log2(obj.TrellisStructure.numStates);
        obj.pNumTails = obj.pMLen*(obj.pN);        
    end

    %% Step
    function y = stepImpl(obj, x, varargin)

        if strcmp(obj.InterleaverIndicesSource,'Input port')
            interlvrIndices = varargin{1};
            blkLen = length(interlvrIndices);
        else
            interlvrIndices = get(obj,'InterleaverIndices');
            blkLen = length(interlvrIndices);
        end
        coder.internal.assert(isequal((1:blkLen).',sort(interlvrIndices(:))),...
            'comm:TurboEncoder:InvalidIndices');
        coder.internal.assert(length(x) == blkLen,...
            'comm:TurboEncoder:InputLengthMismatch', blkLen);

        if strcmp(obj.OutputIndicesSource,'Input port')
            if strcmp(obj.InterleaverIndicesSource,'Input port')
                outIndices = varargin{2};
            else
                outIndices = varargin{1};
            end
        elseif strcmp(obj.OutputIndicesSource,'Property')
            outIndices = get(obj,'OutputIndices');
        else % Auto
            outIndices = getTurboIOIndices(blkLen,obj.pN,obj.pMLen,'LTE');
        end
        % Check outIndices
        fullOutLen = (blkLen+obj.pMLen)*obj.pN*2;
        coder.internal.assert(length(outIndices)>=blkLen,...
            'comm:TurboEncoder:InvalidOutIndicesLength', blkLen);
        
        if (blkLen > 0)            
            coder.internal.assert(max(outIndices(:))<=fullOutLen,...
                'comm:TurboEncoder:InvalidOutIndices', fullOutLen); 
             
            y1 = step(obj.cConvEnc1,x);
            y2 = step(obj.cConvEnc2,x(interlvrIndices));
            
            % Default interlaced full output
            y1r = reshape(y1,obj.pN,blkLen+obj.pMLen);
            y2r = reshape(y2,obj.pN,blkLen+obj.pMLen);
            ycomb = [y1r;y2r];
            allOut = ycomb(:);

            % Apply bit reordering and puncturing
            y = allOut(outIndices,1);        
        else
            % allow empty input
            y =  cast(zeros(0,1),'like',x); 
        end 
    end
    
    %%
    function releaseImpl(obj)
        release(obj.cConvEnc1);
        release(obj.cConvEnc2);
    end
    
    function flag = isInactivePropertyImpl(obj,prop)
        flag = false;
        switch prop
            case 'InterleaverIndices'
                flag = strcmp(obj.InterleaverIndicesSource,'Input port');
            case 'OutputIndices'
                if ( strcmp(obj.OutputIndicesSource,'Input port') || ...
                     strcmp(obj.OutputIndicesSource,'Auto') )    
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
            s.cConvEnc1  = matlab.System.saveObject(obj.cConvEnc1);
            s.cConvEnc2  = matlab.System.saveObject(obj.cConvEnc2);
        end
    end
  
    function loadObjectImpl(obj, s, wasLocked)
        if wasLocked
            obj.pK        = s.pK;
            obj.pN        = s.pN;
            obj.pMLen     = s.pMLen;
            obj.pNumTails = s.pNumTails;
            obj.cConvEnc1  = matlab.System.loadObject(s.cConvEnc1);
            obj.cConvEnc2  = matlab.System.loadObject(s.cConvEnc2);
        end
        % Call the base class method
        loadObjectImpl@matlab.System(obj, s);
    end

    function icon = getIconImpl(~)        
        icon = sprintf('Turbo\nEncoder');
    end
  
    function varargout = getInputNamesImpl(obj)
        varargout = cell(1, getNumInputs(obj));        
        if strcmp(obj.InterleaverIndicesSource,'Input port')
            varargout{1} = 'In';
            varargout{2} = 'IntrInd';
            if strcmp(obj.OutputIndicesSource,'Input port')
                varargout{3} = 'OutInd';
            end
        else
            if strcmp(obj.OutputIndicesSource,'Input port')
                varargout{1} = 'In';
                varargout{2} = 'OutInd';
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

methods(Static, Access=protected)
% System Block customization functions
  function header = getHeaderImpl
        header = matlab.system.display.Header(...
            'comm.TurboEncoder', ...
            'Title','comm:TurboEncoder:TurboEncoderTitle', ...
            'Text','comm:TurboEncoder:TurboEncoderDescription' );
  end    
    % Specify aliases to existing block parameters
    function groups = getPropertyGroupsImpl

        propsList = {'TrellisStructure', 'InterleaverIndicesSource', 'InterleaverIndices', ...
                              'OutputIndicesSource','OutputIndices'};
        paramsList = cell(1,numel(propsList));
            for ind = 1:numel(propsList)
                msgID = ['comm:TurboEncoder:',propsList{ind}];
                paramsList{ind} = matlab.system.display.internal.Property(...
                    propsList{ind},'Description',msgID);
            end
        paramsList{1} = matlab.system.display.internal.Property(propsList{1},'Description',['comm:TurboEncoder:',propsList{1}],'Alias', 'trellis'); 
        paramsList{3} = matlab.system.display.internal.Property(propsList{3},'Description',['comm:TurboEncoder:',propsList{3}],'Alias', 'intrIndices'); 
        
        params = matlab.system.display.Section(...
                'Title', 'comm:TurboEncoder:Parameters', ...
                'PropertyList',paramsList);


        groups = params;
    end

end

end
