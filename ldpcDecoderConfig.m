classdef ldpcDecoderConfig < comm.internal.ConfigBase
%LDPCDECODERCONFIG LDPC decoder configuration
%
%   DECCFG = LDPCDECODERCONFIG returns the default LDPC decoder
%   configuration object based on a rate 5/6 LDPC code from the WLAN
%   standard and the default belief propagation decoding algorithm.
%
%   DECCFG = LDPCDECODERCONFIG(H) returns the LDPC decoder
%   configuration object based on a parity-check matrix H and the default
%   belief propagation decoding algorithm. H must be a sparse logical
%   matrix. The number of columns in H must be greater than the number of
%   rows in H.
%
%   DECCFG = LDPCDECODERCONFIG(ENCCFG) returns the LDPC decoder
%   configuration object based on ENCCFG, an LDPC encoder configuration
%   object of type ldpcEncoderConfig.
%
%   DECCFG = LDPCDECODERCONFIG(H, ALG) returns the LDPC decoder
%   configuration object based on a parity-check matrix H and the decoding
%   algorithm specified by a string ALG. ALG can be 'bp', 'layered-bp',
%   'norm-min-sum', or 'offset-min-sum' and the corresponding algorithms
%   are belief propagation decoding, layered belief propagation decoding,
%   normalized min-sum decoding, and offset min-sum decoding respectively.
%   If ALG is 'layered-bp', 'norm-min-sum', or 'offset-min-sum', the
%   parity-check matrix H should correspond to a quasi-cyclic LDPC code.
%
%   DECCFG = LDPCDECODERCONFIG(ENCCFG, ALG) returns the LDPC
%   decoder configuration object based on an LDPC encoder configuration
%   object ENCCFG and the decoding algorithm specified by a string ALG.
%   
%   LDPCDECODERCONFIG properties:
%  
%   ParityCheckMatrix  - Parity-check matrix of LDPC code
%   Algorithm          - LDPC decoding algorithm
%   BlockLength        - Block length of LDPC code
%   NumInformationBits - Number of information bits in an LDPC codeword
%   NumParityCheckBits - Number of parity-check bits in an LDPC codeword
%   CodeRate           - Code rate of LDPC code
%   NumRowsPerLayer    - Number of rows per layer, when Algorithm is
%                        'layered-bp', 'norm-min-sum', or 'offset-min-sum'
%
%   % Example:
%
%   % Create parity-check matrix of a rate 3/4 LDPC code from the WLAN standard
%   P = [
%        16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
%        25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
%        25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
%         9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
%        24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
%         2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
%       ];
%   blockSize = 27;
%   H = ldpcQuasiCyclicMatrix(blockSize, P);
%   % Create and display LDPC encoder configuration object
%   encoderCfg = ldpcEncoderConfig(H)
%   % Create and display LDPC decoder configuration object
%   decoderCfg = ldpcDecoderConfig(H)
%   % Transmit an LDPC-encoded, QPSK-modulated bit stream through an
%   % AWGN channel, then demodulate, decode, and count errors
%   awgnChan = comm.AWGNChannel(...
%           'NoiseMethod','Signal to noise ratio (SNR)','SNR',6);
%   ber = comm.ErrorRate;
%   for counter = 1:10
%     data           = randi([0 1],encoderCfg.NumInformationBits,1,'int8');
%     encodedData    = ldpcEncode(data, encoderCfg);
%     modSignal      = pskmod(encodedData,4,pi/8,InputType='Bit');
%     receivedSignal = awgnChan(modSignal);
%     demodSignal    = pskdemod(receivedSignal,4,pi/8,OutputType='approxllr',NoiseVariance=1/10^(awgnChan.SNR/10));
%     receivedBits   = ldpcDecode(demodSignal, decoderCfg, 10); % 10 iterations
%     errorStats     = ber(data, receivedBits);
%   end
%   fprintf('Error rate       = %1.2f\nNumber of errors = %d\n', ...
%     errorStats(1), errorStats(2))
%   
%   See also LDPCDECODE, LDPCENCODERCONFIG, LDPCQUASICYCLICMATRIX.

% Copyright 2021-2022 The MathWorks, Inc.

%#codegen

    properties
        %ParityCheckMatrix Parity-check matrix
        %   Specify the parity-check matrix as a sparse logical matrix with
        %   dimension (N-K) by N, where N > K > 0. The default is the
        %   parity-check matrix of a rate 5/6 LDPC code from the WLAN
        %   standard, which is the result of H in the following:
        %   P = [
        %        17 13  8 21  9  3 18 12 10  0  4 15 19  2  5 10 26 19 13 13  1  0 -1 -1
        %         3 12 11 14 11 25  5 18  0  9  2 26 26 10 24  7 14 20  4  2 -1  0  0 -1
        %        22 16  4  3 10 21 12  5 21 14 19  5 -1  8  5 18 11  5  5 15  0 -1  0  0
        %         7  7 14 14  4 16 16 24 24 10  1  7 15  6 10 26  8 18 21 14  1 -1 -1  0
        %       ];
        %   blockSize = 27;
        %   H = ldpcQuasiCyclicMatrix(blockSize, P);
        ParityCheckMatrix
        %Algorithm LDPC decoding algorithm
        %   Specify the LDPC decoding algorithm as one of 'bp' |
        %   'layered-bp' | 'norm-min-sum' | 'offset-min-sum'. The default
        %   is 'bp'.
        Algorithm = 'bp'
    end
    
    properties(Constant, Hidden)
        Algorithm_Values = {'bp','layered-bp','norm-min-sum','offset-min-sum'};
    end

    properties (SetAccess = private)
        %BlockLength Block length
        %   This read-only property indicates the block length of the LDPC
        %   code, i.e. the number of columns in the parity-check matrix.
        BlockLength
        %NumInformationBits Number of information bits
        %   This read-only property indicates the number of information
        %   bits in an LDPC codeword, i.e. the number of columns of the
        %   parity-check matrix minus the number of rows of the
        %   parity-check matrix.
        NumInformationBits
        %NumParityCheckBits Number of parity-check bits
        %   This read-only property indicates the number of parity-check
        %   bits in an LDPC codeword, i.e. the number of rows of the
        %   parity-check matrix.
        NumParityCheckBits
        %CodeRate Code rate
        %   This read-only property indicates the code rate of LDPC code,
        %   i.e. the ratio of NumInformationBits to BlockLength.
        CodeRate
        %NumRowsPerLayer Number of rows per layer
        %   This read-only property is visible when Algorithm is
        %   'layered-bp', 'norm-min-sum', or 'offset-min-sum'. This
        %   property indicates the number of rows per layer when using a
        %   layered decoding algorithm, i.e. the largest integer such that
        %   ParityCheckMatrix can be evenly split into consecutive
        %   submatrices
        %   ParityCheckMatrix(1:NumRowsPerLayer, :)
        %   ParityCheckMatrix((NumRowsPerLayer + 1):2*NumRowsPerLayer, :)
        %   ParityCheckMatrix((2*NumRowsPerLayer + 1):3*NumRowsPerLayer, :)
        %   ...
        %   ParityCheckMatrix((end - NumRowsPerLayer + 1):end, :)
        %   in which there is at most one '1' in any column in any one of
        %   these submatrices.
        NumRowsPerLayer = 0
    end

    properties (Hidden, SetAccess = private)
        derivedParams
        AlgorithmChoice = 0
    end

    methods
        function obj = ldpcDecoderConfig(varargin)
            setAlg = false;
            narginchk(0,2);
            if nargin == 0
                H = ldpcQuasiCyclicMatrix(27, ...
                      [
                       17 13  8 21  9  3 18 12 10  0  4 15 19  2  5 10 26 19 13 13  1  0 -1 -1
                        3 12 11 14 11 25  5 18  0  9  2 26 26 10 24  7 14 20  4  2 -1  0  0 -1
                       22 16  4  3 10 21 12  5 21 14 19  5 -1  8  5 18 11  5  5 15  0 -1  0  0
                        7  7 14 14  4 16 16 24 24 10  1  7 15  6 10 26  8 18 21 14  1 -1 -1  0
                      ]);
            elseif nargin == 1
                H = varargin{1};
            elseif nargin == 2
                H = varargin{1};
                setAlg = true;
            end
            if issparse(H) && islogical(H)
                obj.ParityCheckMatrix = H;
            elseif isa(H, 'ldpcEncoderConfig')
                obj.ParityCheckMatrix = H.ParityCheckMatrix;
            else
                error(message('comm:ldpc:NeedSparseLogicalMatrixOrConfigObj','ldpcEncoderConfig'));
            end
            if setAlg
                obj.Algorithm = varargin{2};                
            end
        end
        
        function obj = set.Algorithm(obj, alg)
            validateattributes(alg,{'char' 'string'},{},'','Algorithm');
            alg = validatestring(alg,{'bp','layered-bp','norm-min-sum','offset-min-sum'},'','Algorithm'); % Ensure that alg is one of the four acceptable strings after partial match
            switch alg
                case 'bp'
                    obj.AlgorithmChoice = 0; %#ok<*MCSUP> 
                case 'layered-bp'
                    obj.AlgorithmChoice = 1;
                case 'norm-min-sum'
                    obj.AlgorithmChoice = 2;
                case 'offset-min-sum'
                    obj.AlgorithmChoice = 3;
                otherwise
                    % Will never happen
                    obj.AlgorithmChoice = 4;
            end
            obj.Algorithm = alg;
        end
        
        function obj = set.ParityCheckMatrix(obj, ParityCheckMatrix)
            if issparse(ParityCheckMatrix) && islogical(ParityCheckMatrix)
                obj.ParityCheckMatrix = ParityCheckMatrix;
            else
                error(message('comm:ldpc:NeedSparseLogicalMatrix'));
            end
            coder.internal.errorIf(full(min(sum(ParityCheckMatrix,1),[],2))==0,...
                'comm:validateLDPCParityCheckMatrix:EmptyColumn');
            coder.internal.errorIf(full(min(sum(ParityCheckMatrix,2),[],1))==0,...
                'comm:validateLDPCParityCheckMatrix:EmptyRow');
            obj.BlockLength = size(ParityCheckMatrix,2);
            obj.NumParityCheckBits = size(ParityCheckMatrix,1);
            obj.NumInformationBits = obj.BlockLength - obj.NumParityCheckBits;
            coder.internal.errorIf(obj.NumInformationBits <= 0, 'comm:validateLDPCParityCheckMatrix:TooFewColumns');
            obj.CodeRate = obj.NumInformationBits/obj.BlockLength;
            obj.derivedParams = CalcDerivedParams(obj);
            obj.NumRowsPerLayer = CalcNumRowsPerLayer(obj);
        end
    end
    
    methods (Access = private)
        function derivedParams = CalcDerivedParams(obj)
            [~, nCols] = size(obj.ParityCheckMatrix);
            if size(obj.ParityCheckMatrix,1) == 1
                [columnIndex, rowIndex] = find(obj.ParityCheckMatrix');
            else
                [rowIndex, columnIndex] = find(obj.ParityCheckMatrix);
            end
            temp = sortrows([(rowIndex-1)*nCols+columnIndex (1:nnz(obj.ParityCheckMatrix))']);
            indexMap = temp(:,2);

            [columnIndex, ~] = find(obj.ParityCheckMatrix');

            rowWeight = full(sum(obj.ParityCheckMatrix,2));
            rowOffset = [0; cumsum(rowWeight(1:end-1,1))];

            columnWeight = full(sum(obj.ParityCheckMatrix,1))';
            columnOffset = [0; cumsum(columnWeight(1:end-1,1))];

            derivedParams.offsetWeight = int32([rowOffset; rowWeight; columnOffset; columnWeight]);
            derivedParams.columnIndexMap = int32([columnIndex; indexMap] - 1);
        end
    
        function numRowsPerLayer = CalcNumRowsPerLayer(obj)
            t = obj.ParityCheckMatrix'; % More efficient to use columns
            totalsum = zeros(size(t,1),1);
            numRowsPerLayer = 0;
            for i = 1:size(t,2)
                totalsum = totalsum + t(:,i);
                if max(totalsum) > 1
                    % Assuming that the prototype matrix of the QC LDPC
                    % code has one value not equal to -1 in the first row,
                    % when totalsum > 1, i is in the beginning of the
                    % second partition
                    break;
                end
                numRowsPerLayer = numRowsPerLayer + 1;
            end
            if mod(size(obj.ParityCheckMatrix,1),numRowsPerLayer) ~= 0
                % Not an integer number of partitions
                % Be conservative and set numRowsPerLayer = 1
                numRowsPerLayer = 1;
            end
        end
    end
    
    methods (Access = protected)
        function flag = isInactiveProperty(obj, prop)
            flag = false;
            if strcmp(prop,'NumRowsPerLayer')
                flag = strcmpi(obj.Algorithm,'bp');
            end
        end
    end
end