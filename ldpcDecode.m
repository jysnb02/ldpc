function [decoderOut, actualNumIter, finalParityChecks] = ldpcDecode(llr, decoderConfig, maxNumIter, varargin)
%LDPCDECODE Decode binary low-density parity-check code
%
%   [OUT, ACTUALNUMITER, FINALPARITYCHECKS] = LDPCDECODE(LLR, DECCFG,
%   MAXNUMITER) decodes log-likelihood ratios LLR by using DECCFG, an
%   LDPC decoder configuration object of type ldpcDecoderConfig, and
%   at most MAXNUMITER iterations. LLR must be a double or single matrix with
%   DECCFG.BlockLength rows. A positive log-likelihood ratio indicates
%   that the corresponding bit is more likely a zero. OUT is an int8 matrix
%   with DECCFG.NumInformationBits rows, representing the decoded bits
%   for LLR(1:DECCFG.NumInformationBits,:). Each column of LLR
%   corresponds to a codeword and is decoded independently. Decoding may stop
%   before MAXNUMITER iterations, if all parity-checks for a codeword are
%   satisfied. ACTUALNUMITER is a row vector of the actual number of
%   iterations executed for the codewords. FINALPARITYCHECKS has
%   DECCFG.NumParityCheckBits rows and each column is the final
%   parity-checks for the corresponding codeword.
%
%   [OUT, ACTUALNUMITER, FINALPARITYCHECKS] = LDPCDECODE(...,Name,Value)
%   specifies additional name-value pair arguments described below:
%
%   'OutputFormat'        - One of 'info', 'whole', specifies the output format.
%                           OUT contains decoded information bits (default) or
%                           whole LDPC codeword bits. For 'info', the number of
%                           rows in OUT is the length of the information bits and
%                           for 'whole', the number of rows in OUT is the codeword
%                           length.
%   'DecisionType'        - One of 'hard', 'soft', specifies the decision type
%                           used for decoding. For 'hard' (default), output is
%                           decoded bits of 'int8' type. For 'soft', output is
%                           log-likelihood ratios with the same type as input.
%   'MinSumScalingFactor' - Specifies the scaling factor for normalized min-sum
%                           decoding algorithm as a real scalar greater than 0 and
%                           less than or equal to 1. The default is 0.75. The
%                           value is only applicable when DECCFG.Algorithm
%                           is set to 'norm-min-sum'.
%   'MinSumOffset'        - Specifies the offset for offset min-sum decoding 
%                           algorithm as a finite real scalar greater than or 
%                           equal to 0. The default is 0.5. The value is only 
%                           applicable when DECCFG.Algorithm is set to
%                           'offset-min-sum'.
%   'Termination'         - One of 'early', 'max', specifies the decoding
%                           termination criteria. For 'early' (default), decoding
%                           is terminated when all parity-checks are satisfied, up
%                           to a maximum number of iterations given by MAXNUMITER.
%                           For 'max', decoding continues till MAXNUMITER
%                           iterations are completed.
%   'Multithreaded'       - True (default) or false. If in MATLAB interpreted mode
%                           and this value is true, the decoding algorithm will be
%                           executed with multiple threads on the CPU. This will
%                           speed up LDPC decoding if the parity-check matrix is big.
%                           This name-value pair argument is ignored if LLR is a
%                           gpuArray.
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
%   See also LDPCENCODE, LDPCDECODERCONFIG, LDPCENCODERCONFIG, LDPCQUASICYCLICMATRIX.

% Copyright 2021-2024 The MathWorks, Inc.

%#codegen

if ~isa(decoderConfig,'ldpcDecoderConfig') || ~isscalar(decoderConfig)
    error(message('MATLAB:validateattributes:expected','the second input argument','an ldpcDecoderConfig object'));
end

validateattributes(llr, {'double','single'}, {'real','nonnan','nrows',decoderConfig.BlockLength}, 1);
% Protects against all non-native double and single types
comm.internal.utilities.isClassDoubleSingleGpuArray(llr,1,"LLR")

validateattributes(maxNumIter, {'numeric'}, {'real','scalar','positive','integer'}, 3);

% Parse
defaults = struct('OutputFormat', 'info', 'DecisionType', 'hard', ...
                  'MinSumScalingFactor', 0.75, 'MinSumOffset', 0.5, ...
                  'Termination', 'early', 'Multithreaded', true);
res = comm.internal.utilities.nvParser(defaults, varargin{:});
validatestring(res.OutputFormat, {'whole','info'}, mfilename, 'OutputFormat');
validatestring(res.DecisionType, {'hard','soft'}, mfilename, 'DecisionType');
validatestring(res.Termination, {'early','max'}, mfilename, 'Termination');
validateattributes(res.Multithreaded,{'logical'},{'scalar'},mfilename,'Multithreaded');

isWholeCodeword = strncmpi(res.OutputFormat, 'w', 1);
isHardDecision = int8(strncmpi(res.DecisionType, 'h', 1));
isEarlyExit = int8(strncmpi(res.Termination, 'e', 1));
useMultithread = int8(res.Multithreaded);

if decoderConfig.AlgorithmChoice == 2
    alphaBeta = res.MinSumScalingFactor;
    validateattributes(alphaBeta,{'numeric'},{'scalar','real','>',0,'<=',1},'ldpcDecode','MinSumScalingFactor');
elseif decoderConfig.AlgorithmChoice == 3
    alphaBeta = res.MinSumOffset;
    validateattributes(alphaBeta,{'numeric'},{'scalar','real','finite','>=',0},'ldpcDecode','MinSumOffset');
else
    alphaBeta = 0;
end

C = size(llr,2);
typeIn = class(llr);
finalParityChecks = coder.nullcopy(zeros(decoderConfig.NumParityCheckBits,C));
actualNumIter = coder.nullcopy(zeros(1,C));
compFinalParityChecks = int8(nargout == 3);
numEdges = length(decoderConfig.derivedParams.columnIndexMap)/2;

if coder.target('MATLAB')
    if isa(llr,'gpuArray')
        [decoderOut,actualNumIter,finalParityChecks] = comm.internal.ldpc.gpuDecode(llr,decoderConfig,maxNumIter,...
                                                           isWholeCodeword,isHardDecision,isEarlyExit,alphaBeta);
    else
        decoderOut = zeros(decoderConfig.BlockLength, C, 'like', llr);
        blockLen      = decoderConfig.BlockLength;
        nPBits        = decoderConfig.NumParityCheckBits;
        nRowsPerLayer = decoderConfig.NumRowsPerLayer;
        oWeight       = decoderConfig.derivedParams.offsetWeight;
        cIndexMap     = decoderConfig.derivedParams.columnIndexMap;
        algChoice     = decoderConfig.AlgorithmChoice;

        for cwIdx = 1:C
            [decoderOut(:,cwIdx),actualNumIter(cwIdx),finalParityChecks(:,cwIdx)] = ...
                mwcomm_ldpcdecode_mt(llr(:,cwIdx),double(gather(maxNumIter)), ...
                blockLen,nPBits,numEdges,nRowsPerLayer,oWeight,cIndexMap, ...
                isHardDecision,isEarlyExit,compFinalParityChecks,algChoice, ...
                cast(alphaBeta,'like',llr),useMultithread);
        end
        % Format output
        if ~isWholeCodeword
            decoderOut = decoderOut(1:decoderConfig.NumInformationBits,:);
        end
        if isHardDecision
            decoderOut = cast(decoderOut,'int8');
        else
            decoderOut = cast(decoderOut,typeIn);
        end
    end
else % codegen
    % Initialize output to fix dimension and data type
    if isHardDecision
        if isWholeCodeword
            decoderOut = coder.nullcopy(zeros(decoderConfig.BlockLength, C, 'int8'));
        else
            decoderOut = coder.nullcopy(zeros(decoderConfig.NumInformationBits, C, 'int8'));
        end
    else
        if isWholeCodeword
            decoderOut = coder.nullcopy(zeros(decoderConfig.BlockLength, C, typeIn));
        else
            decoderOut = coder.nullcopy(zeros(decoderConfig.NumInformationBits, C, typeIn));
        end
    end

    LLR_out = coder.nullcopy(zeros(decoderConfig.BlockLength,C));
    
    switch decoderConfig.AlgorithmChoice
        case 0 % Use BP
            for cwIdx = 1:C
                [LLR_out(:,cwIdx), actualNumIter(:,cwIdx), finalParityChecks(:,cwIdx)] = comm.internal.ldpc.BPDecode(llr(:,cwIdx), ...
                    maxNumIter, decoderConfig.BlockLength, decoderConfig.NumParityCheckBits, numEdges, ...
                    decoderConfig.derivedParams.offsetWeight, decoderConfig.derivedParams.columnIndexMap, isEarlyExit, compFinalParityChecks);
            end
        case 1 % Use Layered BP
            for cwIdx = 1:C
                [LLR_out(:,cwIdx), actualNumIter(:,cwIdx), finalParityChecks(:,cwIdx)] = comm.internal.ldpc.LayeredBPDecode(llr(:,cwIdx), ...
                    maxNumIter, decoderConfig.NumParityCheckBits, numEdges, ...
                    decoderConfig.derivedParams.offsetWeight, decoderConfig.derivedParams.columnIndexMap, isEarlyExit, compFinalParityChecks);
            end
        case 2 % Use Layered normalized min-sum
            for cwIdx = 1:C
                [LLR_out(:,cwIdx), actualNumIter(:,cwIdx), finalParityChecks(:,cwIdx)] = comm.internal.ldpc.LayeredBPNormMSDecode(llr(:,cwIdx), ...
                    maxNumIter, decoderConfig.NumParityCheckBits, numEdges, ...
                    decoderConfig.derivedParams.offsetWeight, decoderConfig.derivedParams.columnIndexMap, isEarlyExit, compFinalParityChecks, alphaBeta);
            end
        otherwise % 3
            % Use Layered offset min-sum
            for cwIdx = 1:C
                [LLR_out(:,cwIdx), actualNumIter(:,cwIdx), finalParityChecks(:,cwIdx)] = comm.internal.ldpc.LayeredBPOffsetMSDecode(llr(:,cwIdx), ...
                    maxNumIter, decoderConfig.NumParityCheckBits, numEdges, ...
                    decoderConfig.derivedParams.offsetWeight, decoderConfig.derivedParams.columnIndexMap, isEarlyExit, compFinalParityChecks, alphaBeta);
            end
    end
    
    % Format output
    if isHardDecision
        if isWholeCodeword
            decoderOut(:) = int8(LLR_out <= 0);
        else
            decoderOut(:) = int8(LLR_out(1:decoderConfig.NumInformationBits,:) <= 0);
        end
    else
        if isWholeCodeword
            decoderOut(:) = cast(LLR_out, typeIn);
        else
            decoderOut(:) = cast(LLR_out(1:decoderConfig.NumInformationBits,:), typeIn);
        end
    end
end