function output = ldpcEncode(informationBits, encoderConfig, varargin)
%LDPCENCODE Encode binary low-density parity-check code
%
%   OUTPUT = LDPCENCODE(INFORMATIONBITS, ENCCFG) encodes INFORMATIONBITS
%   by the LDPC code specified by ENCCFG, an LDPC encoder configuration
%   object of type ldpcEncoderConfig. OUTPUT is a systematic codeword,
%   i.e. OUTPUT(1:ENCCFG.NumInformationBits,:) is equal to INFORMATIONBITS
%   and OUTPUT((ENCCFG.NumInformationBits+1):end,:) is the parity-check
%   bits. OUTPUT has the same data type as INFORMATIONBITS. INFORMATIONBITS must
%   be a double, single, int8, or logical matrix with ENCCFG.NumInformationBits
%   rows. Each column of INFORMATIONBITS is encoded independently. Nonzero values
%   in INFORMATIONBITS are treated as ones.
%
%   OUTPUT = LDPCENCODE(INFORMATIONBITS, ENCCFG, 'OutputFormat',
%   FMT) encodes INFORMATIONBITS and produces OUTPUT according to FMT. FMT
%   can be 'whole' or 'parity'. If FMT is 'whole', OUTPUT is the whole
%   codeword including the information bits and the parity-check bits at
%   the end. If FMT is 'parity', OUTPUT has ENCCFG.NumParityCheckBits
%   rows and contains the parity-check bits only.
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
%   % Create LDPC encoder configuration object
%   encoderCfg = ldpcEncoderConfig(H);
%   % Generate random information bits
%   infoBits = rand(encoderCfg.NumInformationBits,100) < 0.5;
%   % Encode information bits by LDPC code
%   codeword = ldpcEncode(infoBits, encoderCfg);
%
%   See also LDPCDECODE, LDPCENCODERCONFIG, LDPCDECODERCONFIG, LDPCQUASICYCLICMATRIX.

% Copyright 2021-2024 The MathWorks, Inc.

%#codegen

% Parse
defaults = struct('OutputFormat', 'whole');
res = comm.internal.utilities.nvParser(defaults, varargin{:});
fName = 'ldpcEncode';
validatestring(res.OutputFormat, {'whole','parity'}, fName, 'OutputFormat');
isWholeCodeword = strncmpi(res.OutputFormat, 'w', 1);

if ~isa(encoderConfig,'ldpcEncoderConfig') || ~isscalar(encoderConfig)
    error(message('MATLAB:validateattributes:expected','the second input argument','an ldpcEncoderConfig object'));
end
validateattributes(informationBits, {'double','single','int8','logical'}, {'2d','real','finite','nrows',encoderConfig.NumInformationBits}, 1);
comm.internal.utilities.isClassDoubleSingleInt8Logical(informationBits,1,'INFORMATIONBITS')

ParityCheckBits = coder.nullcopy(zeros(encoderConfig.NumParityCheckBits, size(informationBits,2), 'like', informationBits));
if coder.target('MATLAB')
    for i = 1:size(informationBits,2)
        ParityCheckBits(:,i) = mwcomm_ldpcencode(int8(informationBits(:,i)), encoderConfig.NumInformationBits, encoderConfig.NumParityCheckBits, ...
            encoderConfig.derivedParams.EncodingMethod, ...
            encoderConfig.derivedParams.MatrixA_RowIndices, encoderConfig.derivedParams.MatrixA_RowStartLoc, encoderConfig.derivedParams.MatrixA_ColumnSum, ...
            encoderConfig.derivedParams.MatrixB_RowIndices, encoderConfig.derivedParams.MatrixB_RowStartLoc, encoderConfig.derivedParams.MatrixB_ColumnSum, ...
            encoderConfig.derivedParams.MatrixL_RowIndices, encoderConfig.derivedParams.MatrixL_RowStartLoc, encoderConfig.derivedParams.MatrixL_ColumnSum, ...
            encoderConfig.derivedParams.RowOrder)';
    end
else
    for i = 1:size(informationBits,2)
        if encoderConfig.derivedParams.RowOrder(1) >= 0
            % RowOrder[0] >= 0 only if the last (N-K) columns of H are not triangular, or if they are
            % lower/upper triangular along the anti-diagonal and has a full anti-diagonal, e.g.,
            %
            % [ 0 0 1             [ 1 0 1
            %   0 1 0        or     1 1 0      or    a non-triangular matrix
            %   1 1 1 ]             1 0 0 ]

            % Compute the matrix product between first K columns of H and the information bits
            MatrixProductBuffer = GF2MatrixMul(informationBits(:,i), encoderConfig.NumParityCheckBits, encoderConfig.derivedParams.MatrixA_RowIndices, encoderConfig.derivedParams.MatrixA_RowStartLoc, encoderConfig.derivedParams.MatrixA_ColumnSum);
            % Need to perform this substitution if the last (N-K) columns of H are not triangular
            if encoderConfig.derivedParams.EncodingMethod == 0
                % Forward substitution for the lower triangular matrix obtained from factorization in GF(2)
                MatrixProductBuffer = GF2Subst(MatrixProductBuffer, encoderConfig.derivedParams.MatrixL_RowIndices, encoderConfig.derivedParams.MatrixL_RowStartLoc, encoderConfig.derivedParams.MatrixL_ColumnSum, 1);
            end
            % In this case, MatrixProductBuffer = ReorderBuffer
            % Do an extra re-ordering step
            ParityCheckBits(:,i) = MatrixProductBuffer(encoderConfig.derivedParams.RowOrder+1);

            % Solve for the parity-check bits
            %
            % Common step: If object property EncodingAlgorithm is 'Matrix Inverse', do backward substitution;
            % otherwise, do forward or backward substitution according to EncodingAlgorithm
            ParityCheckBits(:,i) = GF2Subst(ParityCheckBits(:,i), encoderConfig.derivedParams.MatrixB_RowIndices, encoderConfig.derivedParams.MatrixB_RowStartLoc, encoderConfig.derivedParams.MatrixB_ColumnSum, encoderConfig.derivedParams.EncodingMethod);

        else
            % RowOrder[0] < 0 only if the last (N-K) columns of H are lower/upper triangular and has
            % a full diagonal, e.g.,
            %
            % [ 1 0 0             [ 1 1 1
            %   1 1 0        or     0 1 0
            %   0 1 1 ]             0 0 1 ]

            % Compute the matrix product between first K columns of H and the information bits
            ParityCheckBits(:,i) = GF2MatrixMul(informationBits(:,i), encoderConfig.NumParityCheckBits, encoderConfig.derivedParams.MatrixA_RowIndices, encoderConfig.derivedParams.MatrixA_RowStartLoc, encoderConfig.derivedParams.MatrixA_ColumnSum);
            ParityCheckBits(:,i) = GF2Subst(ParityCheckBits(:,i), encoderConfig.derivedParams.MatrixB_RowIndices, encoderConfig.derivedParams.MatrixB_RowStartLoc, encoderConfig.derivedParams.MatrixB_ColumnSum, encoderConfig.derivedParams.EncodingMethod);
        end
    end
end

if isWholeCodeword
    output = [informationBits; ParityCheckBits];
else
    output = ParityCheckBits;
end


function dest = GF2MatrixMul(source, destLen, RowIndices, rowloc, ColumnSum)
% GF2MatrixMul computes the modulo-2 matrix product of a vector (specified by source) and
% a matrix (specified by RowIndices, rowloc and ColumnSum).
%
% RowIndices specifies the row indices of the nonzero elements in the matrix.
% ColumnSum specifies the number of nonzero elements in each column.
% rowloc specifies the offset (in RowIndices) of the first nonzero element in each column.
%
% E.g. for this matrix
%      [ 1     0     0
%        0     1     0
%        1     1     1 ]
%
% RowIndices = [0 2 1 2 2]
% rowloc     = [0 2 4]
% ColumnSum  = [2 2 1]

srclen = length(source);
dest = zeros(destLen,1);

% Start from the first column
for columnindex = 1:srclen
    if source(columnindex) ~= 0
        rowindex = RowIndices(rowloc(columnindex) + (1:ColumnSum(columnindex)))+1;
        dest(rowindex) = 1 - dest(rowindex);
    end
end

function srcdest = GF2Subst(srcdest, RowIndices, rowloc, ColumnSum, direction)
% GF2Subst solves a system of linear equations by forward/backward substitution (if source == dest).
%
% RowIndices specifies the row indices of the nonzero elements in the matrix.
% ColumnSum specifies the number of nonzero elements in each column.
% rowloc specifies the offset (in RowIndices) of the first nonzero element in each column.
%
% E.g. for this matrix
%      [ 1     0     0
%        0     1     0
%        1     1     1 ]
%
% RowIndices = [0 2 1 2 2]
% rowloc     = [0 2 4]
% ColumnSum  = [2 2 1]

srclen = length(srcdest);

if direction == 1   % direction = 1 or -1
    % Start from the first column for forward substitution
    for columnindex = 1:srclen
        if srcdest(columnindex) ~= 0
            rowindex = RowIndices(rowloc(columnindex) + (1:ColumnSum(columnindex)),1)+1;
            srcdest(rowindex) = 1 - srcdest(rowindex);
        end
    end
else
    % Start from the last column for backward substitution
    for columnindex = srclen:-1:1
        if srcdest(columnindex) ~= 0
            rowindex = RowIndices(rowloc(columnindex) + (1:ColumnSum(columnindex)),1)+1;
            srcdest(rowindex) = 1 - srcdest(rowindex);
        end
    end
end
