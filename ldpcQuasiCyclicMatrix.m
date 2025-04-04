function H = ldpcQuasiCyclicMatrix(blockSize, P)
%LDPCQUASICYCLICMATRIX Parity-check matrix of a quasi-cyclic LDPC code
%
%   H = LDPCQUASICYCLICMATRIX(BLOCKSIZE, P) builds the parity-check matrix
%   H of a quasi-cyclic LDPC code using BLOCKSIZE and a prototype matrix P.
%   H is a sparse logical matrix.
%
%   BLOCKSIZE must be a positive integer. Each element of P is expanded to
%   a BLOCKSIZE-by-BLOCKSIZE submatrix in H. Each submatrix is either a
%   zero matrix or a cyclically-shifted version of a diagonal matrix. All
%   values of P must be -1, 0, or positive integers less than BLOCKSIZE. -1
%   produces a zero BLOCKSIZE-by-BLOCKSIZE submatrix. Other values indicate
%   the number of columns a BLOCKSIZE-by-BLOCKSIZE diagonal matrix should
%   be cyclically-shifted to the right. The number of columns in P must be
%   greater than the number of rows in P.
%
%   % Example 1:
%   blockSize = 3;
%   P = [0 -1 1 2; 2 1 -1 0];
%   H = ldpcQuasiCyclicMatrix(blockSize,P);
%
%   % Check that H is a sparse logical matrix
%   issparse(H) & islogical(H)
%
%   % Display H as a full matrix. Do so only if H is small
%   full(H)
%
%   % Example 2:
%   % Parity-check matrix of a rate 3/4 LDPC code from the WLAN standard
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
%
%   See also LDPCENCODERCONFIG, LDPCDECODERCONFIG, LDPCENCODE, LDPCDECODE.

% Copyright 2021 The MathWorks, Inc.

%#codegen

validateattributes(blockSize,{'numeric'},{'real','finite','positive','scalar','integer'},'ldpcQuasiCyclicMatrix','blockSize');
validateattributes(P,{'numeric'},{'real','2d','integer','>=',-1,'<',blockSize},'ldpcQuasiCyclicMatrix','P');
coder.internal.errorIf(size(P,1) >= size(P,2), 'comm:validateLDPCParityCheckMatrix:TooFewColumns');

% Each number in P not equal to -1 will produce blockSize ones in H
n = numel(find(P~=-1))*blockSize;
rowIndex = coder.nullcopy(zeros(n,1));
columnIndex = coder.nullcopy(zeros(n,1));

% Expand each number in P into a sub-matrix (blockSize by blockSize)
ind = 0;
[numRows, numCols] = size(P);
for j = 1:numCols
    for i = 1:numRows
        if P(i,j) ~= -1
            % Right-shift a blockSize-by-blockSize diagonal matrix
            % cyclically by P(i,j) columns 
            columnIndex(ind+(1:blockSize)) = (j-1)*blockSize + (1:blockSize);
            rowIndex(ind+(1:blockSize)) = (i-1)*blockSize + [(blockSize-P(i,j)+1):blockSize, 1:(blockSize-P(i,j))];
            ind = ind + blockSize;
        end
    end
end
H = sparse(rowIndex,columnIndex,true,numRows*blockSize,numCols*blockSize);