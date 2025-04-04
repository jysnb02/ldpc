classdef ldpcEncoderConfig < comm.internal.ConfigBase
%LDPCENCODERCONFIG LDPC encoder configuration
%
%   ENCCFG = LDPCENCODERCONFIG returns the default LDPC encoder
%   configuration object based on a rate 5/6 LDPC code from the WLAN
%   standard.
%
%   ENCCFG = LDPCENCODERCONFIG(H) returns the LDPC encoder
%   configuration object based on a parity-check matrix H. H must be a
%   sparse logical matrix. The number of columns in H must be greater than
%   the number of rows in H. If H has N columns and (N-K) rows, the last
%   (N-K) columns of H, H(:,(K+1):end), must be invertible in GF(2).
%
%   ENCCFG = LDPCENCODERCONFIG(DECCFG) returns the LDPC encoder
%   configuration object based on DECCFG, an LDPC decoder configuration
%   object of type ldpcDecoderConfig.
%   
%   LDPCENCODERCONFIG properties:
%  
%   ParityCheckMatrix  - Parity-check matrix of LDPC code
%   BlockLength        - Block length of LDPC code
%   NumInformationBits - Number of information bits in an LDPC codeword
%   NumParityCheckBits - Number of parity-check bits in an LDPC codeword
%   CodeRate           - Code rate of LDPC code
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
%   % Generate random information bits
%   infoBits = rand(encoderCfg.NumInformationBits,1) < 0.5;
%   % Encode information bits by LDPC code
%   codeword = ldpcEncode(infoBits, encoderCfg);
%
%   See also LDPCENCODE, LDPCDECODERCONFIG, LDPCQUASICYCLICMATRIX.

% Copyright 2021-2023 The MathWorks, Inc.

%#codegen

    properties
        %ParityCheckMatrix Parity-check matrix
        %   Specify the parity-check matrix as a sparse logical matrix with
        %   dimension (N-K) by N, where N > K > 0. The last (N-K) columns
        %   of the parity-check matrix must be invertible in GF(2). The
        %   default is the parity-check matrix of a rate 5/6 LDPC code from
        %   the WLAN standard, which is the result of H in the following:
        %   P = [
        %        17 13  8 21  9  3 18 12 10  0  4 15 19  2  5 10 26 19 13 13  1  0 -1 -1
        %         3 12 11 14 11 25  5 18  0  9  2 26 26 10 24  7 14 20  4  2 -1  0  0 -1
        %        22 16  4  3 10 21 12  5 21 14 19  5 -1  8  5 18 11  5  5 15  0 -1  0  0
        %         7  7 14 14  4 16 16 24 24 10  1  7 15  6 10 26  8 18 21 14  1 -1 -1  0
        %       ];
        %   blockSize = 27;
        %   H = ldpcQuasiCyclicMatrix(blockSize, P);
        ParityCheckMatrix
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
    end

    properties (Hidden, SetAccess = private)
        derivedParams
        gpuParams
    end

    methods
        function obj = ldpcEncoderConfig(varargin)
            narginchk(0,1);
            if nargin == 0
                H = ldpcQuasiCyclicMatrix(27, ...
                      [
                       17 13  8 21  9  3 18 12 10  0  4 15 19  2  5 10 26 19 13 13  1  0 -1 -1
                        3 12 11 14 11 25  5 18  0  9  2 26 26 10 24  7 14 20  4  2 -1  0  0 -1
                       22 16  4  3 10 21 12  5 21 14 19  5 -1  8  5 18 11  5  5 15  0 -1  0  0
                        7  7 14 14  4 16 16 24 24 10  1  7 15  6 10 26  8 18 21 14  1 -1 -1  0
                      ]);
            else
                H = varargin{1};
            end
            
            if issparse(H) && islogical(H)
                obj.ParityCheckMatrix = H;
            elseif isa(H, 'ldpcDecoderConfig')
                obj.ParityCheckMatrix = H.ParityCheckMatrix;
            else
                error(message('comm:ldpc:NeedSparseLogicalMatrixOrConfigObj','ldpcDecoderConfig'));
            end
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
            obj.ParityCheckMatrix = ParityCheckMatrix;
            obj.BlockLength = size(ParityCheckMatrix,2); %#ok<*MCSUP> 
            obj.NumParityCheckBits = size(ParityCheckMatrix,1);
            obj.NumInformationBits = obj.BlockLength - obj.NumParityCheckBits;
            coder.internal.errorIf(obj.NumInformationBits <= 0, 'comm:validateLDPCParityCheckMatrix:TooFewColumns');
            obj.CodeRate = obj.NumInformationBits/obj.BlockLength;
            obj.derivedParams = CalcDerivedParams(obj);
            if coder.target('MATLAB')
                if canUseGPU()
                    obj.gpuParams.MatrixL_RowIndices = gpuArray(obj.derivedParams.MatrixL_RowIndices);
                    obj.gpuParams.MatrixL_RowStartLoc = gpuArray(obj.derivedParams.MatrixL_RowStartLoc);
                    obj.gpuParams.MatrixL_ColumnSum = gpuArray(obj.derivedParams.MatrixL_ColumnSum);
                    obj.gpuParams.RowOrder = gpuArray(obj.derivedParams.RowOrder);
                    obj.gpuParams.MatrixA_RowIndices = gpuArray(obj.derivedParams.MatrixA_RowIndices);
                    obj.gpuParams.MatrixA_RowStartLoc = gpuArray(obj.derivedParams.MatrixA_RowStartLoc);
                    obj.gpuParams.MatrixA_ColumnSum = gpuArray(obj.derivedParams.MatrixA_ColumnSum);
                    obj.gpuParams.MatrixB_RowIndices = gpuArray(obj.derivedParams.MatrixB_RowIndices);
                    obj.gpuParams.MatrixB_RowStartLoc = gpuArray(obj.derivedParams.MatrixB_RowStartLoc);
                    obj.gpuParams.MatrixB_ColumnSum = gpuArray(obj.derivedParams.MatrixB_ColumnSum);
                end
            end
        end
    end
    
    methods (Access = private)
        function derivedParams = CalcDerivedParams(obj)
            N = size(obj.ParityCheckMatrix,2);
            K = N - size(obj.ParityCheckMatrix,1);

            PB = obj.ParityCheckMatrix(:,(K+1):end); % Extract the last (N-K) columns of the parity-check matrix.

            % Check if PB is triangular
            switch isfulldiagtriangular(PB)
                case 1
                    % Forward Substitution
                    % PB is lower triangular and has a full diagonal.
                    rowOrder = -1;                  % Don't need to reverse the order of rows in PB
                    derivedParams.EncodingMethod = int8(1);
                    P = tril(PB, -1);  % Remove the diagonal of PB and put it in P.
                    derivedParams.MatrixL_RowIndices = int32(0);
                    derivedParams.MatrixL_RowStartLoc = int32(0);
                    derivedParams.MatrixL_ColumnSum = int32(0);
                    derivedParams.RowOrder = int32(rowOrder-1);
                    [derivedParams.MatrixA_RowIndices, ...
                        derivedParams.MatrixA_RowStartLoc, ...
                        derivedParams.MatrixA_ColumnSum] = ConvertMatrixFormat(obj.ParityCheckMatrix(:,1:K));
                    [derivedParams.MatrixB_RowIndices, ...
                        derivedParams.MatrixB_RowStartLoc, ...
                        derivedParams.MatrixB_ColumnSum] = ConvertMatrixFormat(P);
                case -1
                    % Backward Substitution
                    % PB is upper triangular and has a full diagonal.
                    rowOrder = -1;                  % Don't need to reverse the order of rows in PB
                    derivedParams.EncodingMethod = int8(-1);
                    P = triu(PB, 1);   % Remove the diagonal of PB and put it in P.
                    derivedParams.MatrixL_RowIndices = int32(0);
                    derivedParams.MatrixL_RowStartLoc = int32(0);
                    derivedParams.MatrixL_ColumnSum = int32(0);
                    derivedParams.RowOrder = int32(rowOrder-1);
                    [derivedParams.MatrixA_RowIndices, ...
                        derivedParams.MatrixA_RowStartLoc, ...
                        derivedParams.MatrixA_ColumnSum] = ConvertMatrixFormat(obj.ParityCheckMatrix(:,1:K));
                    [derivedParams.MatrixB_RowIndices, ...
                        derivedParams.MatrixB_RowStartLoc, ...
                        derivedParams.MatrixB_ColumnSum] = ConvertMatrixFormat(P);
                otherwise % 0
                    % Reverse the order of rows in PB, but keep PB, since if PB is not
                    % triangular, we need to factorize it in GF(2).
                    Reversed_PB = PB(end:-1:1,:);
                    switch isfulldiagtriangular(Reversed_PB) % Check if Reversed_PB is triangular.
                        case 1
                            % Forward Substitution
                            rowOrder = ((N-K):-1:1)';
                            PB = Reversed_PB;
                            derivedParams.EncodingMethod = int8(1);
                            P = tril(PB, -1);  % Remove the diagonal of PB and put it in P.
                            derivedParams.MatrixL_RowIndices = int32(0);
                            derivedParams.MatrixL_RowStartLoc = int32(0);
                            derivedParams.MatrixL_ColumnSum = int32(0);
                            derivedParams.RowOrder = int32(rowOrder-1);
                            [derivedParams.MatrixA_RowIndices, ...
                                derivedParams.MatrixA_RowStartLoc, ...
                                derivedParams.MatrixA_ColumnSum] = ConvertMatrixFormat(obj.ParityCheckMatrix(:,1:K));
                            [derivedParams.MatrixB_RowIndices, ...
                                derivedParams.MatrixB_RowStartLoc, ...
                                derivedParams.MatrixB_ColumnSum] = ConvertMatrixFormat(P);
                        case -1
                            % Backward Substitution
                            rowOrder = ((N-K):-1:1)';
                            PB = Reversed_PB;
                            derivedParams.EncodingMethod = int8(-1);
                            P = triu(PB, 1);   % Remove the diagonal of PB and put it in P.
                            derivedParams.MatrixL_RowIndices = int32(0);
                            derivedParams.MatrixL_RowStartLoc = int32(0);
                            derivedParams.MatrixL_ColumnSum = int32(0);
                            derivedParams.RowOrder = int32(rowOrder-1);
                            [derivedParams.MatrixA_RowIndices, ...
                                derivedParams.MatrixA_RowStartLoc, ...
                                derivedParams.MatrixA_ColumnSum] = ConvertMatrixFormat(obj.ParityCheckMatrix(:,1:K));
                            [derivedParams.MatrixB_RowIndices, ...
                                derivedParams.MatrixB_RowStartLoc, ...
                                derivedParams.MatrixB_ColumnSum] = ConvertMatrixFormat(P);
                        otherwise % 0
                            derivedParams.EncodingMethod = int8(0);
                            % If we reach this line, then the last (N-K) columns of the
                            % parity-check matrix (i.e. PB) is not triangular. Let's
                            % try to factorize it in GF(2).
                            [PL, PC, rowOrder, invertible] = gf2factorize(PB);
                            if ~invertible
                                derivedParams.MatrixL_RowIndices = int32(0);
                                derivedParams.MatrixL_RowStartLoc = int32(0);
                                derivedParams.MatrixL_ColumnSum = int32(0);
                                derivedParams.RowOrder = int32(-2);
                                derivedParams.MatrixA_RowIndices = int32(0);
                                derivedParams.MatrixA_RowStartLoc = int32(0);
                                derivedParams.MatrixA_ColumnSum = int32(0);
                                derivedParams.MatrixB_RowIndices = int32(0);
                                derivedParams.MatrixB_RowStartLoc = int32(0);
                                derivedParams.MatrixB_ColumnSum = int32(0);
                                coder.internal.errorIf(~invertible, 'comm:getLDPCEncoderParameters:NonInvertibleParityCheckMatrix');
                            else
                                % PB is invertible in GF(2), and has been modified by gf2factorize.
                                [derivedParams.MatrixL_RowIndices, ...
                                    derivedParams.MatrixL_RowStartLoc, ...
                                    derivedParams.MatrixL_ColumnSum] = ConvertMatrixFormat(tril(PL, -1));
                                PD = PC(rowOrder, :); % Need to do this before the next line.
                                P = triu(PD, 1);      % Now PB is upper triangular. Remove the diagonal and put it in P.
                                derivedParams.RowOrder = int32(rowOrder-1);
                                [derivedParams.MatrixA_RowIndices, ...
                                    derivedParams.MatrixA_RowStartLoc, ...
                                    derivedParams.MatrixA_ColumnSum] = ConvertMatrixFormat(obj.ParityCheckMatrix(:,1:K));
                                [derivedParams.MatrixB_RowIndices, ...
                                    derivedParams.MatrixB_RowStartLoc, ...
                                    derivedParams.MatrixB_ColumnSum] = ConvertMatrixFormat(P);
                            end
                    end
            end

            %------------------------------------------------------------------
            % Local functions
            function shape = isfulldiagtriangular(X)
                % X must be a square logical matrix.
                % shape = 1  if X is lower triangular and has a full diagonal.
                % shape = -1 if X is upper triangular and has a full diagonal.
                % Otherwise, shape = 0.
                % If X is a diagonal matrix, X is considered upper triangular.

                if ~all(diag(X)) % Full diagonal?
                    shape = 0;   % Must have a full diagonal.
                else
                    NumNonzerosInLowerPart = nnz(tril(X));
                    if NumNonzerosInLowerPart == nnz(X)
                        shape = 1;  % X is lower triangular.
                    elseif NumNonzerosInLowerPart == size(X,1)
                        shape = -1; % X is upper triangular.
                    else
                        shape = 0;  % X is not triangular.
                    end
                end
            end

            function [A, B, chosen_pivots, invertible] = gf2factorize(X)
                %GF2FACTORIZE  Factorize a square matrix in GF(2).
                %   [A B chosen_pivots invertible] = gf2factorize(X) factorizes a square matrix X
                %   in GF(2) using Gaussian elimination.
                %
                %   X = A * B    (in GF(2), i.e., using modulo-2 arithmetic)
                %
                %   X may be a sparse matrix. Nonzero elements are treated as ones.
                %   A and B are sparse logical matrices.
                %
                %   A is always lower triangular.
                %   If X is invertible in GF(2), then B(chosen_pivots,:) is upper triangular and
                %   invertible = true.
                %   If X is non-invertible in GF(2), then X = A * B still holds and invertible = false.
                %
                %   To evaluate A * B in GF(2), use mod(A*double(B),2).

                n = size(X,1);
                Y1 = eye(n,'logical');
                Y2 = full(X~=0);
                chosen_pivots = zeros(n,1); % Haven't chosen any pivots yet.
                invertible = true; % Assume that X is invertible in GF(2).

                for col = 1:n
                    candidate_rows = Y2(:,col);
                    candidate_rows(chosen_pivots(1:(col-1))) = 0; % Never use a chosen pivot.
                    candidate_rows = find(candidate_rows); % Convert into row indices.

                    if isempty(candidate_rows)
                        invertible = false; % X is not invertible in GF(2).
                        % Output sparse matrices to save memory.
                        A = sparse(Y1);
                        B = sparse(Y2);
                        % Output the chosen pivots even if X is not invertible in GF(2).
                        d = setdiff((1:n)', sort(chosen_pivots(1:(col-1))));
                        chosen_pivots(col+(0:length(d)-1)) = d;
                        return;
                    else
                        pivot = candidate_rows(1);      % Choose the first candidate row as the pivot row.
                        chosen_pivots(col) = pivot;     % Record this pivot.

                        % Find all nonzero elements in the pivot row.
                        % They will be xor-ed with the corresponding elements in other
                        % candidate rows.
                        columnind = find(Y2(pivot,:));

                        % Subtract the pivot row from all other candidate_rows.
                        % Exploit the fact that we are working in GF(2).
                        % Just use logical NOT.
                        % As we continue in this loop, Y2 will become "psychologically"
                        % upper triangular.

                        if ~isscalar(candidate_rows)
                            Y2(candidate_rows(2:end), columnind) = not(Y2(candidate_rows(2:end), columnind));
                        end

                        % Update the lower triangular matrix Y1.
                        Y1(candidate_rows(2:end), pivot) = 1;
                    end
                end

                % Output sparse matrices to save memory.
                A = sparse(Y1);
                B = sparse(Y2);
            end

            function [RowIndices, RowStartLoc, ColumnSum] = ConvertMatrixFormat(X)
                % Create an alternative representation of a zero-one matrix.

                % Find nonzero elements.
                [i,~] = find(X);

                % Generate zero-based row indices and use int32 for C-MEX function.
                if ~isempty(i)
                    RowIndices = int32(i-1);
                else
                    % In this case, RowIndices' value doesn't matter. But we still return
                    % a nonempty matrix to avoid error from Simulink block mask.
                    RowIndices = int32(0);
                end

                % Find the number of nonzero elements in each column.
                % Use int32 for C-MEX function.
                ColumnSum  = int32(full(sum(X)));

                % For each column, find where the corresponding row indices start in RowIndices.
                % Generate zero-based indices and use int32 for C-MEX function.
                CumulativeSum = cumsum(double(ColumnSum));
                RowStartLoc = int32([0 CumulativeSum(1:end-1)]);
            end
        end
    end
end