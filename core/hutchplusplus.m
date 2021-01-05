function trace_est=hutchplusplus(matVecOracle, num_queries, dimension, args)
% trace_est = hutchplusplus(matVecOracle, num_queries, dimension, args)
% 
% Runs the Hutch++ Algorithm to estimate the trace of the given square
% matrix using matrix-vector queries.
% 
% Required Inputs:
% - matVecOracle: Either a matrix A, or a function_handle that, on input B
% returns the matrix A*B. This function_handle is our "Matrix-Vector Oracle"
% 
% - num_queries: Total number of matrix-vector products to compute
% 
% Optional Inputs:
% 
% - dimension: Dimension of the input matrix. Required if matVecOracle is
% a function_handle, instead of a function.
% 
% Name-Value Pair Inputs:
% 
% - hutch_dist: Function_handle that, given two inputs (m,n), returns a m by n
% matrix that will be used in the Hutchinson estimation. (default value: random sign matrix)
% 
% - sketch_dist: Function_handle that, given two inputs (m,n), returns a m by n
% matrix that will be used in the sketching matrix. (default value: random sign matrix)
% 
% - sketch_fraction: Fraction of queries to commit towards estimating top eigenvalues, as
% opposed to Hutchinson estimation. This is twice the number of columns S in Algorithm 1 of
% the paper. (default value: 2/3)
%
% Examples:
% 
% Let A be a matrix (not a function handle).
% 
% Estimate the trace of A with 100 matrix-vector queries, using default parameters:
%     hutchplusplus(A, 100)
% 
% Estimate the trace of the cube of A with 100 matrix-vector queries, with default parameters.
% We create a function handle for the cube of A as an oracle, and explicitly pass the dimension of A.
%     A_cubed = @(x) A*(A*(A*x))
%     hutchplusplus(A_cubed, 100, size(A,1))
% 
% Estimate the trace of A with 52 matrix-vector queries, using only Gaussian vectors,
% with 1/10 of the vectors used for sketching.
%
%     hutchplusplus(A, 52, 'hutch_dist', @randn, 'sketch_dist', @randn, 'sketch_frac', 1/10)
%

    % Set default values
    arguments
        matVecOracle;
        num_queries;
        dimension = -1;
        args.hutch_dist = @(m,n) 2*randi(2,m,n)-3;  % Random sign matrix
        args.sketch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
        args.sketch_frac = 2/3;
    end

    % Force the input matrix to be a function handle
    if ~isa(matVecOracle,'function_handle')
        dimension = size(matVecOracle,1);
        matVecOracle = @(z) matVecOracle*z;
    elseif (dimension == -1)
        error('hutchplusplus_estimator:HandleWithoutDimension', ...
              'Passed MatVec handle but not the matrix dimension. ' + ...
              '"dimension" is the third (optional) input argument.');
    end
    % From now on, matVecOracle is to be treated as a function_handle to the matrix

    % Calculate which matrices get how many queries
    S_num_queries = round(num_queries * args.sketch_frac / 2);
    Hutch_num_queries = num_queries - S_num_queries;

    % Sketch the matrix, and take the QR
    S = args.sketch_dist(dimension, S_num_queries);
    [Q, ~] = qr(matVecOracle(S), 0); % 0 here means 'use economic qr'

    % Sample a matrix, then project it away from top eigenvalues
    G = args.hutch_dist(dimension, Hutch_num_queries);
    G = G - Q*(Q' * G);

    % Compute Hutch++ Estimate value
    trace_est = trace(Q' * matVecOracle(Q)) + trace(G' * matVecOracle(G)) / Hutch_num_queries;

end  % hutchplusplus
