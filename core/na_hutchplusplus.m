function trace_est=na_hutchplusplus(matVecOracle, num_queries, dimension, args)
% trace_est = na_hutchplusplus(matVecOracle, num_queries, dimension, args)
%
% Runs the NA-Hutch++ Algorithm to estimate the trace of the given square
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
% matrix that will be used in the sketching matrices. (default value: random sign matrix)
% 
% - c1: The constant c1 in Algorithm 2 of the paper: the fraction of queries
% devoted to matrix S. (default value: 1/6)
% 
% - c2: The constant c2 in Algorithm 2 of the paper: the fraction of queries
% devoted to matrix R. (default value: 1/3)
%
% Note that c3 = 1 - c1 - c2 is calculated by the script. Further, note that
% c2 / c1 should be a large enough constant (e.g. 2, like in the default values)
% 
% Examples:
% 
% Let A be a matrix (not a function handle).
% 
% Estimate the trace of A with 100 matrix-vector queries, using default parameters:
%     na_hutchplusplus(A, 100)
% 
% Estimate the trace of the cube of A with 100 matrix-vector queries, with default parameters.
% We create a function handle for the cube of A as an oracle, and explicitly pass the dimension of A.
%     A_cubed = @(x) A*(A*(A*x))
%     na_hutchplusplus(A_cubed, 100, size(A,1))
% 
% Estimate the trace of A with 52 matrix-vector queries, using only Gaussian vectors, and with constants
% c1 = 1/4, c2 = 1/2, c3 = 1/2
% 
%     na_hutchplusplus(A, 52, 'hutch_dist', @randn, 'sketch_dist', @randn, 'c1', 1/4, 'c2', 1/2)
% 

    % Set default values
    arguments
        matVecOracle;
        num_queries;
        dimension = -1;
        args.hutch_dist = @(m,n) 2*randi(2,m,n)-3;  % Random sign matrix
        args.sketch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
        args.c1 = 1/6
        args.c2 = 1/3
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
    R_num_queries = round(num_queries * args.c1);
    S_num_queries = round(num_queries * args.c2);
    Hutch_num_queries = num_queries - R_num_queries - S_num_queries;

    % Sample from the corresponding distributions
    R = args.sketch_dist(dimension, R_num_queries);
    S = args.sketch_dist(dimension, S_num_queries);
    G = args.hutch_dist(dimension, Hutch_num_queries);

    % Compute matrix-vector products for sketching
    Z = matVecOracle(R);
    W = matVecOracle(S);

    % Compute the full trace estimate
    trace_est = trace(pinv(S'*Z)*(W'*Z)) + 1/Hutch_num_queries*[trace(G'*matVecOracle(G)) - trace((G'*Z)*pinv(S'*Z)*(W'*G))];

end  % na_hutchplusplus
