function trace_est=hutchinson(matVecOracle, num_queries, dimension, args)
% trace_est = hutchinson(matVecOracle, num_queries, dimension, args)
% 
% Compute the Hutchinson Estimator for the trace of a given square
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
% Examples:
% 
% Let A be a matrix (not a function handle).
% 
% Estimate the trace of A with 100 matrix-vector queries, using default parameters:
%     hutchinson(A, 100)
% 
% Estimate the trace of the cube of A with 100 matrix-vector queries, with default parameters.
% We create a function handle for the cube of A as an oracle, and explicitly pass the dimension of A.
%     A_cubed = @(x) A*(A*(A*x))
%     hutchinson(A_cubed, 100, size(A,1))
% 
% Estimate the trace of A with 52 matrix-vector queries, using only Gaussian vectors,
% 
%     hutchinson(A, 52, 'hutch_dist', @randn)
% 

    % Set default values
    arguments
        matVecOracle;
        num_queries;
        dimension = -1;
        args.hutch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
    end

    % Force the input matrix to be a function handle
    if ~isa(matVecOracle,'function_handle')
        dimension = size(matVecOracle,1);
        matVecOracle = @(z) matVecOracle*z;
    elseif (dimension == -1)
        error('hutchinson_estimator:HandleWithoutDimension', ...
              'Passed MatVec handle but not the matrix dimension. ' + ...
              '"dimension" is the third (optional) input argument.');
    end
    % From now on, matVecOracle is to be treated as a function_handle to the matrix

    % Sample random queries
    G = args.hutch_dist(dimension, num_queries);
    trace_est = trace(G' * matVecOracle(G)) / num_queries;
    
end  % hutchinson
