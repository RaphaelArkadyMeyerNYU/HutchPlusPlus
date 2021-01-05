function trace_est=subspace_projection(matVecOracle, num_queries, dimension, args)
% trace_est = subspace_projection(matVecOracle, num_queries, dimension, args)
% 
% Runs the subspace projection algorithm from [SAI17] to estimate the trace of a given
% square matrix using matrix-vector queries.
% 
% [SAI17]: Arvind K. Saibaba, Alen Alexanderian, and Ilse C. F. Ipsen. Randomized matrix-
% free trace and log-determinant estimators. Numerische Mathematik, 137(2):353–395, 2017.
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
% - sketch_dist: Function_handle that, given two inputs (m,n), returns a m by n
% matrix that will be used to generate random vectors to start with. (default value: random sign matrix)
% 
% - sketch_iterations: Number of iterations of subspace iteration to perform. In the paper,
% we only used one iteration. For fixed num_queries, increasing sketch_iterations decreases
% the size of the subspace we estimate, since it takes k(q+1) queries to estimate the top
% k dimensions of A with q iterations. (default value: 1)
% 
% Examples:
% 
% Let A be a matrix (not a function handle).
% 
% Estimate the trace of A with 100 matrix-vector queries, using default parameters:
%     subspace_projection(A, 100)
% 
% Estimate the trace of the cube of A with 100 matrix-vector queries, with default parameters.
% We create a function handle for the cube of A as an oracle, and explicitly pass the dimension of A.
%     A_cubed = @(x) A*(A*(A*x))
%     subspace_projection(A_cubed, 100, size(A,1))
% 
% Estimate the trace of A with 52 matrix-vector queries and 3 steps of subspace iteration, using only
% Gaussian vectors.
% 
%     subspace_projection(A, 52, 'sketch_dist', @randn, 'sketch_iterations', 3)
% 

    % Set default values
    arguments
        matVecOracle;
        num_queries;
        dimension = -1;
        args.sketch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
        args.sketch_iterations = 1;
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

    % Dimension of the subspace we estimate
    subspace_dimension = floor(num_queries / (1 + args.sketch_iterations));
    
    % Generate initial matrix
    Q = args.sketch_dist(dimension, subspace_dimension);

    % Iterate
    for t=1:args.sketch_iterations
        [Q, ~] = qr(matVecOracle(Q), 0); % 0 here means 'use economic qr'
    end

    % Estimate
    trace_est = trace(Q' * matVecOracle(Q));
    
end  % subspace_projection
