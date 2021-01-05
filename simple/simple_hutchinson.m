function trace_est=simple_hutchinson(A, num_queries)
% Estimates the trace of square matrix A with num_queries many matrix-vector products.
% Implements Hutchinson's estimator. Random sign vectors are used.

	% Generate a random sign matrix
    S = 2*randi(2,size(A,1),num_queries)-3;

    % Apply Hutchinson's estimator
    trace_est = trace(S'*(A*S)) / num_queries;
    
end  % simple_hutchinson
