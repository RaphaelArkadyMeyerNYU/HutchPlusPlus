function trace_est=simple_subspace_projection(A, num_queries)
% Estimates the trace of square matrix A with num_queries many matrix-vector products.
% Implements the Subspace Projection Algorithm from 
% 	"Randomized matrixfree trace and log-determinant estimators".
% Random sign vectors are used.
    
    % Generate random sign matrix
    S = 2*randi(2,size(A,1),ceil(num_queries/2))-3;

    % Subspace Project
    [Q,R] = qr(A*S,0);
    trace_est = trace(Q'*A*Q);
    
end %simple_subspace_projection
