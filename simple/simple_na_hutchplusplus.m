function trace_est=simple_na_hutchplusplus(A, num_queries)
% Estimates the trace of square matrix A with num_queries many matrix-vector products.
% Implements the Non-Adaptive Hutch++ Algorithm, where c1 = 1/6, c2 = 1/3, c3 = 1/2
% Random sign vectors are used.

	% Generate random sign matrices
	% c1, c2, and c3 appear in the floor/ceil expression below:
    S = 2*randi(2,size(A,1),floor(num_queries/6))-3;
    R = 2*randi(2,size(A,1),floor(num_queries/3))-3;
    G = 2*randi(2,size(A,1),ceil(num_queries/2))-3;
    
    % Compute NA-Hutch++ Estimator
    Z = A*R;
    W = A*S;
    trace_est = trace(pinv(S'*Z)*(W'*Z)) + 1/size(G,2)*[trace(G'*A*G) - trace((G'*Z)*pinv(S'*Z)*(W'*G))];

end  % simple_na_hutchplusplus
