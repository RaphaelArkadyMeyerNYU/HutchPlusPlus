function trace_est=simple_hutchplusplus(A, num_queries)
% Estimates the trace of square matrix A with num_queries many matrix-vector products
% Implements the Hutch++ algorithm. Random sign vectors are used.
    
    % Calculate which matrices get how many queries, and generate random sign matrices
    S = 2*randi(2,size(A,1),ceil(num_queries/3))-3;
    G = 2*randi(2,size(A,1),floor(num_queries/3))-3;
    
    % Compute only the Q of the QR decomposition
    [Q,_] = qr(A*S,0);
    G = G - Q*(Q'*G);

    trace_est = trace(Q'*A*Q) + 1/size(G,2)*trace(G'*A*G);
    
end  % simple_hutchplusplus
