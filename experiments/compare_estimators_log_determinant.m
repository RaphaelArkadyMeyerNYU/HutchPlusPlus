function compare_estimators_log_determinant(matrix, num_queries, num_trials, lanczos_iterations)
% compare_estimators_log_determinant(matrix, num_queries, num_trials, lanczos_iterations)
% 
% Compares the estimates achieved by the 4 trace estimators in the `core`
% folder on repeated trials, when estimating the log-determinants of a matrix.
% 
% The log-determinant is estimated as the trace of the logarithm of a matrix,
% and uses Lanczos iteration to approximate matrix-vector product with this
% logarithm matrix.
% 
% Required Inputs:
% - matrix: the matrix whose log-determinant we are estimating. This should be 
% a square matrix, not a function_handle.
% 
% - num_queries: the number of matrix-vector products computed with `log(matrix)`. This
% is NOT the total number of matrix-vector products computed with `matrix`, because
% `lanczos_iterations` many matrix-vector products with `matrix` are needed to compute a
% single matrix-vector product with `log(matrix)`.
%
% Optional Inputs:
% 
% - num_trials: The number of time to run each estimator, allowing us to analyze
% the variance in the errors. (default value: 50)
% 
% - lanczos_iterations: The number of lanczos iterations used to approximate a
% matrix-vector product with `log(matrix)`. (default value: 50)
% 

	arguments
		matrix;
		num_queries;
		num_trials = 50;
		lanczos_iterations = 20;
	end

	logHandle = @(B) lanczos(matrix, B, @log, lanczos_iterations);
	compare_estimators_on_matvec_oracle(logHandle, num_queries, size(matrix,1), num_trials, 'objective_name', 'log-determinant')
end
