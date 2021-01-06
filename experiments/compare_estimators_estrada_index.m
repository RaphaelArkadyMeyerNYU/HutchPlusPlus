function compare_estimators_estrada_index(matrix, num_queries, num_trials, lanczos_iterations)
% compare_estimators_estrada_index(matrix, num_queries, num_trials, lanczos_iterations)
% 
% Compares the relative errors achieved by the 4 trace estimators in the `core`
% folder on repeated trials, when estimating the Estrada Index of a matrix.
% 
% The Estrada index is estimated as the trace of the exponential of a matrix,
% and uses Lanczos iteration to approximate matrix-vector product with this
% exponential matrix.
% 
% Required Inputs:
% - matrix: the matrix whose Estrada index we are estimating. This should be 
% a square matrix, not a function_handle.
% 
% - num_queries: the number of matrix-vector products computed with `e^matrix`. This
% is NOT the total number of matrix-vector products computed with `matrix`, because
% `lanczos_iterations` many matrix-vector products with `matrix` are needed to compute a
% single matrix-vector product with `e^matrix`.
%
% Optional Inputs:
% 
% - num_trials: The number of time to run each estimator, allowing us to analyze
% the variance in the errors. (default value: 50)
% 
% - lanczos_iterations: The number of lanczos iterations used to approximate a
% matrix-vector product with `e^matrix`. (default value: 50)
% 

	arguments
		matrix;
		num_queries;
		num_trials = 50;
		lanczos_iterations = 50;
	end

	expHandle = @(B) lanczos(matrix, B, @exp, lanczos_iterations);
	compare_estimators_on_matvec_oracle(expHandle, num_queries, size(matrix,1), num_trials, 'objective_name', 'Estrada Index')
end
