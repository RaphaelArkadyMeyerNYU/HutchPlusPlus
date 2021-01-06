function compare_estimators_number_of_triangles(matrix, num_queries, num_trials)
% compare_estimators_number_of_triangles(matrix, num_queries, num_trials)
% 
% Compares the estimates achieved by the 4 trace estimators in the `core`
% folder on repeated trials, when estimating the trace of the cube of a matrix.
% 
% In order to output relative errors, the true trace of the cube of the given
% matrix is computed exactly, which may be prohibitively expensive for large matrices.
% 
% Required Inputs:
% - matrix: the matrix whose trace-of-cube we are estimating. This should be 
% a square matrix, not a function_handle.
% 
% - num_queries: the number of matrix-vector products computed with `matrix^3`. This
% is NOT the total number of matrix-vector products computed with `matrix`, because
% 3 matrix-vector products with `matrix` are needed to compute a single matrix-vector
% product with `matrix^3`.
%
% Optional Inputs:
% 
% - num_trials: The number of time to run each estimator, allowing us to analyze
% the variance in the errors. (default value: 50)
% 

	arguments
		matrix;
		num_queries;
		num_trials = 50;
	end

	% Possibly very expensive
	true_trace = trace(matrix * matrix * matrix);

	matrix_cubed = @(x) matrix*(matrix*(matrix*x));
	compare_estimators_on_matvec_oracle(matrix_cubed, num_queries, size(matrix,1), num_trials, 'true_trace', true_trace, 'objective_name', 'number of triangles')
end
