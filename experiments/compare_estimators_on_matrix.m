function compare_estimators_on_matrix(matrix, num_queries, num_trials, args)
% compare_estimators_on_matrix(matrix, num_queries, num_trials, args)
% 
% Compares the relative errors achieved by the 4 trace estimators in the `core`
% folder on repeated trials, for a given square matrix.
% 
% This function requires an explicit matrix, not a function handle. See
% `compare_estimators_on_matvec_oracle` if you need to pass a function handle.
% 
% Required Inputs:
% - matrix: a square matrix
% 
% - num_queries: Total number of matrix-vector products to compute
% 
% Optional Inputs:
% 
% - num_trials: The number of time to run each estimator, allowing us to analyze
% the variance in the errors. (default value: 50)
% 
% Name-Value Pair Inputs:
% 
% - hutch_dist, sketch_dist, sketch_frac, c1, c2, sketch_iterations: optional parameters
% that are passed down to hutchinson(), hutchplusplus(), na_hutchplusplus(), and
% subspace_projection().
% 
% Examples:
% 
% Let A be a matrix (not a function handle).
% 
% Compare the estimators on A with 40 matrix-vector products
%     compare_estimators_on_matrix(A, 40)
% 
% Compare the estimators on A with 36 matrix-vector products, using 100 trials
%     compare_estimators_on_matrix(A, 36, 100)
%
% Compare the estimators on A with 52 matrix-vector products, using 100 trials,
% only using Gaussian vectors. Hutch++ sketches with 15% of its vectors, NA-Hutch++
% uses c1=1/4 and c2=1/2, and subspace projection uses 2 iterations.
% 
%     compare_estimators_on_matrix(A, 52, 100, 'hutch_dist', @randn, 'sketch_dist', @randn, 'sketch_frac', 0.15, 'c1', 1/4, 'c2', 1/2, 'sketch_iterations', 2)
% 

    arguments
        matrix;
        num_queries;
        num_trials = 50;
        args.hutch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
        args.sketch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
        args.sketch_frac = 2/3;
        args.c1 = 1/6
        args.c2 = 1/3
        args.sketch_iterations = 1;
	end

	true_trace = trace(matrix);
	trials = zeros(num_trials, 4); % 4 = number of trace estimators we consider

	for t=1:num_trials
		trials(t,1) = hutchinson(matrix, num_queries, 'hutch_dist', args.hutch_dist);
		trials(t,2) = hutchplusplus(matrix, num_queries, 'hutch_dist', args.hutch_dist, 'sketch_dist', args.sketch_dist, 'sketch_frac', args.sketch_frac);
		trials(t,3) = na_hutchplusplus(matrix, num_queries, 'hutch_dist', args.hutch_dist, 'sketch_dist', args.sketch_dist, 'c1', args.c1, 'c2', args.c2);
		trials(t,4) = subspace_projection(matrix, num_queries, 'sketch_dist', args.sketch_dist, 'sketch_iterations', args.sketch_iterations);
	end

	relative_errors = abs(trials - true_trace) / true_trace;
	bottom_quartiles = quantile(relative_errors, 0.25);
	medians = quantile(relative_errors, 0.5);
	top_quartiles = quantile(relative_errors, 0.75);
	iqrs = iqr(relative_errors);

	real2str = @(x) sprintf('%0.3f',x);

	fprintf("After " + num_trials + " trials using " + num_queries + " queries, we find that:\n")
	fprintf("Hutchinson      achieves relative error in [" + real2str(bottom_quartiles(1)) + "," + real2str(top_quartiles(1)) + "] 50%% of the time (iqr=" + real2str(iqrs(1)) + ", median=" + real2str(medians(1)) +")\n");
	fprintf("Hutch++         achieves relative error in [" + real2str(bottom_quartiles(2)) + "," + real2str(top_quartiles(2)) + "] 50%% of the time (iqr=" + real2str(iqrs(2)) + ", median=" + real2str(medians(2)) +")\n");
	fprintf("NA-Hutch++      achieves relative error in [" + real2str(bottom_quartiles(3)) + "," + real2str(top_quartiles(3)) + "] 50%% of the time (iqr=" + real2str(iqrs(3)) + ", median=" + real2str(medians(3)) +")\n");
	fprintf("SubspaceProject achieves relative error in [" + real2str(bottom_quartiles(4)) + "," + real2str(top_quartiles(4)) + "] 50%% of the time (iqr=" + real2str(iqrs(4)) + ", median=" + real2str(medians(4)) +")\n");

end  % compare_estimators_on_matrix
