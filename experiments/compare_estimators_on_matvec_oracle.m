function compare_estimators_on_matvec_oracle(matVecOracle, num_queries, dimension, num_trials, args)
% compare_estimators_on_matvec_oracle(matVecOracle, num_queries, dimension, num_trials, args)
% 
% Compares the relative errors achieved by the 4 trace estimators in the `core`
% folder on repeated trials, for a given function handle to a matrix.
% 
% This function requires a function handle, not an explicit matrix. See
% `compare_estimators_on_matrix` if you need to pass an explicit matrix.
% 
% Required Inputs:
% - matVecOracle: A function_handle that, on input B returns the matrix A*B
% 
% - num_queries: Total number of matrix-vector products to compute
% 
% - dimension: dimension of the matrix inside matVecOracle
% 
% Optional Inputs:
% 
% - num_trials: The number of time to run each estimator, allowing us to analyze
% the variance in the errors. (default value: 50)
% 
% Name-Value Pair Inputs:
% 
% - true_trace: The true trace of the matrix inside matVecOracle. If specified, the
% output will be measured in terms of relative error: |estimate_trace - true_trace|/true_trace
% If not specified, the output will be described in terms of the output value of the estimators,
% making it harder to see if one estimator is more biased than another estimator.
% 
% - hutch_dist, sketch_dist, sketch_frac, c1, c2, sketch_iterations: optional parameters
% that are passed down to hutchinson(), hutchplusplus(), na_hutchplusplus(), and
% subspace_projection().
% 
% Examples:
% 
% Let A be a not a function handle (not a matrix) for a matrix with dimension 101 and trace 3.14
% 
% Compare the estimators on A with 40 matrix-vector products
%     compare_estimators_on_matvec_oracle(A, 40, 101)
%
% Compare the estimators on A with 40 matrix-vector products, in terms of relative error
%     compare_estimators_on_matvec_oracle(A, 40, 101, 'true_trace', 3.14)
% 
% Compare the estimators on A with 36 matrix-vector products, using 75 trials
%     compare_estimators_on_matvec_oracle(A_cubed, 40, 101, 75);
% 
% Compare the estimators on A with 36 matrix-vector products, using 75 trials, in terms of relative error
%     compare_estimators_on_matvec_oracle(A_cubed, 40, 101, 75, 'true_trace', 3.14);
%
% Compare the estimators on A with 52 matrix-vector products, using 100 trials,
% only using Gaussian vectors, with output in terms of relative error. Hutch++
% sketches with 15% of its vectors, NA-Hutch++ uses c1=1/4 and c2=1/2, and
% subspace projection uses 2 iterations.
% 
%     compare_estimators_on_matvec_oracle(A, 52, 101, 100, 'true_trace', 3.14, 'hutch_dist', @randn, 'sketch_dist', @randn, 'sketch_frac', 0.15, 'c1', 1/4, 'c2', 1/2, 'sketch_iterations', 2)
%

    arguments
        matVecOracle;
        num_queries;
        dimension;
        num_trials = 50;
        args.true_trace = NaN;
        args.hutch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
        args.sketch_dist = @(m,n) 2*randi(2,m,n)-3; % Random sign matrix
        args.sketch_frac = 2/3;
        args.c1 = 1/6
        args.c2 = 1/3
        args.sketch_iterations = 1;
	end

	trials = zeros(num_trials, 4); % 4 = number of trace estimators we consider

	for t=1:num_trials
		trials(t,1) = hutchinson(matVecOracle, num_queries, dimension, 'hutch_dist', args.hutch_dist);
		trials(t,2) = hutchplusplus(matVecOracle, num_queries, dimension, 'hutch_dist', args.hutch_dist, 'sketch_dist', args.sketch_dist, 'sketch_frac', args.sketch_frac);
		trials(t,3) = na_hutchplusplus(matVecOracle, num_queries, dimension, 'hutch_dist', args.hutch_dist, 'sketch_dist', args.sketch_dist, 'c1', args.c1, 'c2', args.c2);
		trials(t,4) = subspace_projection(matVecOracle, num_queries, dimension, 'sketch_dist', args.sketch_dist, 'sketch_iterations', args.sketch_iterations);
	end


	if isnan(args.true_trace)
		% If the true trace is unkown, then report the statistics of the raw outputs

		real2str = @(x) sprintf('%0.2e',x);

		bottom_quartiles = quantile(trials, 0.25);
		medians = quantile(trials, 0.5);
		top_quartiles = quantile(trials, 0.75);
		iqrs = iqr(trials);

		fprintf("After " + num_trials + " trials using " + num_queries + " queries, we find that:\n")
		fprintf("Hutchinson      estimates trace to be in [" + real2str(bottom_quartiles(1)) + "," + real2str(top_quartiles(1)) + "] 50%% of the time (iqr=" + real2str(iqrs(1)) + ", median=" + real2str(medians(1)) +")\n");
		fprintf("Hutch++         estimates trace to be in [" + real2str(bottom_quartiles(2)) + "," + real2str(top_quartiles(2)) + "] 50%% of the time (iqr=" + real2str(iqrs(2)) + ", median=" + real2str(medians(2)) +")\n");
		fprintf("NA-Hutch++      estimates trace to be in [" + real2str(bottom_quartiles(3)) + "," + real2str(top_quartiles(3)) + "] 50%% of the time (iqr=" + real2str(iqrs(3)) + ", median=" + real2str(medians(3)) +")\n");
		fprintf("SubspaceProject estimates trace to be in [" + real2str(bottom_quartiles(4)) + "," + real2str(top_quartiles(4)) + "] 50%% of the time (iqr=" + real2str(iqrs(4)) + ", median=" + real2str(medians(4)) +")\n");

	else
		% If the true trace is given, then report relative errors
		relative_errors = abs(trials - args.true_trace) / args.true_trace;
		real2str = @(x) sprintf('%0.3f',x);

		bottom_quartiles = quantile(relative_errors, 0.25);
		medians = quantile(relative_errors, 0.5);
		top_quartiles = quantile(relative_errors, 0.75);
		iqrs = iqr(relative_errors);

		fprintf("After " + num_trials + " trials using " + num_queries + " queries, we find that:\n")
		fprintf("Hutchinson      achieves relative error in [" + real2str(bottom_quartiles(1)) + "," + real2str(top_quartiles(1)) + "] 50%% of the time (iqr=" + real2str(iqrs(1)) + ", median=" + real2str(medians(1)) +")\n");
		fprintf("Hutch++         achieves relative error in [" + real2str(bottom_quartiles(2)) + "," + real2str(top_quartiles(2)) + "] 50%% of the time (iqr=" + real2str(iqrs(2)) + ", median=" + real2str(medians(2)) +")\n");
		fprintf("NA-Hutch++      achieves relative error in [" + real2str(bottom_quartiles(3)) + "," + real2str(top_quartiles(3)) + "] 50%% of the time (iqr=" + real2str(iqrs(3)) + ", median=" + real2str(medians(3)) +")\n");
		fprintf("SubspaceProject achieves relative error in [" + real2str(bottom_quartiles(4)) + "," + real2str(top_quartiles(4)) + "] 50%% of the time (iqr=" + real2str(iqrs(4)) + ", median=" + real2str(medians(4)) +")\n");
	end
end  % compare_estimators_on_matvec_oracle
