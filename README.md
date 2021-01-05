This is the public code created for [Hutch++, a paper on trace estimation, found online here](https://arxiv.org/abs/2010.09649).


This repo contains MATLAB code to recreate the paper's experiments, as well as simplified code that makes it easier to use Hutch++ for your next application (or translate it to your language of choice).

## Simple Code

The `simple` directory contains minimized scripts for the algorithms presented, in appropriately named files:

- `simple\simple_hutchinson.m`: Hutchinson's Estimator
- `simple\simple_hutchplusplus.m`: The Hutch++ Algorithm (Algorithm 1 from the paper)
- `simple\simple_na_hutchplusplus.m`: The NA-Hutch++ Algorithm (Algorithm 2 from the paper)
- `simple\simple_subspace_projection.m`: The Subspace Projection Algorithm (with exactly one iteration; from ["Randomized matrix-free trace and log-determinant estimators"](https://arxiv.org/abs/1605.04893))


## Core Code

The `core` directory contains more flexible and better commented scripts for the algorithms presented, in appropriately named files.
These scripts accept both matrices and function handles as inputs.
They also allow the user to specify certain optional parameters, like the exact distributions used in sketching or Hutchinson steps, and the number of subspace projection iterations to run.

- `core\hutchinson.m`: Hutchinson's Estimator
	- <i>Optional</i>: Specify the distribution of random vectors
- `core\hutchplusplus.m`: The Hutch++ Algorithm (Algorithm 1 from the paper)
	- <i>Optional</i>: Specify the fraction of queries used for sketching, the distribution of random Hutchinson vectors, and the distribution of random sketching vectors.
- `core\na_hutchplusplus.m`: The NA-Hutch++ Algorithm (Algorithm 2 from the paper)
	- <i>Optional</i>: Specify the constants c<sub>1</sub> and c<sub>2</sub>, the distribution of random Hutchinson vectors, and the distribution of random sketching vectors.
- `core\subspace_projection.m`: The Subspace Projection Algorithm (with exactly one iteration; from ["Randomized matrix-free trace and log-determinant estimators"](https://arxiv.org/abs/1605.04893))
	- <i>Optional</i>: Specify the number of subspace iteration steps, and the distribution of random sketching vectors.


## Experiment Code

The `experiments` directory has scripts that measure the error of these four stochastic trace estimators.

- `experiments\compare_estimators_on_mat.m`: Outputs the median, IQR, 25th quartile, and 75th quartile of the relative errors achieved by the four estimators on repeated trials. Only works for explicit matrices.
- `experiments\compare_estimators_on_matvec_oracle.m`: Outputs the median, IQR, 25th, and 75th quartiles of the relative errors achieved by the four estimators on repeated trials. Only works for function handles.
	- If the true trace of the underlying matrix is unknown, this script will output the median, IQR, 25th quartile, and 75th quartile of the actual outputs of the estimators (instead of the errors).

We measure the relative error of estimated trace as the ratio |trace<sub>true</sub> - trace<sub>estimate</sub>| / trace<sub>true</sub>

## LaTeX Styling

Aesthetically, LaTeX does not type "Hutch++" in a very nice way.
It makes the plus symbols too large.
For aesthetic reasons, we use the following script to define the Hutch++ LaTeX command, which shrinks and raises the + symbols:
```
\usepackage{relsize} % Relsize package needed

\newcommand{\hutchpp}{\text{Hutch\raisebox{0.35ex}{\relscale{0.75}++}}}
```
