function f = lanczos(A, B, func, iter)
% f = lanczos(A, B, func, iter)
% 
% Approximates func(A) * B by Lanczos iteration, where func(A)
% denotes a matrix with the same eigenvectors as A, but whose
% i^th eigenvalue is 'func' applied to the i^th eigenvalue of A.
% 
% Notably, func is a function on the real line. It maps
% eigenvalues to eigenvalues, and is not an arbitrary map
% from matrices to matrices.
% 
% This is a modification of the Lanczos code available at
% 		github.com/cpmusco/fast-pcr
% where this code allows for block matrix-vector products
% (i.e. B is allowed to me a matrix, not just a column vector)
% 
% Required Inputs:
% - A: Square input matrix.
% 
% - B: Rectangular input matrix. A*B should be a valid computation.
% 
% - func: function_handle that given a real eigenvalue, returns a new eigenvalue.
%
% Examples:
% 
% Let A be a matrix with n rows and n columns.
% Let x be a column vector with n entries.
% Let B be a matrix with n rows and k columns.
% 
% Estimate e^A * x with 10 iterations of Lanczos:
% 	lanczos(A, x, @exp, 10)
% 
% Estimate log(A)*B with 23 iterations of Lanczos. In order to avoid numerical issues
% with eigenvalues close to zero, slightly offset the input eigenvalues to the log function:
%   lanczos(A, B, @(t) log(t + 0.008), 23);
%

	% Count dimensions. Allocate Space.
	n = size(B,1);
	d = size(B,2); 
	beta = zeros(iter,d);
	alpha = zeros(iter,d);

	% Build Krylov subspace
	K = zeros(n,d,iter+1);
	K(:,:,1) = normc(B);
	K(:,:,2) = A*K(:,:,1);
	alpha(1,:) = dot(K(:,:,2),K(:,:,1));
	for i = 1:iter-1
		K(:,:,i+1) = K(:,:,i+1) - bsxfun(@times,K(:,:,i),alpha(i,:));
		beta(i+1,:) = sqrt(dot(K(:,:,i+1),K(:,:,i+1)));
		K(:,:,i+1) = bsxfun(@times,K(:,:,i+1),1./beta(i+1,:));
		K(:,:,i+2) = A*K(:,:,i+1) - bsxfun(@times,K(:,:,i),beta(i+1,:));
		alpha(i+1,:) = dot(K(:,:,i+2),K(:,:,i+1));
	end

	% Allocate and Compute output
	f = zeros(n,d);
	for z=1:d
		T = spdiags([[beta(2:iter,z)',0]' alpha(:,z) beta(:,z)], -1:1, iter, iter);
		[U,S] = eig(full(T));
		fS = diag(arrayfun(func, (diag(S))));
		xall = norm(B(:,z))*reshape(K(:,z,1:iter),n,iter)*U*fS*U';
		f(:,z) = xall(:,1);
	end

end
