function out = generate_weighted_lambda(S, alpha, tau, epfs, opts)
% generate the weighted regularization parameter
% lamba_ij = alpha / ( | X_ij + epfs |^tau ), for any i \neq j;
% lambda_ij = 0, for any i = j;
% X is the maximum likelihood estimator

if nargin < 5
   opts.max_iter = 30;
   opts.max_time = 300;
   opts.tol = 1e-5;

   if nargin < 4
       epfs = 1e-3;
       if nargin < 3
          tau = 1;
       end
   end
else

    if ~isfield(opts,'tol')
        opts.tol = 1e-5;
    end
    if ~isfield(opts,'max_time')
        opts.max_time = 300;
    end
    if ~isfield(opts,'max_iter')
        opts.max_iter = 30;
    end

end

opts.display = 0;

%% compute the maximum likelihood estimator
lambda = 0;
mle = solver_fpn(S, lambda, opts);
initial = abs(mle.X_est);

H = 1./((initial + epfs).^tau);
H_off = H - diag(diag(H));
out = alpha*H_off;
end

