function [out] = solver_fpn(S, lmbd, opts)
% Estimate precision matrices with nonnegative partial correlations 
% using fast projected Newton-like algorithm

% Inputs:
%    'S'    : sample covariance matrix
%    'lmbd' : regularization parameter (a scalar for L1 norm; a matrix for weighted L1 norm) 
%    'opts' : stopping criterion options, passed as a structure
%             possible field name values:
%                'max_iter'  : maximum number of iterations  
%                'max_time'  : maximum time 
%                'tol'       : tolerance, the algorithm stops if ||X_k+1 - X_k||_F / ||X_k||_F < tol   
%                'beta'      : parameter used in backtracking line search
%                'X_opt'     : exists if the optimal solution is available for computing the objective function error        
%                'edge'      : exists if a disconnectivity set is imposed
%                'display'   : "1" if display the iteration results, "0" otherwise.
%
%
% Outputs:
%     'out' : outputs as a structure
%             possible field name values:
%                'time'        : run time cost before algorithm stopped
%                'X_est'       : the estimated precision matrix
%                'objective'   : the objective function value when algorithm stopped
%                'obj_itr'     : store the objective function value for each iteration
%                'time_itr'    : store the cpu time for each iteration
%                'iterate'     : the number of iterations cost before algorithm stopped 
%                'relobj_itr'  : exists if 'X_opt' exists, store the relative error 
%                               of the objective function value for each iteration, i.e., |fk - fopt|/|fopt|     
%                'relerr_iter' : exists if 'X_opt' exists, store the relative error 
%                               of each iteration, i.e., ||X_k - X_opt||_F / ||X_opt||_F 
%                'converge'    : "1": algorithm converges; "0" algorithm does not converge.

%% Initialize CPU time    
t0 = tic;

%%
p  = size(S,1);
iter = 0;
delta = 1e-15;
alpha = 0.5;

%%
if size(lmbd, 1)>1
   Lamb = lmbd;       % weighted L1 norm 
else
   Lamb = lmbd * (ones(p,p) - eye(p));  % L1 norm
end

%% set default parameters
if nargin < 3
   opts.max_iter = 1e3;
   opts.max_time = 1e3;
   opts.tol      = 1e-6;
end

if isfield(opts,'edge')
    flag_edge = 1;
    edgeset   = find(~opts.edge); 
else
    flag_edge = 0;
end

if isfield(opts,'beta')
    beta = opts.beta; 
else
    beta = 0.5;
end

if isfield(opts,'X_opt')
    flag_opt    = 1;
    X_opt       = opts.X_opt;
    fopt        = objective_function(X_opt, S - Lamb);
    f_opt       = fopt.value;
    relobj_iter = [];
    relerr_iter = [];
else
    flag_opt    = 0;
end

if isfield(opts,'display')
    display     = opts.display;
else
    display     = 1;
end

%%
proj = @(T) min(0, T - diag(diag(T))) + diag(diag(T)); 

%% Initilization
X    = diag(1./diag(S));
Fcur = objective_function(X, S - Lamb);
if Fcur.flag
    X    = diag(1./(diag(S)+ones(p,1)*1e-3));
    Fcur = objective_function(X, S - Lamb);
end
objcur    = Fcur.value;
obj_iter  = objcur;
time_iter = toc(t0);

%%
if flag_opt
    rel_object   = abs(f_opt - objcur)/abs(f_opt);
    relobj_iter  = [relobj_iter; rel_object];
    rel_err      = norm(X_opt - X, 'fro')/norm(X_opt, 'fro');
    relerr_iter  = [relerr_iter; rel_err];
end

%%
grad  = @(T) -inv(T) + S - Lamb;
gradf = grad(X);

%%
check = check_stopping(opts, X, X + 1e16, iter, t0);

if check.stop
    X_new  = X;
    objnew = objcur;
end
%%
while ~check.stop
    step_size = 1;
    iter      = iter + 1;
%%    
    rstset = find( X - diag(diag(X)) - eye(p)*1e3 > -delta & gradf < 0);
    
    if flag_edge
       rstset = union(rstset,edgeset);
    end
           
    X_up            = X;
    X_up(rstset)    = 0;
        
    grady           = gradf;
    grady(rstset)   = 0;
    
    descent         = X_up*grady*X_up;    
    descent(rstset) = 0;
    
    %%
    Theta_f         = @(gamma) proj(X_up - gamma * descent);
    X_new           = Theta_f(step_size);
        
    %%                
     [A, flag_pd]   = chol(X_new,'lower');
     if ~flag_pd
        lg_det      = 2*sum(log(diag(A)));
        objnew      = -1*lg_det + dot(X_new(:), S(:) - Lamb(:));
     else
        objnew      = 1e8;
     end
     
    %% select step size
    gd  = abs(dot(gradf(:), descent(:)));
    gdI = abs(dot(gradf(rstset), X(rstset)));
    
    while  (objnew > objcur - alpha*step_size*gd - alpha * gdI) && step_size > eps 
         step_size = beta*step_size;
         X_new     = Theta_f(step_size);
            %%
         [A, flag_pd] = chol(X_new,'lower');
         if ~flag_pd
            lg_det    = 2*sum(log(diag(A)));
            objnew    = -1*lg_det + dot(X_new(:), S(:) - Lamb(:));
         else
            objnew    = 1e8;
         end
    end

%%
   check     = check_stopping(opts, X_new, X, iter, t0);
%%
   obj_iter  = [obj_iter; objnew];    
   time_iter = [time_iter; toc(t0)];
   
   if flag_opt
      rel_object  = abs(f_opt - objnew)/abs(f_opt);
      relobj_iter = [relobj_iter; rel_object];
      rel_err     = norm(X_opt - X_new, 'fro')/norm(X_opt, 'fro');
      relerr_iter = [relerr_iter; rel_err];
   end

   if display
       if mod(iter,5) == 0
          disp(['iter: ', num2str(iter), '   ', 'objective: ', num2str(obj_iter(end),11), '   ', 'cpu time: ', num2str(time_iter(end)) ]);
       end
   end
       
%% update variable
     gradf  = grad(X_new);
     X      = X_new;
     objcur = objnew;        
end

%%
run_time = toc(t0);

%% optput
if flag_opt
    out.relobj_itr = relobj_iter;
    out.relerr_itr = relerr_iter;
end

out.time       = run_time;
out.X_est      = X_new;
out.objective  = objnew; 
out.obj_itr    = obj_iter;
out.time_itr   = time_iter;
out.iterate    = iter;
out.converge.  = check.converge;
end
