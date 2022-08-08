function check = check_stopping(opts, x, y, iter, t0)
% check the stopping condition of the algorithm

if ~isfield(opts,'max_time')
    opts.max_time = 1e3;
end

if ~isfield(opts,'max_iter')
    opts.max_iter = 1e3;
end

if ~isfield(opts,'tol')
    opts.tol = 1e-6;
end

stop = 0;
converge = 0;

if toc(t0) > opts.max_time
    stop = 1;     % maximum cpu time exceeded  
end

if iter > opts.max_iter
    stop = 1;     % maximum iterations exceeded    
end

if norm(x - y,'fro')/norm(y,'fro')< opts.tol
   stop = 1;      % condition on successive iterations holds
   converge = 1; 
end

check.stop = stop;
check.converge = converge;
