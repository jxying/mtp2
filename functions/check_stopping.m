% t0 -- cputime at the beginning
% options_general

function check = check_stopping(opts, x, y, iter, t0)

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
reason=[];

if toc(t0) > opts.max_time
    stop = 1;
    reason = 0; % reason = 0: maximum cpu time exceeded
end

if iter > opts.max_iter
    stop = 1;
    reason = 1; % reason = 1: maximum iterations exceeded
end

if norm(x - y,'fro')/norm(y,'fro')< opts.tol
   stop = 1;
   reason = 2; % reason = 2: Condition on successive iterations holds
end

check.stop = stop;
check.reason = reason;






