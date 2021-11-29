function F = objective_function(T,S)
% compute the objective function value
[A, FLAG]=chol(T,'lower');
lg_det =2*sum(log(diag(A)));
fun = -1*lg_det + dot(T(:), S(:));
F.value = fun;
F.flag = FLAG;
end
