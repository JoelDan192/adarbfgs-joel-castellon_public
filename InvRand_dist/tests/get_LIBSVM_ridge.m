function Prob = get_LIBSVM_ridge( title,prefix,lambda )
% generate a ridge regression matrix from LIBSVM data.
% lambda         regression parameter, set to 1 for default
[~, A] = libsvmread([prefix title]);
A = A'*A; 
n = length(A);
if(isempty(lambda))
    Prob.A = A + eye(n);
else
    Prob.A = A + lambda*eye(n);
end
Prob.title = title;
Prob.n = n;
end

