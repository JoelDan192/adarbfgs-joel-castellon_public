n=20;
Q =  orth(rand(n,n));
A = randn(n,n);
eig(A)
eig(Q'*A*Q)