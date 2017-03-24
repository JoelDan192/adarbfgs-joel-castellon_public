function Problem = generate_gaussian_kernal( n )
% generates an n X n Gaussianl kernal matrix

coord = 0.1*randn(n,1);
A = rbf(coord,1);
A = A + 0.1*eye(n,n);  
Problem.title =[ 'gaussKern-' num2str(n)];
Problem.cond  =cond(A);
Problem.A = A;
end

