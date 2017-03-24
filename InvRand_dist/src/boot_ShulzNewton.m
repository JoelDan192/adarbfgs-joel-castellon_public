function M0 = boot_ShulzNewton(A,options)

options.name = 'NewtonShulz';
%  M0 = eye(size(A))*trace(A)/trace(A*A');
M0 = M0_select(A, options);
n= length(A);
options.flopsperiter = nnz(A)*n +n^3 +2*n^2;
assignin('caller', 'options', options);
end