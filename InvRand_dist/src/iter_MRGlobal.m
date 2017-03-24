function M = iter_MRGlobal(A,M,options)
R = -A*M;
% Add a one to diagonal to form the residual
r= size(R,1);
R(1:r+1:end) = R(1:r+1:end) + 1.0;
Z = M*R;
AZ =A*Z;
a = trace(R'*AZ)/trace(AZ'*AZ);
M =M+a*Z;

%flopcount:   2*nnz(A)*n +3*n +n^3 +4*n^2
end