function M0 = boot_MRGlobal(A,options)
options.name = 'MR';

M0 = M0_select(A, options);
n= length(A);
options.flopsperiter =   2*nnz(A)*n +3*n +n^3 +4*n^2
assignin('caller', 'options', options);
end