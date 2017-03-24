%% Testing AIM inverse

%% AIM - Approximate Inverse
options.factored =0;
%options.sample_method = 'factored_columns';  % factored_gaussian ,  factored_columns
[M, qNoutput] = invert_matrix(A, @iter_AIM, @boot_AIM,iter,test, options );
OUTPUTS = [ OUTPUTS ; qNoutput];