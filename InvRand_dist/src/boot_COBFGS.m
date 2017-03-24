function L0 = boot_COBFGS(A,options)

options = set_quNac_standard_options(A,options);
options.name = ['AdaRBFGS\_' options.sample_method ] ;
% options.error_function = 
%L0 =eye(size(A))*sqrt(trace(A)/trace(A*A'));
M0 = M0_select(A, options);
L0 = sqrt(M0);
%L0 =  eye(size(A));%speye(size(A)); %*(options.p/trace(A)); %A; %eye(size(A)); zeros(size(A));
n= options.n;
options.flopsperiter = 4*n^2 *options.p +5*n*options.p^2 +4*options.p^3 +n*options.p +n^2;

assignin('caller', 'options', options);
end
