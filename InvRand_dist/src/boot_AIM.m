function M0 = boot_AIM(A,options)
options = set_quNac_standard_options(A,options);
options.name = 'AIM';
% options.error_function = 
options.flopsperiter = 7*options.n^2 *options.p +2*options.n*options.p^2 +options.p^3+3*options.n^2+options.n;
M0 = eye(size(A));    %eye(size(A));   zeros
assignin('caller', 'options', options);
end