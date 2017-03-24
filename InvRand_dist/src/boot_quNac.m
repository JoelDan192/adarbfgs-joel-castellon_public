function M0 = boot_quNac(A,options)

options = set_quNac_standard_options(A,options);
%options.name = ['quNac\_H\_sample\_size\_' num2str(options.p) ];
options.name = 'ssBFGS';
% options.error_function = 
options.flopsperiter = 7*options.n^2 *options.p +2*options.n*options.p^2 +options.p^3+3*options.n^2+options.n;
M0 = eye(size(A));  
assignin('caller', 'options', options);
end