function M0 = boot_quNacPSB(A,options)

options = set_quNac_standard_options(A,options);
%options.name = ['quNacPSB\_H\_sample\_size\_' num2str(options.p) ];
options.name = 'RPSB';
assignin('caller', 'options', options);
M0 = eye(size(A));  
end