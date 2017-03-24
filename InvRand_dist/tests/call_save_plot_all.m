%% Call all methods, plot and save
%% Setup tests
test_prob_setup;
options.factored =0;
OUTPUTS ={};
test=1;
%% MR Global self-conditioned method
options.M0 ='identity';
[M, MRoutput] = invert_matrix(Prob, @iter_MRGlobal, @boot_MRGlobal,iter,test, options );
OUTPUTS = [ OUTPUTS ; MRoutput];
%% Schulz-Newton method
[M, SNoutput] = invert_matrix(Prob, @iter_ShulzNewton, @boot_ShulzNewton,iter,test, options );
OUTPUTS = [ OUTPUTS ; SNoutput];
%% CO-BFGS COnditioned ChOlesky BFGS 
% Choose to match time to inverse cost pXp and 9 nXn nXp matrix product  p^3 = 3pn^2 => p =3\sqrt(n)  
options.factored =1;
options.sample_method = 'cols';  % factored_gaussian ,  factored_columns
[M, colsoutput] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
OUTPUTS = [ OUTPUTS ; colsoutput];
%% CO-BFGS COnditioned ChOlesky BFGS 
% Choose to match time to inverse cost pXp and 9 nXn nXp matrix product  p^3 = 3pn^2 => p =3\sqrt(n)  
options.factored =1;
options.sample_method = 'gauss';  % factored_gaussian ,  factored_columns
[M, gaussoutput] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
OUTPUTS = [ OUTPUTS ; gaussoutput];
%% plotting 
close all;
h =  figure('visible','off');
plotdata = extract_plot_data(OUTPUTS,Prob,'flopsperiter');
prettyPlot_plotdata(plotdata)