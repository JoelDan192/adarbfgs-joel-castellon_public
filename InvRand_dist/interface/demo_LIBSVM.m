%% demo that inverts the Hessian of ridge regression problem using the data set "vowel" from LIBSVM
% InvRand Copyright (C) 2016, Robert Gower 
%% LIBSVM
% test ridge regression matrices formed using data from LIBSVM
title = 'vowel';   %gisette_scale, protein, real-sim, letter.scale (easy), aloi (easy), covtype (rank def), connect-4 (rank def), mushrooms rank def, m7a (rank def), connect-4
prefix = './tests/LIBSVM/data/';
Prob = get_LIBSVM_ridge( title,prefix,[] )
%% Setup tests
%test_prob_setup;
Prob.n =length(Prob.A);
iter = 5*Prob.n;
options =[];  options.factored =0;
options.max_time = 1000;
OUTPUTS ={};
test=1;
options.M0 = 1; %;'identity_proj'; %% NewtonShulz
%% MR Global self-conditioned method
options.M0 = 'identity_proj';
[M, MRoutput] = invert_matrix(Prob, @iter_MRGlobal, @boot_MRGlobal,iter,test, options );
OUTPUTS = [ OUTPUTS ; MRoutput];
%% AdaRBFGS_cols
options.M0 = 1;
options.factored =1;
options.sample_method = 'cols';  
[M, colsoutput] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
OUTPUTS = [ OUTPUTS ; colsoutput];
%% AdaRBFGS_gauss
options.M0 = 1;
options.factored =1;
options.sample_method = 'gauss';  
[M, gaussoutput] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
OUTPUTS = [ OUTPUTS ; gaussoutput];
%% Schulz-Newton method
options.M0=  'NewtonShulz';
options.factored =0;
[M, SNoutput] = invert_matrix(Prob, @iter_ShulzNewton, @boot_ShulzNewton,iter,test, options );
if(~SNoutput.fail)
    OUTPUTS = [ OUTPUTS ; SNoutput];
end
%% plotting
close all;
h =  figure('visible','off');
Prob.title = [ Prob.title ];
plotdata = extract_plot_data(OUTPUTS,Prob,'flopsperiter');
prettyPlot_plotdata(plotdata,options)


