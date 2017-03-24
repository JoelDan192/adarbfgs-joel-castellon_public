%setup_InvRand
%% main for testing convergence
n = 1000;    %select dimension
A = randn(n,n);   % randn(n,n)
Prob.A = (A')*A;      % symmetric postive definite matrix
Prob.title =[ 'randn-' num2str(n)];
%% Setup tests
%test_prob_setup;
Prob.n =length(Prob.A);
iter = 5*Prob.n;
options =[];  options.factored =0;
options.max_time = 1000;
OUTPUTS ={};
test=1;
options.M0 = 1;%;'identity_proj'; %% NewtonShulz
% options.M0 = 10^(-5);%
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

%figure;
%subplot(2,2,1);
%plot(MRrate.nums, colsrate.nums, gaussrate.nums,SNrate.nums )
%subplot(2,2,2);
%plot(MRrate.dens, colsrate.dens, gaussrate.dens,SNrate.dens)
%subplot(2,2,3);
%plot(MRrate.nums./MRrate.dens, colsrate.nums./colsrate.dens, gaussrate.nums./gaussrate.dens, SNrate.nums./gaussrate.dens)
%subplot(2,2,4);
%plot(MRrate.ogs, colsrate.ogs, gaussrate.ogs, SNrate.ogs)
%quit
%%



