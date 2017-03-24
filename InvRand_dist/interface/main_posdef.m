%% Setup paths
setup_InvRand;
%% main for testing convergence
n =500;    %select dimension
%% Uniform Random matrix
% A = randn(n,n);   % randn(n,n)
% Prob.A = (A')*A;      % symmetric postive definite matrix
% Prob.title =[ 'randn-' num2str(n)];
%% Guassian kernal
% Prob = get_gaussian_kernal(n);
%% get UF matrices  exe: 315 which is the Bai/mhdb416 matrix
% Prob = UFget(315)  
% Prob.title = Prob.name;
%% Harwell-Boeing format
% exe: Prob.title='1138_bus.rsa';
%[ Prob.A, ~] = hb_to_msm(Prob.title);  Prob.n=length(Prob.A);
%% Matrix Market format   % s3rmt3m3.mtx
% exe:  Prob.title= 's3rmt3m3.mtx';
%[Prob.A,Prob.n,~,~,~,~,~] = mmread(Prob.title);
%% sparse psd
%density = 1/log(n); rc = 1/n; type =1;
%[Prob.A,Prob.title] = sparse_symmetric_matlab(n,density,rc,1);
% A= full(A);
%% LIBSVM
% test ridge regression matrices formed using data from LIBSVM
title = 'gisette_scale';   %gisette_scale, protein, real-sim, letter.scale (easy), aloi (easy), covtype (rank def), connect-4 (rank def), mushrooms rank def, m7a (rank def), connect-4
prefix = '../libsvm/';
Prob = get_LIBSVM_ridge( title,prefix,[] )
%% Setup tests
%test_prob_setup;
Prob.n =length(Prob.A);
iter = 5*Prob.n;
options =[];  options.factored =0;
options.max_time = 1000;
OUTPUTS ={};
test=1;
options.M0 = 1;%;'identity_proj'; %% NewtonShulz % options.M0 = 10^(-5);%
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
options.M0= 'NewtonShulz';
[M, SNoutput] = invert_matrix(Prob, @iter_ShulzNewton, @boot_ShulzNewton,iter,test, options );
if(~SNoutput.fail)
    OUTPUTS = [ OUTPUTS ; SNoutput];
end
%% plotting
close all;
h =  figure('visible','off');
Prob.title = [ Prob.title '_M0_' options.M0  ];
plotdata = extract_plot_data(OUTPUTS,Prob,'flopsperiter');
prettyPlot_plotdata(plotdata)
%quit
%%



