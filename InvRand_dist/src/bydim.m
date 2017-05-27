

options =[];  options.factored =0;
options.max_time = 1000;
OUTPUTS ={};
test=1;
options.M0 = 1;%;'identity_proj'; %% NewtonShulz

% options.M0 = 10^(-5);%

%% main for testing convergence
n = 2*1000;    %select dimension
A = randn(n,n);   % randn(n,n)
Prob.A = (A')*A;      % symmetric postive definite matrix
Prob.title =[ 'randn-' num2str(n)];
%% Setup tests
%test_prob_setup;
Prob.n =length(Prob.A);
iter = 5*Prob.n;
options.n=length(Prob.A);

%% AdaRBFGS_cols
options.M0 = 1;
options.factored =1;
options.p=ceil(log(options.n));
options.sample_method = 'gauss_{logn}';  
[M, colsoutput1] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter1,test, options );
OUTPUTS = [ OUTPUTS ; colsoutput1];
%% AdaRBFGS_cols
options.M0 = 1;
options.factored =1;
options.p=ceil(log(options.n)^2);
options.sample_method = 'gauss_{log2n}';  
[M, colsoutput1] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter1,test, options );
OUTPUTS = [ OUTPUTS ; colsoutput1];

%% AdaRBFGS_cols
options.M0 = 1;
options.factored =1;
options.p=ceil(options.n^(1/4));
options.sample_method = 'gauss_{0.25}';  
[M, colsoutput2] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter2,test, options );
OUTPUTS = [ OUTPUTS ; colsoutput2];
%% AdaRBFGS_cols
options.M0 = 1;
options.factored =1;
options.p=ceil(options.n^(1/2));
options.sample_method = 'gauss_{0.5}';  
[M, colsoutput2] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter2,test, options );
OUTPUTS = [ OUTPUTS ; colsoutput2];
%% plotting
close all;
h =  figure('visible','off');
Prob.title = [ Prob.title ];
plotdata = extract_plot_data(OUTPUTS,Prob,'flopsperiter');
prettyPlot_plotdata(plotdata,options)
%quit