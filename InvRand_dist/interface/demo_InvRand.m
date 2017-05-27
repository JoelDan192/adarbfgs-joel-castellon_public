%setup_InvRand
%% main for testing convergence
n = 600;    %select dimension
%A = randn(n,n);   % randn(n,n)
A=randn(n);
Prob.A=0.5*(A'+A); % M^T=M.
L=eig(Prob.A);
Nbins=60;
% [U,S,V] = svd(A);
% diag_s = diag(S)';
% cut_idx = round(0.7*length(diag_s));
% upper_s = diag_s(1:cut_idx-1);
% lower_s = diag_s(cut_idx:end);
% lower_s = lower_s/lower_s(1);
% upper_s = lower_s(1)+(lower_s(2)*wrev(1:length(upper_s)).^1.4);
% size(upper_s)
% size(lower_s)
% diag_s = [upper_s lower_s];
% S=diag(diag_s);
% A=U*S*V';
%A=U*diag(diag_s)*V';
Prob.A = (A')*A; % symmetric postive definite matrix
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
% %% MR Global self-conditioned method
% options.M0 = 'identity_proj';
% [M, MRoutput] = invert_matrix(Prob, @iter_MRGlobal, @boot_MRGlobal,iter,test, options );
% OUTPUTS = [ OUTPUTS ; MRoutput];
%% AdaRBFGS_cols
options.M0 = 1;
options.factored =1;
options.sample_method = 'cols';
[M, colsoutput1] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
OUTPUTS = [ OUTPUTS ; colsoutput1];

% %% AdaRBFGS_cols
% options.M0 = 1;
% options.factored =1;
% options.sample_method = 'cols';
% p_order = ceil((Prob.n)^(0.5));%min(3*ceil(sqrt(options.n)),ceil(options.n/2));
% rdivs = 1:Prob.n;
% rdivs = rdivs(rem(Prob.n,rdivs)==0);
% [res, idx_min] = min(abs(rdivs - p_order));
% options.p = rdivs(idx_min);
% [M, colsoutput1] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
% OUTPUTS = [ OUTPUTS ; colsoutput1];
% 
% %% AdaRBFGS_cols
% options.M0 = 1;
% options.factored =1;
% options.sample_method = 'cols';
% p_order = ceil((Prob.n)^(0.6));%min(3*ceil(sqrt(options.n)),ceil(options.n/2));
% rdivs = 1:Prob.n;
% rdivs = rdivs(rem(Prob.n,rdivs)==0);
% [res, idx_min] = min(abs(rdivs - p_order));
% options.p = rdivs(idx_min);
% [M, colsoutput1] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
% OUTPUTS = [ OUTPUTS ; colsoutput1];
%% AdaRBFGS_gauss
%options.M0 = 1;
%options.factored =1;
%options.sample_method = 'gauss';  
%[M, gaussoutput] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
%OUTPUTS = [ OUTPUTS ; gaussoutput];
% %% Schulz-Newton method
% options.M0=  'NewtonShulz';
% options.factored =0;
% [M, SNoutput] = invert_matrix(Prob, @iter_ShulzNewton, @boot_ShulzNewton,iter,test, options );
% if(~SNoutput.fail)
%     OUTPUTS = [ OUTPUTS ; SNoutput];
% end
%% plotting
%close all;
%h =  figure('visible','off');
%Prob.title = [ Prob.title ];
%plotdata = extract_plot_data(OUTPUTS,Prob,'flopsperiter');
%prettyPlot_plotdata(plotdata,options)
%quit
%%
%options = set_quNac_standard_options(Prob.A,options);
%figure
%hist(OUTPUTS{1}.errors1, Nbins);
%title('Eigenvalue distribution at convergence (tol=1e-5)')
%hist(L, Nbins)



% subplot(1,2,1)
% ax1=OUTPUTS{1}.times*1000;
% % % 
% plot(OUTPUTS{1}.errors2,'DisplayName','q=n^{1/3}');
% hold on
% plot(OUTPUTS{2}.errors2,'DisplayName','q=n^{1/2}');
% plot(OUTPUTS{3}.errors2,'DisplayName','q=n^{3/5}');
% xlabel('time (ms)')
% %title('p_i = \lambda_{max}(S_i^TAS_i)/\lambda_{max}(S^TAS)')
% set(gca, 'YScale', 'log')
% 
% % 
% hold off
% % xlabel('time (ms)')
% % legend('show','Location','southeast')
% % %ylim([0, 0.9e-1])
% % 
% % 
% subplot(1,2,2)
% %ax1=OUTPUTS{1}.flopsperiter*(1:size(OUTPUTS{1}.errors1,2));
% ax1=1:size(OUTPUTS{1}.errors1,2);
% plot(OUTPUTS{1}.errors2,'DisplayName','q=n^{1/3}');
% hold on
% plot(OUTPUTS{2}.errors2,'DisplayName','q=n^{1/2}');
% plot(OUTPUTS{3}.errors2,'DisplayName','q=n^{3/5}');
% xlabel('iterations')
% set(gca, 'YScale', 'log')
% % 
% hold off
% % xlabel('iterations')

%legend('show','Location','southeast')




subplot(1,2,1)
ax1=OUTPUTS{1}.times*1000;
% % 
plot(ax1,OUTPUTS{1}.errors1,'DisplayName','1-\rho_{lower}');
hold on
plot(ax1,OUTPUTS{1}.errors2,'DisplayName','1-\rho');
plot(ax1,OUTPUTS{1}.errors3,'DisplayName','1-\rho_{upper}');
title('p_i = tr(S_i^TAS_i)/tr(S^TAS)')
xlabel('time (ms)')
set(gca, 'YScale', 'log')

% 
hold off
% xlabel('time (ms)')
% legend('show','Location','southeast')
% %ylim([0, 0.9e-1])
% 

subplot(1,2,2)
ax1=OUTPUTS{1}.flopsperiter*(1:size(OUTPUTS{1}.errors1,2));
figure
ax1=1:size(OUTPUTS{1}.errors1,2);
plot(ax1,OUTPUTS{1}.errors1,'DisplayName','above 1');
hold on
plot(ax1,OUTPUTS{1}.errors2,'DisplayName','at 1');
plot(ax1,OUTPUTS{1}.errors3,'DisplayName','below 1');
title('Count of eigenvalues (n=600)')
xlabel('iterations')
set(gca, 'YScale', 'log')
% 
hold off
% % xlabel('iterations')

legend('show','Location','southeast')





