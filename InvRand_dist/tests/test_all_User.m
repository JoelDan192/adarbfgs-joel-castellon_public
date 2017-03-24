function test_all_User()
%% Tests three matrices that are the square of a uniformly random matrix. Their dimensions are 
ns = [10^4, 2*10^4, 5*10^4];

for i =1:length(ns)
    %% load problem
    n = ns(i);
    A = randn(n,n);   % randn(n,n)
    Prob.A = (A')*A;      % symmetric postive definite matrix
    Prob.title =[ 'randn-' num2str(n)];
    %% setup
    %test_prob_setup;
    Prob.n =length(Prob.A);
    iter = 5*Prob.n;
    options =[];  options.factored =0;
    options.max_time = 1000;
    OUTPUTS ={};
    test=1;
    options.M0 = 1;%'identity_proj';
    %% call all methods
    call_all_posdef_methods;
    %% plotting
    name_plot_save;
end


end