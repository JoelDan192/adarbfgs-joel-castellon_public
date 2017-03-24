function test_all_LIBSVM(problems,prefix)
%% problems is a cell array of file names, located in prefix
%% exe: To run the experiments in the paper [1], choose 
% problems =  {'aloi', 'protein','gisette_scale', 'real-sim' }
% prefix = './tests/LIBSVM/data/';
%% other examples: gisette_scale, protein, real-sim, letter.scale (easy), aloi (easy), covtype (rank def), connect-4 (rank def), mushrooms rank def, m7a (rank def), connect-4
for i =1:length(problems)
    %% load problem
    Prob = get_LIBSVM_ridge( problems{i},prefix,[] );
    %% setup
    Prob.n =length(Prob.A);
    iter = 5*Prob.n;
    options =[];  options.factored =0;
    options.max_time = 1000;
    OUTPUTS ={};
    test=1;
    options.M0 = 1;%  all methods use the identity as a initial guess 
   %% call all methods
    call_all_posdef_methods;
    %% plotting
    name_plot_save;
end


end