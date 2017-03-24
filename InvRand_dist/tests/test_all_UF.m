function test_all_UF(problems)
%% problems is an array of UF ids
%% exe: problems =  [873,  440  ,  359 	,  40	, 937, 315	 , 1401, 1611, 1331, 868, 1288];
% where
% 873 'mesh1em1' ,  440 'FIDAP/ex9'  ,  359 'Boeing/msc01440'	,  40	'HB/bcsstk18', 315	'Bai/mhdb416' , 1401 'Bates/Chem97ZtZ'
% 	1611	'Cylshell/s3rmt3m3' ,  1331 'MathWorks/Muu', 868
% 	'Pothen/bodyy4',  GHS_psdef 1288,  Nasa/nasa4704 760
for i =1:length(problems)
    %% load problem
    Prob = UFget(problems(i)) ;
    Prob.title = Prob.name;
    %% setup
    test_prob_setup;
    Prob.n =length(Prob.A);
    iter = 5*Prob.n;
    options =[];  options.factored =0;
    options.max_time = 1000;
    OUTPUTS ={};
    test=1;
    options.M0 = 1;
    %% call all methods
    call_all_posdef_methods;
    %% plotting
    name_plot_save;
end


end