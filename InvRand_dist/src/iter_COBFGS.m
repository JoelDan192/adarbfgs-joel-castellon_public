function L= iter_COBFGS(A,L, options)
%%  Updates Factored form using Sampled action of Factored form.
if(strcmp(options.sample_method,'cols'))  % Sample Columns of Factored form
    s = randsample(options.n,options.p);
    R = L(:,s);    AR =A*R;
    L=chol_quNac_L_cols(s,AR,L);
else % Gaussian samples 
    S = randn(options.n,options.p);
    R = L*S;    AR =A*R;
    L= chol_quNac_L_sample(S,R,AR,L);
end
%% Test against cholseky factorization
%  Mc = chol(M)';
%  R = Mc(:,s);

%  M= quNac(R,AR,M);
% options.L= chol_quNac(s,R,AR,options.L);
assignin('caller', 'options', options);
%
end