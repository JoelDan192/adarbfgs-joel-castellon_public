Prob.n =length(Prob.A); 
iter = 5*Prob.n;
if( tril( full(Prob.A) ) ==full(Prob.A) )
    Prob.A = Prob.A +Prob.A' -diag(diag(Prob.A));
    spy(Prob.A); 
end
