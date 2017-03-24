%% Testing convergence proofs

n =6;
S= rand(n,n);

invS = inv(S);

for i=1:n
   ss(i) = (S(:,i)')*S(:,i) ; 
end

invS'*diag(ss)*invS

A'*A