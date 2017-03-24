function M = iter_AIM(A,M, options)
  S = randsample(options.n,options.p);
     R = M(:,S); 
%   Msqrt = M^(1/2); R=Msqrt(:,S);
 % Mchol = chol(M)';
 % R = Mchol(:,S);

%R = rand(options.n,options.p);
AR =A*R;
invA =  inv((R')*AR);
RinvA = R*invA;
RAMAR = (((AR')*M)*AR);
% M = M + RinvA*(R') - RinvA*RAMAR*(RinvA'); % MAR*(RinvA')+RinvA*((MAR')*Aproj);

%% Testing - Comparing
%eye(size(A))
M = M + RinvA*(R')*(A -A*M*A )*(RinvA*(R')); % MAR*(RinvA')+RinvA*((MAR')*Aproj);
end