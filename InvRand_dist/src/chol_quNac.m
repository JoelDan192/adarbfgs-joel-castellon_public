function Ln= chol_quNac(S,AS,L)
% Linv,LinvT,D
% Update the Cholseky factor L to Ln of a Block BFGS update
% tic;
%% Type I, Gratton,  Ilunga, Tshimanga
% Cholseky in Matlab R^TR = chol(A)
% 
% SAS = (S')*AS;
% invA =  inv(SAS);
% R = chol(invA);
% 
% Ln =  L - (S*invA)*((AS')*L);
% LinvS = L\S; 
% SMinvS = (LinvS')*LinvS;
% 
% G= chol(SMinvS);
% Ln = Ln + S*(R')*(inv(G'))*LinvS';
%% **** Type II,  use self conditioned simplification - nested expression - L^(-T) cancellation
SAS = (S')*AS;
invA =  inv(SAS);
R = chol(invA);

sL =size(L);
Rs =  pick_columns (sL(1),s,R');
Ln = L + S*(Rs- R'*R*((AS')*L));
%% Type III,  use self conditioned simplification
% SAS = (S')*AS;
% invA =  inv(SAS);
% R = chol(invA);
% 
% %Ln =  L - (S*invA)*((AS')*L);
% Ln =  L - (S*R')*R*((AS')*L);
% LinvS = L\S; 
% Ln = Ln + S*(R')*LinvS';


%% Type IV,  use self conditioned simplification - nested expression
% SAS = (S')*AS;
% invA =  inv(SAS);
% R = chol(invA);
% LinvS = L\S; 
% Ln = L + S*R'*(LinvS'- R*((AS')*L));

% toc
end
