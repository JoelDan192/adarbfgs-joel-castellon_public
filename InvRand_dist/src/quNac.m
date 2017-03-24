function M= quNac(R,AR,M)
     
% tic;
%% Type I, Stable, full matrix products
% I_n= eye(size(A));
% invRAR = inv((R')*A*R);
% proj = (R*invRAR*(R'));
% projA = proj*A;
% M = proj + (I_n -projA)*M*(I_n -projA');
%% Type II, direct diagonal manip
invA =  inv((R')*AR);
RinvA = R*invA;
MAR = M*(AR);
Aproj = AR*(RinvA');
r= size(Aproj,1);
Aproj(1:r+1:end) = Aproj(1:r+1:end) - 1.0;
M = M + RinvA*(R') -MAR*(RinvA')+RinvA*((MAR')*Aproj);
%% Type III, save on products
% AR = A*R;
% invA =  inv((R')*AR);
% invAR = invA*(R');
% MAR = M*(AR);
% MAproj = MAR*invAR;
% M = M + R*invAR- MAproj'+(((AR*invAR)'-eye(size(A)))*MAproj) ;

%% Type IV, Gram-schimdt to invA, S. Gratton?

% toc
end


% % Type IV, Cheaper but numerically unstable!
%     Z = R*(inv(A(S,S)));
%     AS= A(:,S);
%     ZSAM = Z*((AS')*M);
%     M = M+Z*R'   -ZSAM -ZSAM' + ZSAM*AS*Z';

%% Type III, save on products, fast and unstable
% AR = A*R;
% invA =  inv((R')*AR);
% invAR = invA*(R');
% MAR = M*(AR);
% MAproj = MAR*invAR;
% M = M + R*invAR- MAproj- MAproj'+(invAR')*((AR')*MAproj);