function M= quNac_PSB(R,A,M)
 % Is this numerically stable? Carry out tests
I_n= eye(size(A));
AR= A*R;
invAS = inv((R')*AR);
proj = (R*invAS*(R'));
%
projAS = AR*(inv(AR'*AR))*AR';
M = proj + (I_n -projAS)*(M-proj)*(I_n -projAS);
% M = proj - (I_n -projAS)*proj*(I_n -projAS);
% projAS*proj + 
%   + (I_n -projAS)*M*(I_n -projAS);

I= eye(size(A));
V= A*R;
D = inv(V'*V);
M1= (I-V*D*V')*M*(I-V*D*V') + V*D*R'+ R*D*V' - V*D*R'*V*D*V';
norm(M-M1)
end

