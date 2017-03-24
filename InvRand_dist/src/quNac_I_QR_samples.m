A = randn(n,n);   % randn(n,n)
%A = normrnd(0,1, n);
As = A + A';
%%
invAs = pinv(As);
invAsnorm = norm(invAs);
invAsnormfro = norm(invAs,'fro');
test.invAnorms =invAsnorm;
test.invAs = invAs;

%% preparing row-weighted probability dist
[pna]= row_weight_prob(A);
%% Convergence test of H quNac W=I, QR samples
% Convergence occurs with almost any random generator for S
p = 10;
H =0; % also works with  H = 0 or H =eye(n)
EIGS =[];

diffqu = zeros(n,1);
[Q,R] = qr(As);
S= inv(R);
for i=1:n/2
e =S(:,i); %randn(n,p); %mvnrnd(zeros(p,1),eye(p,p),n); mvnrnd(zeros(n,1),eye(n,n),p)' % mvnrnd(zeros(n,1),invA,p)'
H= quNac_PSB(e,As,H); % Calculate estimate update
%EIGS(:,i) = eig(H*As);
diffqu(i) = 100*norm(invAs -H,'fro')/invAsnormfro;
end
%subplot(2,2,p) % first subplot
plot(diffqu)
title(['Sample PSB' num2str(p) ' randn.'])
xlabel('updates')
ylabel('|invA -H| / |invA|')
sprintf('%f %% difference\n',diffqu )
% plot(real(EIGS'))