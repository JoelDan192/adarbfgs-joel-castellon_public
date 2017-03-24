% main for tower matrices
cd  /home/s1065527/Dropbox/Matlab/stochastic_newton
%% main for testing convergence
n =6;    %select dimension
m = 6;
R = randn(n,m);   % random matrix
pinvR = pinv(R);
%% Schulz-Newton method
invRnorm = norm(pinvR,'fro');
test.invAnorm =invRnorm;
test.invA = pinvR;
[M, diffqu] = ShulzNewtonInverse(R,[],n, test);
plot(diffqu)
%% preparing row-weighted probability dist
sR  =size(R);
nb = zeros(1,sR(2));
for i = 1: sR(2)
nb(i) = norm(R(:,i))^2;
end
pnb= nb/sum(nb);
%Ek= 2*log(0.1/norm(pinvR,'fro')^2)*(1/(log(1-1/kA^2))) ;
size(gendist(pnb,m,p) )
size(randsample(m,p))
%% Broyden W = I
p = 1;
num_samp = ceil(30*n);% 30*(ceil(n/p))^(4);
H =eye(m,n);  % also works with  H = 0 or H =eye(n)
diffqu = zeros(num_samp,1);
for i=1:num_samp
S =  randsample(m,p); % gendist(pnb,p,1) ; %randsample(m,p);
%e = imerse_matrix(m,S);
%Re = R(:,S);
e = H(:,S);
Re = R*H(:,S);
%H=  H+ (e-H*A(:,S))*(A(:,S)')./na(S);
%H= H +(e-H*A(:,S))*((e'*A*e)^(-1))*(e');
H= H +(e-H*Re)*((Re'*Re)^(-1))*(Re');
diffqu(i) = 100*norm(pinvR -H,'fro')/norm(pinvR,'fro');
end
figure();
plot(diffqu)
%% Good Broyden (Inverse of direct with W =I)
p = 1;
num_samp = ceil(3);% 30*(ceil(n/p))^(4);
H =eye(m,n);  % also works with  H = 0 or H =eye(n)
B = eye(m,n); 
diffqu = zeros(num_samp,1);
for i=1:num_samp
S =  randsample(m,p); % gendist(pnb,p,1) ; %randsample(m,p); S = i;
e = imerse_matrix(m,S);
Re = R(:,S);
div = e'*H*(Re);
if (det(div) ~= 0)
       H= H +(e-H*Re)*((div)^(-1))*(e'*H);
       B = B+(Re-B*e)*(e'*e)^(-1)*e';
end
diffqu(i) = 100*norm(pinvR -H,'fro')/norm(pinvR,'fro');
directdiffqu(i) = 100*norm(R -B,'fro')/norm(R,'fro');
end
figure();
subplot(2,2,1)
plot(directdiffqu)
subplot(2,2,2)
plot(diffqu)