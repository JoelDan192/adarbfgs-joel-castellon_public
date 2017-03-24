%% Numeric test lab

% Comparing inv against many system solves
p =10;
n = 1000;
repeat = 100;
for i =1:100
B = rand(p,p);
R= rand(n,p);
% inverting
tic;
(R)*inv(B);
cinv(i)= toc
tic;
%system solves
R/B;
cbar(i) = toc
end

mean(cinv)
mean(cbar)

subplot(2,2,1)
hist(cinv)
hold on 
subplot(2,2,2)
hist(cbar)
