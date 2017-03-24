function [M,diffqu]= quNacInverse(p,A,M0,iter,test)

% Calculate update
sA= size(A); 
n= sA(1); 
I_p = eye(p);  

if (isempty(M0))
    M0 = zeros(size(A));
end
M = M0;
if (~isempty(test))
    diffqu = zeros(1,iter);
else 
    diffqu = [];
end
for i = 1:iter
    S = randsample(n,p);
    R = imerse_matrix(n,S);
    MR= M*R;
    MRR = MR*R';
    M = M +R*(I_p + (R')*MR)*(R') - MRR -MRR';
    if (~isempty(test))
        diffqu(i) = 100*norm(test.invA -M)/test.invAnorm;
    end
end

end

