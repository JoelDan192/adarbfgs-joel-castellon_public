function [pna]= row_weight_prob(A)

sA  =size(A);
for i = 1: sA(1)
na(i) = norm(A(i,:))^2;
end
pna= na/sum(na);

end