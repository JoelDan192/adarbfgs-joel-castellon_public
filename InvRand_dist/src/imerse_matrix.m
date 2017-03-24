function R =  imerse_matrix (n,S)

R = zeros(n,length(S));
for i =1:length(S)
   R(S(i),i) =1 ;
end


end