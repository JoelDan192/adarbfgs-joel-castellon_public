function R =  pick_columns (n,s,T)

R = zeros(length(s),n);
for i =1:length(s)
   R(:,s(i)) =T(:,i) ;
end


end