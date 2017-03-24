function M= Broyden(S,A,M)
     
AS= A*S;
M= M +(S-M*AS)*((AS'*AS)^(-1))*(AS');

end

