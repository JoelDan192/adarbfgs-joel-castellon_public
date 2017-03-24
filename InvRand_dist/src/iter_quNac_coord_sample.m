function M = iter_quNac_coord_sample(A,M, options)
    S = randsample(n,p);
    R = imerse_matrix(n,S);
    Z = R*(inv(A(S,S)));
    ZSA = Z*(A(:,S)');
    M = Z*R' + (I_n -ZSA)*M*(I_n -ZSA');  
    
end