function [probs, D] = complete_dicrete_sampling(A,M,p,mode)
    r = size(M,2)/p;
    D = [];
    if(contains(mode,'tr'))
        normalization = trace(M'*A*M);
    else
        normalization = real(eigs(M'*A*M,1));
    end
    
    probs = [];
    for j=0:r-1
        from = 1+j*p;
        to = (j+1)*p;
        Sj = M(:,from:to);
        SAS = Sj'*A*Sj;
        if(contains(mode,'tr'))
            pj=sqrt(trace(SAS)/normalization);
        else
            pj=sqrt(real(eigs(SAS,1))/normalization);
        end
        probs(end+1)=pj^2;
        D = [D;pj*(sqrtm(inv(SAS)))];
    end

end