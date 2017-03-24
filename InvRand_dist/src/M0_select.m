function M0 = M0_select(A, options)

if(isfield(options, 'M0'))
    if(strcmp(options.M0,'identity_proj'))
        M0 = eye(size(A))*trace(A)/trace(A*A');
    elseif(strcmp(options.M0,'AAT_proj'))
        B = A*A';
        M0 = A'*trace(B)/trace(B*B');
    elseif(strcmp(options.M0,'NewtonShulz'))
        if(issparse(A) )
            M0 = A'/(0.99*normest(A'*A));
        else
            M0 = A'/(0.99*norm(A'*A));
        end
    elseif(isnumeric(options.M0)) % direct assignment
        if(length(options.M0)==1)
            M0 = options.M0*eye(size(A));
        else
        M0 = options.M0;
        end
    end
else % just use identity
    M0= eye(size(A));
end


end