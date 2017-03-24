 %% MR Global self-conditioned method
    [M, MRoutput] = invert_matrix(Prob, @iter_MRGlobal, @boot_MRGlobal,iter,test, options );
    OUTPUTS = [ OUTPUTS ; MRoutput];
    %% CO-BFGS COnditioned ChOlesky BFGS
    options.factored =1;
    options.sample_method = 'cols';  % factored_gaussian ,  factored_columns
    [M, colsoutput] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
    OUTPUTS = [ OUTPUTS ; colsoutput];
    %% CO-BFGS COnditioned ChOlesky BFGS
    options.factored =1;
    options.sample_method = 'gauss';  % factored_gaussian ,  factored_columns
    [M, gaussoutput] = invert_matrix(Prob, @iter_COBFGS, @boot_COBFGS,iter,test, options );
    OUTPUTS = [ OUTPUTS ; gaussoutput];
    %% Schulz-Newton method
    options.factored =0;
    [M, SNoutput] = invert_matrix(Prob, @iter_ShulzNewton, @boot_ShulzNewton,iter,test, options );
    if(~SNoutput.fail)
        OUTPUTS = [ OUTPUTS ; SNoutput];
    end