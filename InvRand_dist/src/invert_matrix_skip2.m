function [M, output] = invert_matrix(Prob, iter_func, boot_func, iter,test,options )
% A wrapper function for testing and timing iterative methods for
% inverting a matrix - 2015 - Robert M. Gower
tic;
M = boot_func(Prob.A,options); display(options.name); 
times(1) = toc;
timeac =0;
Ident = speye(size(Prob.A));
%Asqr = A^(1/2);
ecount = 1;
if(options.factored)
    initial_error= norm(Ident -M'*Prob.A*M,'fro');
    %initial_error= norm(Ident -M*(M')*A,'fro');
else
    initial_error= norm(Ident -M*Prob.A,'fro');
end
errors(1) = 1;
if(~isfield(options,'max_time'))
    options.max_time = 60*10; % 10 minutes
end
for i = 1:iter
    tic;
    M = iter_func(Prob.A,M,options);
    timeac = timeac+ toc;
skip =10;
if( mod(i,skip)==0)
    if (mod(i,30*skip)==1)
        fprintf('-------------------\n'); fprintf('It | Error%%   |  Time   \n'); fprintf('-------------------\n');
    end
    fprintf('%3.0d  | %3.2f%%   |  %3.4f \n',i,100*errors(ecount),times(ecount) );
    times(ecount+1)= times(ecount) +  timeac;
    timeac = 0;
    ecount = ecount+1;
    if (~isempty(test))
        if(options.factored)
            errors(ecount)= norm(Ident -M'*Prob.A*M,'fro')/initial_error;
            %errors(i+1) = norm(Ident -M*(M')*A,'fro')/initial_error;
        else
            errors(ecount) = norm(Ident -M*Prob.A,'fro')/initial_error;
        end
        if(errors(ecount)/100 < 10^(-4))
            break;
        end
    end
    if(times(ecount) >options.max_time )
        break;
    end
end % end mod if

end
clear M;
M =0;
if (~isempty(test)) fprintf('%3.0d  | %3.2f%%  |  %3.4f \n',i,100*errors(end),times(end) ); end %#ok<SEPEX>
output.flopsperiter = options.flopsperiter;
output.times = [ 0 times];
output.errors = [1 errors];
output.name = options.name;
end
