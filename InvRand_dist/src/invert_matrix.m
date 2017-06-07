% A wrapper function for testing and timing iterative methods for
% inverting a matrix - 2015 - Robert M. Gower
% InvRand Copyright (C) 2016, Robert Gower 
function [M, output, eigs] = invert_matrix(Prob, iter_func, boot_func, iter,test,options )
tic;
M = boot_func(Prob.A,options);  %M = full(M);

Ident = speye(size(Prob.A));
%Asqr = A^(1/2);
if(options.factored)
    initial_error= norm(Ident -M'*Prob.A*M,'fro');
    %initial_error= norm(Ident -M*(M')*A,'fro');
else
    initial_error= norm(Ident -M*Prob.A,'tr');
end
errors(1) = 1;
if(~isfield(options,'max_time'))
    options.max_time = 60*10; % 10 minutes
end

lo_e_ax = [];
mi_e_ax = [];
up_e_ax = [];
spec1_z = [];
spec2_z = [];
max_ps = [];
min_ps = [];
all_eigs = [];
[probs D] = complete_dicrete_sampling(Prob.A,M,options.p,'tr');
times(1) = toc;
%options.Sk = randsample(options.n,options.p);
for i = 1:iter
    if (mod(i,30)==1)
        display(options.name);
        fprintf('-------------------\n'); fprintf('It | Error%%   |  Time   \n'); fprintf('-------------------\n');
    end
    fprintf('%3.0d  | %3.2f%%   |  %3.4f \n',i,100*errors(i),times(i) );
    tic;
    M = iter_func(Prob.A, M,options,probs);
    times(i+1)= times(i) +  toc;
    
    [probs D] = complete_dicrete_sampling(Prob.A,M,options.p,'tr');
    r = size(M,2)/options.p;
    min_ps(end+1)=min(probs);
    max_ps(end+1)=max(probs);
    tmp = mat2cell(D,repmat(options.p,r,1),size(D,2));
    
    D = blkdiag(tmp{:});
    Xk = M*M';
    AX = ((Prob.A)^(1/2))*Xk*((Prob.A)^(1/2));
    options1.tol = 1e-3;
    Z = ((Prob.A)^(1/2))*M*(D^2)*(M')*((Prob.A)^(1/2));
    e_z = real(eig(Z));
    e_ax = real(eig(AX));
    all_eigs = [all_eigs; sort(e_ax')];
    inv_zmd = sqrtm(inv(Prob.A))*inv(M')*inv(M)*sqrtm(inv(Prob.A));
    
    band = 1e-3;
    lo_e_ax(end+1) = sum((1-e_ax)>band);%<1
    mi_e_ax(end+1) = sum(abs(e_ax-1)< band);%1
    up_e_ax(end+1) = sum((e_ax-1)>band); %>1
    
    %lo_e_ax(end+1) = real(min(diag(Z)));
    %e_ax_i = wrev(eig(inv_zmd));
    %e_prod = e_z.*e_ax_i;
    %lo_e_ax(end+1) = (trace(D^2)-sum(e_prod(2:end)))/e_ax_i(1);
    
    %lo_e_ax(end+1) = trace(D^2)/trace(inv_zmd);
    %mi_e_ax(end+1) = min(e_z);
    %up_e_ax(end+1) = min(e_ax)/sum(e_ax);
    %up_e_ax(end+1) = min(e_ax)/max(e_ax);
    
    spec1_z(end+1) = max(e_z);
    spec2_z(end+1) = max_ps(end);
    
    
    if (~isempty(test))
        if(options.factored)
            errors(i+1)= norm(Ident -M'*Prob.A*M,'fro')/initial_error;
            %errors(i+1) = norm(Ident -M*(M')*A,'fro')/initial_error;
        else
            errors(i+1) = norm(Ident -M*Prob.A,'fro')/initial_error;
        end
        if(errors(i+1)/100 < 10^(-4))
            output.fail =0;
            break;
        end
        if(isnan(errors(i+1)))
            output.fail =1;
            return;
        end
        if(isinf(errors(i+1)))
            output.fail =1;
            return;
        end
    end
    if(times(i+1) >options.max_time )
        output.fail ='times_up';
        break;
    end
end
if(i==iter) output.fail ='max_iter';  end %#ok<SEPEX>
if (~isempty(test)) fprintf('%3.0d  | %3.2f%%  |  %3.4f \n',i,100*errors(i+1),times(i+1) ); end %#ok<SEPEX>
clear M;
M =0;

lo_e_ax =[lo_e_ax, lo_e_ax(end)];
lo_e_ax = [0 lo_e_ax];

mi_e_ax =[mi_e_ax, mi_e_ax(end)];
mi_e_ax = [0 mi_e_ax];

up_e_ax =[up_e_ax, up_e_ax(end)];
up_e_ax = [0 up_e_ax];

spec1_z = [spec1_z, spec1_z(end)];
spec1_z = [spec1_z(1), spec1_z];

spec2_z = [spec2_z, spec2_z(end)];
spec2_z = [spec2_z(1), spec2_z];

min_ps = [min_ps min_ps(end)];
min_ps = [0 min_ps];

max_ps = [max_ps max_ps(end)];
max_ps = [1 max_ps];

output.flopsperiter = options.flopsperiter;
output.times = [ 0 times];
%output.errors = [1 errors];
%output.errors1 = 1-(nums./dens);
%output.errors2 = 1-ogs;
%output.errors3 = 1-lowb;
output.errors1 = lo_e_ax;
output.errors2 = mi_e_ax;
output.errors3 = up_e_ax;
eigs = all_eigs;
% uno=up_e_ax(end)
% dos=1/length(Prob.A)
%output.errors1=spec1_z;
%output.errors2=spec2_z;
%output.errors1 = e_ax;
%output.errors1=min_ps
%output.errors2=max_ps

output.name = options.name;
end
