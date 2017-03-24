function [A, title] = sparse_symmetric_matlab(n,density,rc,type)

if(isempty(type))
    A = sprandsym(n,density,rc) ;  % sprandsym(n,density,rc,kind)
    title =['sprandsym-' num2str(n)  '-' num2str(density) '-' num2str(rc)  ];
else
    A = sprandsym(n,density,rc,type);
    title =['sprandsym-' num2str(n) '-' num2str(density) '-' num2str(rc) '-' type  ];
end
% Prob.b = rand(n,1);
% Prob.sol = Prob.A\Prob.b;
end