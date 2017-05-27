function options = set_quNac_standard_options(A,options)

if(~isfield(options,'n'))
    options.n = length(A);
end
if(~isfield(options,'p'))
    p_order = ceil((options.n)^(0.5));%min(3*ceil(sqrt(options.n)),ceil(options.n/2));
    rdivs = 1:options.n;
    rdivs = rdivs(rem(options.n,rdivs)==0);
    [res, idx_min] = min(abs(rdivs - p_order));
    options.p = rdivs(idx_min);
end

end
