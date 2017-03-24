function options = set_quNac_standard_options(A,options)

if(~isfield(options,'n'))
    options.n = length(A);
end
if(~isfield(options,'p'))
    options.p = ceil(sqrt(options.n));%min(3*ceil(sqrt(options.n)),ceil(options.n/2));
end

end
