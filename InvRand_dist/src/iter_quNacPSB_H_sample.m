function M = iter_quNacPSB_H_sample(A,M, options)
  S = randsample(options.n,options.p);
  R = M(:,S);
  M= quNac_PSB(R,A,M);
end