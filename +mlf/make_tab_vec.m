function [tab,xpt,dim] = make_tab_vec(H,p_c,p_r)

ok  = mlf.check(p_c,p_r);
if ~ok
    error('Not disjoint left/right data sets')
end
n   = length(p_c);
for ii = 1:n
    dim(ii) = length(p_c{ii}) + length(p_r{ii});
end
idx = mlf.combinations_dim(dim);
xpt = zeros(size(idx));
for ii = 1:n
    
    tmp         = [p_c{ii} p_r{ii}].';
    xpt(:,ii)   = tmp(idx(:,ii));
end
tab = H(xpt);