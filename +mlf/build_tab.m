function [tab,xpt,dim] = build_tab(H,ip)

n   = length(ip);
dim = zeros(1,n);
for ii = 1:n
    dim(ii) = length(ip{ii});
end
idx = mlf.combinations_dim(dim);
xpt = zeros(size(idx));
for ii = 1:n
    tmp         = ip{ii};
    xpt(:,ii)   = tmp(idx(:,ii));
end
tab = H(xpt);