% normalized Lagrange basis
% (2) ell_j^lam 

function lk = ellsj(s,lam,j)

k   = length(lam);
lk  = 1;
for i = 1:k
    if i ~= j
        lk = lk*(s-lam(i))/(lam(j)-lam(i));
    end
end