% Lagrange basis
% (2) l_j^lam 

function lk = lsj(s,lam,j)

k   = length(lam);
lk  = 1;
for i = 1:k
    if i ~= j
        lk = lk*(s-lam(i));
    end
end