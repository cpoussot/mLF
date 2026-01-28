function B = basis_lag(s,lam)

for ii = 1:numel(lam)
    B(ii) = s-lam(ii);
end
