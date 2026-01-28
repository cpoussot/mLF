function B = basis_mon_inv(s,n)

for ii = 1:n+1
    B(ii) = s^-(ii-1);
end
