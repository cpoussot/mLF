function X = X_mon_inv(s,n)

for ii = 1:n
    X(ii,ii+1)  = s;
    X(ii,ii)    = -1;
end
X(n+1,1) = 1;