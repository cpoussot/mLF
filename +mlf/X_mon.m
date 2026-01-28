function X = X_mon(s,n)

for ii = 1:n
    X(ii,ii)    = s;
    X(ii,ii+1)  = -1;
end
X(n+1,1) = 1;