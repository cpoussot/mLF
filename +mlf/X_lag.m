function X = X_lag(s,lam)

    t_u     = unimodular(lam);
    X       = [Jlag(s-lam); t_u.'];

end

function Js = Jlag(xi)
    k   = length(xi);
    Js  = [ones(k-1,1)*xi(1) diag(-xi(2:end))];
end

function p = unimodular(lam)
    q = length(lam);
    for ii = 1:q
        tmp     = lam(setdiff(1:q,ii));
        tmp     = prod(lam(ii)-tmp);
        p(ii,1) = 1/prod(tmp);
    end
end