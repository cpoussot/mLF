function LL = loewnerMatrix(lam,mu,w,v)

%%% Compute combinations
comb_lam    = mlf.combinations(w);
comb_mu     = mlf.combinations(v);
if mean(comb_lam(:,1)) == 1
    comb_lam = comb_lam(:,2);
    comb_mu  = comb_mu(:,2);
end

%%% Compute Loewner lmatrix
n   = numel(lam);
LL  = zeros(size(comb_mu,1),size(comb_lam,1));
for co = 1:size(LL,2)
    tmp_lam = num2cell(comb_lam(co,:));
    for li = 1:size(LL,1)
        tmp_mu      = num2cell(comb_mu(li,:));
        num         = v(tmp_mu{:})-w(tmp_lam{:});
        den         = 1;
        for ii = 1:n
            den = den * (mu{ii}(tmp_mu{ii}) - lam{ii}(tmp_lam{ii}));
        end
        LL(li,co)   = num/den;
    end
end

