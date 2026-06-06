function [s_gam_opt,s_del_opt] = compute_min_real_order(p_c)

maxPerm = 1e2;

n = length(p_c);
if n == 1
    s_gam_opt = 1;
    s_del_opt = [];
elseif n > 10
    n2          = ceil(n/2);
    s_gam_opt   = 1:n2;
    s_del_opt   = n2+1:n;
else
    for ii = 1:n
        ord(ii) = length(p_c{ii});
    end
    Mopt        = inf;
    ord_perm    = perms(ord);
    s_perm      = perms(1:n);
    if size(ord_perm,1) > maxPerm
        fprintf('%d configurations tested only (out of %d)\n',maxPerm,size(ord_perm,1))
        ord_perm    = ord_perm(1:maxPerm,:);
    end
    for ii = 1:size(ord_perm,1)
        for jj = 1:n-1
            set1    = 1:jj;
            set2    = setdiff(1:n,set1);
            s_gam   = ord_perm(ii,set1);
            s_del   = ord_perm(ii,set2);
            kappa   = prod(s_gam);
            ell     = prod(s_del);
            M       = 2*ell + kappa - 1;
            if M < Mopt
                Mopt        = M;
                s_gam_opt   = s_perm(ii,set1);
                s_del_opt   = s_perm(ii,set2);
            end
        end
    end
end
