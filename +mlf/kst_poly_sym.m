function [H,PHI,W] = kst_poly_sym(ip,tab)
    %%% Generic 
    %info = [];
    n = length(ip);
    for l = 1:n
        k(l) = numel(ip{l});
        eval(['syms x' num2str(l) ';'])
    end
    %
    k1 = [1 k];
    k2 = [k 1];
    %%% Lagrangian basis
    for l = 1:n
        eval(['xl = x' num2str(l) ';'])
        for j = 1:k(l)
            Ls{l}(j,1) = ellsj(xl,ip{l},j);
        end
        poly        = repmat(Ls{l},prod(k1(1:l)),1);
        PHI(:,l)    = kron(poly,ones(prod(k2(l+1:n)),1));
    end
    %%% W idx
    comb    = mlf.combinations_dim(k);
    N       = size(comb,1);
    for l = 1:n
        k0          = prod(k(l+1:n));
        w_idx{l}    = k0:k0:N;
        k0n         = prod(k(l:n));
        wn_idx{l}   = kron(k0n:k0n:N,ones(1,k(l)));
    end
    wn_idx{1}   = kron(1,ones(1,k(1)));
    %%% W normalized
    W = sym(zeros(N,n));
    for l = 1:n
        nj  = size(w_idx{l},2);
        w   = zeros(nj,1);
        for j = 1:nj
            tmp     = comb(w_idx{l},:);
            tmp     = num2cell(tmp(j,:));
            w(j,1)  = tab(tmp{:});
            % normalization
            if l ~= 1
                tmp     = comb(wn_idx{l},:);
                tmp     = num2cell(tmp(j,:));
                wsn     = tab(tmp{:});
                w(j,1)  = w(j,1)/wsn;
            end
        end
        W(:,l) = kron(w,ones(prod(k2(l+1:n)),1));
        %eval(['info.Phi_x'  num2str(l) ' = sym(w);'])
    end
    %%% Reconstruction KST
    H = simplify(sum(prod(W.*PHI,2)));

    %%% 
    
end

%%% 
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
end