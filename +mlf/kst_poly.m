function [H,PHI,W] = kst_poly(ip,tab)
    % Generic 
    n = length(ip);
    for l = 1:n
        k(l) = numel(ip{l});
    end
    k1 = [1 k];
    k2 = [k 1];
    % PHI(x_l)
    PHI = cell(n,1);
    for l = 1:n
        poly0 = [];
        for j = 1:k(l)
            poly0 = [poly0; ellsj_num(ip{l},j)];
        end
        poly    = repmat(poly0,prod(k1(1:l)),1);
        poly    = kron(poly,ones(prod(k2(l+1:n)),1));
        PHI{l}  = @(x) poly*(x.^(k(l)-1:-1:0)).';
    end
    % W idx
    comb    = mlf.combinations_dim(k);
    N       = size(comb,1);
    for l = 1:n
        k0          = prod(k(l+1:n));
        w_idx{l}    = k0:k0:N;
        k0n         = prod(k(l:n));
        wn_idx{l}   = kron(k0n:k0n:N,ones(1,k(l)));
    end
    wn_idx{1}   = kron(1,ones(1,k(1)));
    % W normalized
    W = zeros(N,n);
    for l = 1:n
        kk = 0;
        for j = 1:size(w_idx{l},2)
            kk          = kk+1;
            tmp         = comb(w_idx{l},:);
            tmp         = num2cell(tmp(j,:));
            ws{l}(kk,1) = tab(tmp{:});
            % normalization
            if l ~= 1
                tmp         = comb(wn_idx{l},:);
                tmp         = num2cell(tmp(j,:));
                wsn         = tab(tmp{:});
                ws{l}(kk,1) = ws{l}(kk,1)/wsn;
            end
        end
        W(:,l) = kron(ws{l},ones(prod(k2(l+1:n)),1));
    end
    H = @(x) kst(W,PHI,x);
end

%%% Normalized Lagrange basis
function lk = ellsj_num(lam,j)
    k   = length(lam);
    lk  = 1;
    for i = 1:k
        if i ~= j
            lk = conv(lk,[1 -lam(i)])/(lam(j)-lam(i));
        end
    end
end

%%% KST 
function Hr = kst(W,PHI,x0)
    [N,n] = size(W);
    PHIx0 = zeros(N,n);
    for l = 1:numel(PHI)
        PHIx0(:,l) = PHI{l}(x0(l));
    end
    Hr = sum(prod(W,2).*prod(PHIx0,2));
end