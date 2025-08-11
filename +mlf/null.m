function [c,sig,km] = null(LL,METHOD)

    if nargin < 2
        METHOD = 'svd';
    end
    %METHOD = strcat(METHOD,'XXXXXXXXXXXXXXXX');
    %warning('off','backtrace')
    c   = [];
    sig = [];
    METHOD_ = METHOD;
    if str2double(METHOD(end)) >= 0
        METHOD_ = METHOD(1:end-1);
    end

    switch lower(METHOD_)
        case 'qr'
            [Q,~]       = qr(double(LL).');
            sig         = [];
            c           = Q(:,end);
            c           = normalize(c,str2double(METHOD(3:end)));
        case 'svd'
            [~,sig,V]   = svd(double(LL),'econ');
            sig         = diag(sig);
            sig         = sig/sig(1);
            c           = V(:,end);
            [c,km]      = normalize(c,str2double(METHOD(4:end)));
        case 'rsvd'
            [~,sig,V]   = svd_random(double(LL),[]);
            sig         = diag(sig);
            sig         = sig/sig(1);
            c           = V(:,end);
            c           = normalize(c,str2double(METHOD(5:end)));
        case 'null_symr'
            LL          = sym(LL);
            c           = null(LL,'r');
            %c           = c(:,end);
            %c           = normalize(c,str2double(METHOD(10:end)));
        case 'null_sym'
            LL          = sym(LL);
            c           = null(LL);
            c           = c(:,end);
            c           = normalize(c,str2double(METHOD(9:end)));
        case 'null'
            c           = null(LL);
            c           = c(:,end);
            sig         = [];
            c           = normalize(c,str2double(METHOD(5:end)));
        case 'mldivide'
            c_div       = LL(:,1:end-1)\LL(:,end);
            c           = [-c_div; 1];
            c           = normalize(c,str2double(METHOD(9:end)));        
        otherwise %isnumeric(METHOD)
            [U,S,V]     = svd(double(LL),'econ');
            sig         = diag(S);
            sig         = sig/sig(1);
            idx         = sig<METHOD;
            sig(idx)    = 0;
            sum(idx);
            fprintf('%d normalized singular value(s) < %1.2g -> to 0\n',sum(idx),METHOD)
            S           = diag(sig);
            [~,S,V]     = svd(double(U*S*V'),'econ');
            sig         = diag(S);
            sig         = sig/sig(1);
            k           = sum(sig>METHOD);
            c_svd       = V(:,k);
            c_svd       = c_svd/c_svd(end);
            c           = c_svd;
    end
    %warning('on','backtrace')
end

function [U,S,V] = svd_random(A,r)
    [m,n] = size(A);
    if isempty(r)
        r = n;
        p = n;
    else
        p = min(2*r,n);
    end
    X        = randn(n,p);
    Y        = A*X;
    W1       = orth(Y);
    B        = W1'*A;
    [W2,S,V] = svd(B,0);
    U        = W1*W2;
    r        = min(r,size(U,2));
    U        = U(:,1:r);
    S        = S(1:r,1:r);
    V        = V(:,1:r);
end

function [cn,km] = normalize(c,k)
    if k > 0
        cn      = c/c(k);
        km      = k;
    elseif k == 0
        [~,km]  = max(abs(c));
        cn      = c/c(km);
    else
        cn      = c/c(end);
        km      = size(c,1);
    end
end


