function [r,info] = alg1(tab,p_c,p_r,opt)

%%% Option
ord_tol     = opt.ord_tol;
ord_obj     = opt.ord_obj;
ord_N       = opt.ord_N;
ord_show    = opt.ord_show;
min_data    = opt.data_min;
method      = opt.method;
METHOD_NULL = opt.method_null;

%%% Order detection & Points selection
n                   = length(p_c);
ok                  = mlf.check(p_c,p_r);
if ord_tol == -1
    for i = 1:n
        dim(i)      = length(p_c{i});
    end
    ord             = dim-1;
else
    ord             = mlf.compute_order(p_c,p_r,tab,ord_tol,ord_obj,ord_N,ord_show);
end
[pc,pr,W,V,tabr]    = mlf.points_selection(p_c,p_r,tab,ord,min_data);
w                   = mlf.mat2vec(W);

%%% Lagrange via Loewner
if strcmp(method,'full')
    LL      = mlf.loewnerMatrix(pc,pr,W,V);
    size(LL)
    c       = mlf.null(LL,METHOD_NULL);
    flop    = size(LL,1)*size(LL,2)^2;
elseif strcmp(method,'rec')
    [c,lag] = mlf.loewner_null_rec(pc,pr,tabr,METHOD_NULL);
    flop    = lag.nflop;
else
    error('Unknown method')
end

%%% Output
r           = @(x) mlf.eval_lagrangian(pc,w,c,x,false);
info.flop   = flop;
info.method = method;
info.pc     = pc;
info.w      = w;
info.c      = c;
info.ord    = ord;

% %%% Basis and X canonical lagrangian basis
% for ii = 1:length(p_c)
%     tmp     = sym(['s' num2str(ii)]);
%     B{ii}   = mlf.basis_lag(tmp,p_c{ii});
%     S{ii}   = mlf.X_lag(tmp,p_c{ii});
% end
% 
% %%% Realization minimality
% [s_gam,s_del]   = mlf.compute_min_real_order(p_c);
% opt.s_gam       = s_gam;
% opt.s_del       = s_del;
% 
% %%% Realization
% C               = mlf.vec2mat(c,size(W));
% [Hr,info_real]  = mlf.realization(C,W,S,opt);
% info_real