function [Hr,info_real] = make_realization_lag(p_c,W,C,opt)

if nargin < 4
    opt         = [];
    MAKE_MIN    = true;
else
    if isfield(opt,'s_gam') && isfield(opt,'s_del')
        MAKE_MIN = false;
    end
    if isempty(opt)
        MAKE_MIN = false;
    end
end
if length(p_c) == 1
    MAKE_MIN = false;
end

%%% Basis and X canonical lagrangian basis
for ii = 1:length(p_c)
    tmp     = sym(['s' num2str(ii)]);
    B{ii}   = (mlf.basis_lag(tmp,p_c{ii})).';
    S{ii}   = mlf.X_lag(tmp,p_c{ii});
end

%%% Realization minimality
if MAKE_MIN 
    [s_gam,s_del]   = mlf.compute_min_real_order(p_c);
    opt.s_gam       = s_gam;
    opt.s_del       = s_del;
end

%%% Realization
[Hr,info_real]  = mlf.realization_lag(C,W,S,opt);

%%% Add informations
info_real.S     = S;
info_real.basis = B;
