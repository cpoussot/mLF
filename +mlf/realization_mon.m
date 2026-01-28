function [Hr,info] = realization_mon(C,W,S,opt)

%%% Options
n = numel(size(W));
if n == 2
    if min(size(W)) == 1
        n = 1;
    end
end

if isempty(opt)
    if n == 1
        s_gam = 1;
        s_del = [];
    else
        s_gam = 1:floor(n/2);
        s_del = floor(n/2)+1:n;
    end
else
    s_gam = opt.s_gam;
    s_del = opt.s_del;
end

%%% Parameters
for ii = 1:n
    si(ii)  = sym(['s' num2str(ii)]);
end

%%% Build GAMMA & DELTA
for ii = 1:length(s_gam)
    if ii == 1
        GAMMA   = S{s_gam(1)};
    else
        GAMMA   = kron(GAMMA,S{s_gam(ii)});
    end
end
if isempty(s_del)
    DELTA   = [];
end
for ii = 1:length(s_del)
    if ii == 1
        DELTA   = S{s_del(1)};
    else
        DELTA   = kron(DELTA,S{s_del(ii)});
    end
end

%%% Build AA & BB
Bs_conf = [s_gam s_del];
if n > 1
    W   = permute(W,Bs_conf);
    C   = permute(C,Bs_conf);
end
w       = mlf.mat2vec(W);
c       = mlf.mat2vec(C);
kappa   = length(GAMMA);
ell     = length(DELTA);
if ell == 0
    ell = 1;
    AA  = reshape(c,ell,kappa);
    BB  = [];
else
    AA  = reshape(c,ell,kappa);
    BB  = reshape(w,ell,kappa);%.*AA;
end

%%% Build realization
if n == 1
    Phi = [GAMMA(1:kappa-1,:); ...
           -AA              ];
    Gr  = [zeros(1,kappa-1) -1].';
    Wr  = w.';
else
    Phi = [GAMMA(1:kappa-1,:)  zeros(kappa-1,ell-1)  zeros(kappa-1,ell); ...
           AA                  (DELTA(1:ell-1,:)).'  zeros(ell,ell); ...
           BB                  zeros(ell,ell-1)      DELTA.'];
    Gr  = [zeros(kappa-1,1); ... 
           (DELTA(ell,:)).'; ... 
           zeros(ell,1)];
    Wr  = [zeros(1,kappa+2*ell-2) -1];
end
%
Gr = double(Gr);
Wr = double(Wr);

%%% Compute transfer function
%%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
Phi_f   = matlabFunction(Phi,'vars',si);
tmp     = [];
for ii = 1:n
    tmp = [tmp ['s' num2str(ii) ',']];
end
tmp = tmp(1:end-1);
eval(['Hr = @(' tmp ') Wr*(double(Phi_f(' tmp '))\Gr);']);
%%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS

%%% Some additional informations
info.vars   = si;
info.Phi_s  = Phi;
info.Phi    = Phi_f;
info.Gr     = Gr;
info.Wr     = Wr;
info.dim    = 2*ell + kappa - 1;
info.dim1   = kappa;
info.dim2   = ell;
info.c      = c;
info.w      = w;
info.AA     = AA;
info.BB     = BB;
info.GAMMA  = GAMMA;
info.DELTA  = DELTA;
info.s_gam  = s_gam;
info.s_del  = s_del;

