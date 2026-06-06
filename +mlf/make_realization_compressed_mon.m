function [Hr,info_real] = make_realization_compressed_mon(p_c,w,c,opt)

%%% Nominal
[~,itmp0]   = mlf.make_realization_mon(p_c,w,c,opt);
[~,itmp]    = mlf.make_realization_compressed(itmp0);

%%% Parameters
n   = numel(p_c);
for ii = 1:n
    si(ii)  = sym(['s' num2str(ii)]);
end

%%% Compression for monomial form
N       = length(itmp0.S{1});
Phi_s   = itmp.Phi_s;
B       = itmp.Gr;
C       = itmp.Wr_s;
A11     = Phi_s(1:N,1:N);
A12     = Phi_s(1:N,N+1:end);
A21     = Phi_s(N+1:end,1:N);
A22     = Phi_s(N+1:end,N+1:end);
As      = (A22\A21);
Phi_s   = A11-A12*As;
Gr      = B(1:N,:);
C1      = C(:,1:N);
C2      = C(:,N+1:end);
Wr      = C1-C2*As;
%Hr      = matlabFunction((Wr*(Phi_s\Gr)));

%%%
Phi_f   = matlabFunction(Phi_s,'vars',si);
Wr_f    = matlabFunction(Wr,'vars',si);
tmp     = [];
for ii = 1:n
    tmp = [tmp ['s' num2str(ii) ',']];
end
tmp = tmp(1:end-1);
eval(['Hr = @(' tmp ') double(Wr_f(' tmp '))*(double(Phi_f(' tmp '))\Gr);']);
%
info_real.Phi   = Phi_f;
info_real.Wr    = Wr_f;

info_real.Phi_s = Phi_s;
info_real.Wr_s  = Wr;
info_real.Gr    = Gr;

