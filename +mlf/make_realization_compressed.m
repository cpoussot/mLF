function [Hr,info] = make_realization_compressed(info_real)

%%% Parameters
n   = numel(info_real.vars);
for ii = 1:n
    si(ii)  = sym(['s' num2str(ii)]);
end

Phi_s   = info_real.Phi_s;
B       = info_real.Gr;
C       = info_real.Wr;
N       = info_real.dim - size(info_real.BB,1);

%%% Compression
A11     = Phi_s(1:N,1:N);
A12     = Phi_s(1:N,N+1:end);
A21     = Phi_s(N+1:end,1:N);
A22     = Phi_s(N+1:end,N+1:end);
As      = (A22\A21);  
Phi_s   = A11-A12*As;
Gr      = B(1:size(A11,1),:);
C2      = C(:,size(A11,2)+1:end);
Wr      = -C2*As;

%%% Compute transfer function
Phi_f   = matlabFunction(Phi_s,'Vars',si);
Wr_f    = matlabFunction(Wr,'Vars',si(info_real.s_del));
% vars in Phi
tmp1    = [];
for ii = 1:n
    tmp1 = [tmp1 ['s' num2str(ii) ',']];
end 
tmp1    = tmp1(1:end-1);
% vars in W
tmp2    = [];
for ii = info_real.s_del
    tmp2 = [tmp2 ['s' num2str(ii) ',']];
end
tmp2 = tmp2(1:end-1);
% Hr
eval(['Hr = @(' tmp1 ') Wr_f(' tmp2 ')*(double(Phi_f(' tmp1 '))\Gr);']);

%%% Some additional informations
info.vars   = si;
info.Phi_s  = Phi_s;
info.Phi    = Phi_f;
info.Gr     = Gr;
info.Wr_s   = Wr;
info.Wr     = Wr_f;
info.dim    = length(Gr);
info.dim1   = info_real.dim1;
info.dim2   = info_real.dim2;


%info.c      = c;
%info.w      = w;
%info.AA     = AA;
%info.BB     = BB;
%info.GAMMA  = GAMMA;
%info.DELTA  = DELTA;
%info.s_gam  = s_gam;
%info.s_del  = s_del;
