function [c,info] = loewner_null_rec2_sym(p_c,p_r,tab)

%warning('off','backtrace')
c       = [];
nflop   = 0;
kk      = 0;
% Along p1
p2      = length(p_c{2});
W       = tab(1:length(p_c{1}),p2);
V       = tab(1+length(p_c{1}):end,p2);
LL_p1   = mlf.loewnerMatrix_sym({p_c{1}},{p_r{1}},W,V);
if any(~isreal(p_c{1})) % /!\ here you should check complex conjugation also
    warning('Enforce realness')
    J0  = sqrt(2)/2 * [1 -1i; 1 1i];
    T1  = [];
    for ii = 1:length(p_c{1})/2
        T1 = blkdiag(T1,J0);
    end
    c_p1 = T1*null(T1'*LL_p1*T1);
else
    c_p1 = null(LL_p1);
end
nflop   = nflop + size(LL_p1,1)*size(LL_p1,2)^2;
% % Check rank theorem
% if (min(size(c_p1)) > 1) || (isempty(c_p1))
%     %warning(['rank(LL_p1) < ' num2str(length(LL_p1)-1) ' (scaling is wrong)'])
%     warning('Go numeric (svd): null(LL_p1) issue')
%     [~,~,V] = svd(double(LL_p1)); 
%     c_p1    = V(:,end)/V(end,end);
% end
% Along p2
for p1 = 1:length(p_c{1})
    kk      = kk + 1;
    W       = tab(p1,1:length(p_c{2}));
    V       = tab(p1,1+length(p_c{2}):end);
    LLp2    = mlf.loewnerMatrix_sym({p_c{2}},{p_r{2}},W,V);
    c_it    = null(LLp2);
    nflop   = nflop + size(LLp2,1)*size(LLp2,2)^2;
    % Null-space
    c       = [c; c_it(:,1)*c_p1(kk,1)];
    % % Check rank theorem
    % if (min(size(c_pn{kk})) > 1) || (isempty(c_pn{kk}))
    %     %warning(['rank(LL_p2{' num2str(kk) '}) < ' num2str(length(LLp2)-1) ' (c_p1{' num2str(kk) '} is certainly wrong)'])
    %     warning(['Go numeric (svd): null(LL_p2{' num2str(kk) '}) issue'])
    %     [~,~,V]     = svd(double(LLp2)); 
    %     c_pn{kk}    = V(:,end)/V(end,end);
    % end
end
%
info.nflop  = nflop;

%warning('on','backtrace')
