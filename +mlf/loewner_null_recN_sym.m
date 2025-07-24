function [c_null,info] = loewner_null_recN_sym(p_c,p_r,tab,Nc);%nflop0)

c_null  = [];
nflop   = 0;
kk      = 0;
n       = length(p_c);
% Along p(1)
for ii = 2:n
    p(ii-1) = length(p_c{ii});
end
p_comma = regexprep(num2str(p), '\s*', ',');
eval(['W = tab(1:length(p_c{1}),' p_comma ');']);
eval(['V = tab(1+length(p_c{1}):end,' p_comma ');']);
LL_p1   = mlf.loewnerMatrix_sym({p_c{1}},{p_r{1}},W,V);
if any(~isreal(p_c{1})) && (sum(imag(p_c{1}))==0) % /!\ here you should check complex conjugation also
    warning('Enforce realness')
    J0  = sqrt(2)/2 * [1 -1i; 1 1i];
    T1  = [];
    ii  = 1;
    while ii <= length(p_c{1})
        if isreal(p_c{1}(ii))
            T1 = blkdiag(T1,1);
            ii = ii + 1;
        else
            T1 = blkdiag(T1,J0);
            ii = ii + 2;
        end
    end
    c_p1 = T1*null(T1'*LL_p1*T1);
else
    c_p1 = null(LL_p1);
end
%c_p1    = null(LL_p1);
nflop1  = size(LL_p1,1)*size(LL_p1,2)^2;
nflop   = nflop + nflop1;
% Check rank theorem
if min(size(c_p1)) > 1
    warning(['rank(LL_p1) < ' num2str(length(LL_p1)-1) ' (scaling is wrong)\n'])
end
% Along p(2...n)
col = {p_c{2:end}};
row = {p_r{2:end}};
%
tab_str1 = ones(1,n-1)*Inf;
tab_str1 = regexprep(num2str(tab_str1), '\s*', ',');
tab_str1 = regexprep(tab_str1, 'Inf', ':');
%
tab_str2 = size(tab);
tab_str2 = regexprep(num2str(tab_str2(2:end)), '\s*', ',');
%
eval(['info.c' num2str(length(col)+1) '=c_p1;'])
eval(['info.w' num2str(length(col)+1) '=W(:);'])
eval(['info.d' num2str(length(col)+1) '=W(end);'])
%
for p1 = 1:length(p_c{1})
    kk  = kk + 1;
    eval(['tab_p1 = reshape(tab(p1,' tab_str1 '),' tab_str2 ');'])
    if length(col) > 2 
        [c_it,ii]   = mlf.loewner_null_recN_sym(col,row,tab_p1,Nc);
    else 
        [c_it,ii]   = mlf.loewner_null_rec2_sym(col,row,tab_p1);
    end
    nflop   = nflop + ii.nflop;
    c_null  = [c_null; c_it(:,1)*c_p1(kk,1)];
    if length(p_c) == numel(Nc)
        fprintf('Recusive null-space construction : %d %%\n',floor(100*length(c_null)/prod(Nc)))
    end
    ii_p1{p1}   = ii;
end
%
eval(['info.i' num2str(length(col)) '=ii_p1;'])
info.nflop  = nflop;
