function [c_null,info] = loewner_null_recNH_sym(p_c,p_r,H,Nc)

Ndigit  = 32;
c_null  = [];
nflop   = 0;%nflop0;%0;
kk      = 0;
n       = length(p_c);
%%% Create symbolic variables
for ii = 1:n
    eval(['s' num2str(ii) ' = sym(''s' num2str(ii) ''');']);
end
%fprintf('Iteration %d\n',n)
% Along p(1)
for ii = 2:n
    p(ii-1) = p_c{ii}(end);
end
%p_comma = regexprep(num2str(p,Ndigit), '\s*', ',');
p       = cellstr(string(p));
p_comma = [];
for ii = 1:numel(p)
    p_comma = [p_comma p{ii} ','];
end
p_comma = p_comma(1:end-1);
eval(['W = mlf.eval1_sym(H,{p_c{1},' p_comma '},false);'])
eval(['V = mlf.eval1_sym(H,{p_r{1},' p_comma '},false);'])
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
param_str = [];
for ii = 1:n-1
    param_str = [param_str 's' num2str(ii) ','];
end
param_str = param_str(1:end-1);
%
eval(['info.c' num2str(length(col)+1) '=c_p1;'])
eval(['info.w' num2str(length(col)+1) '=W(:);'])
%eval(['info.d' num2str(length(col)+1) '=W(1);'])
eval(['info.d' num2str(length(col)+1) '=W(end);'])
for p1 = 1:length(p_c{1})
    kk  = kk + 1;
    %eval(['H_p1 = @(' param_str ')H(p_c{1}(p1),' param_str ');']);
    H_p1 = subs(H,s1,p_c{1}(p1));
    for ll = 1:length(col)
        H_p1 = subs(H_p1,['s' num2str(ll+1)],['s' num2str(ll)]);
    end
    if length(col) > 2 
        [c_it,ii]   = mlf.loewner_null_recNH_sym(col,row,H_p1,Nc);
        %eval(['info.c' num2str(length(col)) '=ii;'])
        %eval(['info.c' num2str(length(col)) '=ii_p1;'])
    else 
        [c_it,ii]   = mlf.loewner_null_rec2H_sym(col,row,H_p1);
    end
    nflop    = nflop + ii.nflop;
    c_null   = [c_null; c_it(:,1)*c_p1(kk,1)];
    if length(p_c) == numel(Nc)
        fprintf('Recusive null-space construction : %d %%\n',floor(100*length(c_null)/prod(Nc)))
    end
    ii_p1{p1}   = ii;
    %cit(:,p1)   = c_it;
    %sca(:,p1)   = W(end);
end

%
%eval(['info.cit' num2str(length(col)+1) '=cit;'])
eval(['info.i' num2str(length(col)) '=ii_p1;'])
info.nflop  = nflop;
