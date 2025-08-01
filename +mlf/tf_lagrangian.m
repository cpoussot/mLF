function [g,info] = tf_lagrangian(p_c,w,c,VERBOSE)

Ndigit  = 32;
info    = [];
for ii = 1:length(p_c)
    p_c{ii}  = conj(p_c{ii});
end

%%% Compute combinations
n = numel(p_c);
for ii = 1:n
    dim(ii) = length(p_c{ii});
    eval(['s' num2str(ii) ' = sym(''s' num2str(ii) ''');']);
end
if numel(dim) == 1
    tab     = sym('tab',[1,dim]);
    comb    = mlf.combinations(tab);
    comb    = comb(:,2);
else
    tab     = sym('tab',dim);%zeros(dim);
    comb    = mlf.combinations(tab);
end

%%% Build tensor
N   = size(comb,1);
for ii = 1:N
    % Show advances
    if VERBOSE 
        fprintf('Advancement %0.1f %%\n',min(floor(100*ii/N),100))
    end
    %
    tmp = num2cell(comb(ii,:));
    bas = 1;
    for jj = 1:n
        pts = p_c{jj}(tmp{jj});
        if isa(pts,'sym')
            eval(strcat(['bas = bas*(s' num2str(jj) ' - '],string(pts),');'))
        else
            eval(['bas = bas*(s' num2str(jj) ' - ' num2str(pts,Ndigit) ');'])
        end
    end
    B(ii)   = bas;
    num(ii) = w(ii)*c(ii)/bas;
    den(ii) = c(ii)/bas;
    NUM(ii) = w(ii)*c(ii);
    DEN(ii) = c(ii);
end
g               = sum(num)/sum(den);
info.basis      = B.';
info.num        = num.';
info.den        = den.';
info.num_coeff  = NUM.';
info.den_coeff  = DEN.';
% info.c          = c;
% info.w          = w;
