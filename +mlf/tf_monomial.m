function [g,info] = tf_monomial(p_c,w,c,VERBOSE)

% Ndigit  = 32;
info    = [];
% for ii = 1:length(p_c)
%     p_c{ii}  = conj(p_c{ii});
% end

%%% Compute combinations
n = numel(p_c);
for ii = 1:n
    dim(ii) = length(p_c{ii});
    eval(['s' num2str(ii) ' = sym(''s' num2str(ii) ''');']);
end
if numel(dim) == 1
    dim     = [dim 1];
end
%     tab     = sym('tab',[1,dim]);
%     comb    = mlf.combinations(tab);
%     comb    = comb(:,2);
% else
%     tab     = sym('tab',dim);%zeros(dim);
%     comb    = mlf.combinations(tab);
% end

%%% Projectors
for jj = 1:n
    % Show advances
    if VERBOSE 
        fprintf('Advancement %0.1f %%\n',min(floor(100*jj/n),100))
    end
    S = [];
    B = [];
    for ii = 1:numel(p_c{jj})
        tmp = setdiff(p_c{jj},p_c{jj}(ii));
        x   = poly(double(tmp));
        S   = [S x.'];
        eval(['B=[s' num2str(jj) '^(' num2str(ii-1) ') B];']);
    end
    Sj{jj}  = S;
    Bj{jj}  = B;
end

%%% Basis
B   = 1;
S   = 1;
for jj = 1:n
    S   = kron(S,Sj{jj});
    B   = kron(B,Bj{jj});
end
DEN     = (S*c);
NUM     = (S*(w.*c));
%
norm_max    = max(abs(DEN));
NUM         = NUM/norm_max;
DEN         = DEN/norm_max;
DEN(abs(DEN)<1e-12) = 0;
NUM(abs(NUM)<1e-12) = 0;
%
den     = B*DEN;%(S*c);
num     = B*NUM;%(S*(w.*c));
g       = num/den;

%%% 
info.Sj         = Sj;
info.S          = S;
info.basis      = B.';
info.num        = num;
info.den        = den;
info.num_coeff  = NUM;
info.den_coeff  = DEN;
info.C          = mlf.vec2mat(flipud(DEN),dim);
info.W          = mlf.vec2mat(flipud(NUM),dim);
% info.c          = (S*c);
% info.w          = (S*w);