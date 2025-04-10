function C = vec2mat(c,dim)

%%% Compute combinations
comb    = mlf.combinations(rand(dim));

%%% Build vector
N   = size(comb,1);
if isa(c,'sym')
    C   = sym('C',dim);
else
    C   = zeros(dim);
end
for ii = 1:N
    tmp         = num2cell(comb(ii,:));
    C(tmp{:})   = c(ii);
end
