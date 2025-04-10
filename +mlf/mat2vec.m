function c = mat2vec(C)

%%% Compute combinations
comb    = mlf.combinations(C);

%%% Build vector
N   = size(comb,1);
if isa(C,'sym')
    c   = sym('c',[N 1]);
else
    c   = zeros(N,1);
end
for ii = 1:N
    tmp     = num2cell(comb(ii,:));
    c(ii,1) = C(tmp{:});
end
