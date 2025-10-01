function tab = eval(H,Z,VERBOSE)

%%% Compute combinations
n = numel(Z);
for ii = 1:n
    dim(ii) = length(Z{ii});
end
if numel(dim) == 1
    tab     = zeros(1,dim);
    comb    = mlf.combinations(tab);
    comb    = comb(:,2);
else
    tab     = zeros(dim);
    comb    = mlf.combinations(tab);
end

%%% Build tensor
N   = size(comb,1);
for ii = 1:N
    % Show advances
    if VERBOSE 
        fprintf('Evaluation %0.1f %%\n',min(floor(100*ii/N),100))
    end
    %
    tmp = num2cell(comb(ii,:));
    for jj = 1:n
        tmp_z(jj) = Z{jj}(tmp{jj});
    end
    input_tmp_z = num2cell(tmp_z);
    tab(tmp{:}) = H(input_tmp_z{:});
end
