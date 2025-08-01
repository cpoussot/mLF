function g = eval_lagrangian(z,w,c,zi,VERBOSE)

% for ii = 1:length(z)
%     z{ii}  = conj(z{ii});
% end

%%% Compute combinations
n = numel(z);
for ii = 1:n
    dim(ii) = length(z{ii});
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
        fprintf('Evaluation %d %%\n',min(floor(100*ii/N),100))
    end
    %
    tmp = num2cell(comb(ii,:));
    for jj = 1:n
        tmp_z(jj) = z{jj}(tmp{jj});
    end
    % 
    den = 1;
    for jj = 1:n
        den = den * (zi(jj)-tmp_z(jj));
    end
    %
    g_den(ii) = c(ii)/den;
    g_num(ii) = w(ii)*g_den(ii);
end
g = sum(g_num)/sum(g_den);


