function tab = make_tab_sym(H,p_c,p_r,VERBOSE)

n       = length(p_c);
p{1}    = p_c;
p{2}    = p_r;
allDim  = cell(n,1);
for ii = 1:n
    allDim{ii}  = 1:2;
end
comb        = cell(1,n);
[comb{:}]   = ndgrid(allDim{:});
comb        = reshape(cat(n+1,comb{:}),[],n);
% 
N   = size(comb,1);
for ii = 1:size(comb,1)
    % Show advances
    if VERBOSE 
        fprintf('Tableau construction %0.1f %%\n',min(floor(100*ii/N),100))
    end
    tmp = num2cell(comb(ii,:));
    for jj = 1:n
        p_tmp{jj}   = p{tmp{jj}}{jj};
    end
    tab{tmp{:}} = mlf.eval_sym(H,p_tmp,false);
end
tab = cell2sym(tab);