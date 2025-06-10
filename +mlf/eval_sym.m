function tab = eval_sym(H,Z,VERBOSE)

Ndigit = 32;

%%% Compute combinations
n = numel(Z);
for ii = 1:n
    dim(ii) = length(Z{ii});
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
        fprintf('Evaluation %0.1f %%\n',min(floor(100*ii/N),100))
    end
    %
    tmp = num2cell(comb(ii,:));
    expr_left   = [];
    expr_right  = [];
    for jj = 1:n
        expr_left   = [expr_left 'subs('];
        tmp_z       = Z{jj}(tmp{jj});
        expr_right  = [expr_right 's' num2str(jj) ',' string(tmp_z) '),'];
        expr_right = strjoin(expr_right,'');
    end
    expr_right = expr_right{1}(1:end-1);
    
    %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
    eval(['tab(tmp{:}) = ' expr_left 'H,' expr_right ';']);
    %%subs(subs(H,s1,ss),t,tt);
    %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
end
