function tab = eval1_sym(H,Z,VERBOSE)

Ndigit = 32;
%%% Compute combinations
n = numel(Z);

%%% Symbolic
for ii = 1:n
    dim(ii) = length(Z{ii});
    eval(['s' num2str(ii) ' = sym(''s' num2str(ii) ''');']);
end
tab = sym('tab',[1,dim]);

%%% Build tensor
N   = length(Z{1});
for ii = 1:N
    % Show advances
    if VERBOSE 
        fprintf('Evaluation %d %%\n',min(floor(100*ii/N),100))
    end
    %
    tmp_z(1) = Z{1}(ii);
    for jj = 2:n
        tmp_z(jj) = Z{jj}(1);
    end
    tmp = tmp_z;
    %
    %tmp = num2cell(comb(ii,:));
    expr_left   = [];
    expr_right  = [];
    for jj = 1:n
        expr_left   = [expr_left 'subs('];
        tmp_z       = tmp(jj);%Z{jj}(tmp{jj});
        expr_right  = [expr_right 's' num2str(jj) ',' string(tmp_z) '),'];
        expr_right = strjoin(expr_right,'');
    end
    expr_right = expr_right{1}(1:end-1);

    %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
    eval(['tab(ii) = ' expr_left 'H,' expr_right ';']);
    %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS

    % %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
    % tmp_z_comma = regexprep(num2str(tmp_z,Ndigit), '\s*', ',');
    % eval(['tab(ii) = H(' tmp_z_comma ');']);
    % %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
end
