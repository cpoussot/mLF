function tab = eval1(H,Z,VERBOSE)

Ndigit = 32;

%%% Compute combinations
n = numel(Z);

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
    %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
    tmp_z_comma = regexprep(num2str(tmp_z,Ndigit), '\s*', ',');
    eval(['tab(ii) = H(' tmp_z_comma ');']);
    %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
end
