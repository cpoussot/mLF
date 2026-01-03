function ord = alg2_order_bnd(p_c,p_r,tab,tol,N)

%Ndigit  = 32;
%ord     = [];
%%% Compute combinations
n   = numel(p_c);
%col = parula(n+1);
%%% Frozen Loewner matrices
% Case 1D
if n == 1
    W           = tab(1:length(p_c{1}));
    V           = tab(1+length(p_c{1}):end);
    LL1         = mlf.loewnerMatrix({p_c{1}},{p_r{1}},W,V);
    % 
    [~,sig,~]   = svd(double(LL1));
    sig_n       = diag(sig)/sig(1);
    rank_s(1)   = sum(sig_n>tol);
else
    comb_tab    = mlf.combinations_dim(size(tab));
    for kk = 1:n
        lam_k   = p_c{kk};
        mu_k    = p_r{kk};
        rank_si = 0;
        % Combination for kk-th configuration
        comb_tab_kk         = comb_tab(:,setdiff(1:n,kk));
        % Remove same rows
        [comb_tab_kk,ia]    = unique(comb_tab_kk,'rows','legacy');
        Ncomb               = size(comb_tab_kk,1);
        if isempty(N)
            Nspace  = 1:Ncomb;
        else
            %Nspace  = unique(floor(linspace(1,Ncomb,N*Ncomb/100)));
            Nspace  = unique(floor(linspace(1,Ncomb,N)));
            %Nspace  = 1:N:Ncomb;
        end
        
        for ii = Nspace
            tmp = num2cell(comb_tab(ia(ii),:));
            for jj = 1:n
                if jj == kk
                    tmp_z(jj) = Inf;
                else
                    tmp_z(jj) = tmp{jj};
                end
            end
            %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
            tmp_z_comma = regexprep(num2str(tmp_z), '\s*', ',');
            tmp_z_comma = regexprep(tmp_z_comma, num2str(Inf), ':');
            eval(['WV = tab(' tmp_z_comma ');']);
            %%%% HOW TO ELEGANTLY EVALUATE H AT THESE POINTS
            WV          = WV(:);
            W           = WV(1:length(lam_k));
            V           = WV(1+length(lam_k):end);
            LL1         = mlf.loewnerMatrix({p_c{kk}},{p_r{kk}},W,V);
            % 
            [~,sig,~]   = svd(double(LL1));
            sig_n       = diag(sig)/sig(1);
            rank_s_kk   = sum(sig_n>tol);
            rank_si     = max(rank_si,rank_s_kk);
        end
        rank_s(kk)      = rank_si;
    end
end
ord = rank_s(:);
ord = ord.';

% % if order = number of data, then not enought data and reduce order of -1
% for ii = 1:n
%     if ord(ii) == length(p_c{ii})
%         ord(ii) = min(ord(ii),length(p_c{ii})-1);
%         warning(sprintf('Not enought data along the %d-th variable',ii))
%     end
% end

