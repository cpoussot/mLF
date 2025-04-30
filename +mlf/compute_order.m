function ord = compute_order(p_c,p_r,tab,tol,ord_obj,N,SHOW)

Ndigit  = 32;
ord     = [];
%%% Compute combinations
n   = numel(p_c);
col = parula(n+1);
%%% Frozen Loewner matrices
if SHOW
    figure, hold on, grid on, axis tight
end
% Case 1D
if n == 1
    W           = tab(1:length(p_c{1}));
    V           = tab(1+length(p_c{1}):end);
    LL1         = mlf.loewnerMatrix({p_c{1}},{p_r{1}},W,V);
    % 
    [~,sig,~]   = svd(double(LL1));
    sig_n       = diag(sig)/sig(1);
    rank_s(1)   = sum(sig_n>tol);
    if SHOW
        h(1)    = plot(sig_n,'-s','Color',col(1,:));
        set(gca,'YScale','log'), drawnow
        xlabel('\# singular value','Interpreter','latex')
        ylabel('Normalized singular value','Interpreter','latex')
        legend_str{1} = '$\mathstrut^1 s$';
        ylim([1e-17 1])
    end
else
    comb_tab    = mlf.combinations(tab);
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
        %Nspace
        
        for ii = Nspace
           tmp = num2cell(comb_tab(ia(ii),:));
        % for ii = 1:size(comb_tab,1)
        %     tmp = num2cell(comb_tab(ii,:));
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
            if SHOW
                h(kk)       = plot(sig_n,'-s','Color',col(kk,:));
                %title(sprintf(['$s_' num2str(kk) '$']),'Interpreter','latex')
                title(sprintf(['$^' num2str(kk) 's$']),'Interpreter','latex')
                xlabel('\# singular value','Interpreter','latex')
                ylabel('Normalized singular value','Interpreter','latex')
                set(gca,'YScale','log'), drawnow
                ylim([1e-17 1])
            end
        end
        if SHOW
            legend_str{kk}  = ['$\mathstrut^' num2str(kk) 's$'];
        end
        rank_s(kk)      = rank_si;
    end
end

ord = rank_s(:);
if ~isempty(ord_obj)
    ord = min(ord_obj(:),rank_s(:));
end
ord = ord.';

% if order = number of data, then not enought data and reduce order of -1
for ii = 1:n
    if ord(ii) == length(p_c{ii})
        ord(ii) = min(ord(ii),length(p_c{ii})-1);
        warning(sprintf('Not enought data along the %d-th variable',ii))
    end
end

if SHOW
    hh          = get(gca); 
    ht          = plot(hh.XLim,tol*[1 1],'k--');
    legend([h, ht(1)],{legend_str{:}, '$\epsilon$'},'Location','best','interpreter','latex');
    title_str_left{1}   = 'd_1';
    title_str_right{1}  = '%d';
    for ii = 2:n
        title_str_left{ii}  = [',d_' num2str(ii)];
        title_str_right{ii} = ',%d';
    end
    %title_str = ['$(' title_str_left{:} ')=(' title_str_right{:} ')$'];
    %title(sprintf(title_str,rank_s),'Interpreter','latex')
    %title_str = ['$(' title_str_left{:} ')=(' title_str_right{:} ') \\rightarrow (' title_str_right{:} ')$'];
    set(gca,'TickLabelInterpreter','latex')
    title_str = ['$(' title_str_left{:} ')=(' title_str_right{:} ') \\rightarrow (' title_str_right{:} ') \\rightarrow N=%d$'];
    title(sprintf(title_str,rank_s,ord,prod(ord+1)),'Interpreter','latex')
    drawnow
end

