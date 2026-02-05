function [r,info] = alg2(tab,p_c,p_r,opt)

if nargin < 4
    opt = struct();
end

%%% Options
n   = length(p_c);
% Tolerance
if ~isfield(opt,'tol')
    tol = 1e-10;
else
    tol = opt.tol;
end
% Method for null space computation
if ~isfield(opt,'method')
    method = 'rec';
else
    method = opt.method;
end
% Method for null space computation
if ~isfield(opt,'method_null')
    method_null = 'svd0';
else
    method_null = opt.method_null;
end
% "Estimate" if l-th variable is rational 
if ~isfield(opt,'max_itpl')
    tol1 = [1e-9 1e-12];
    for i = 1:numel(tol1)
        ord(i,:) = mlf.alg2_order_bnd(p_c,p_r,tab,tol1(i),2);
    end
    for i = 1:n
        ip_est(1,i) = length(p_c{i});
        if min(ord(:,i)) == max(ord(:,i))
            ip_est(1,i) = max(ord(:,i)) + 1;
        end
    end
    opt.max_itpl = ip_est;
end
% Initial interpolation partition (if any)
if ~isfield(opt,'ip_keep')
    ip_col = cell(1,n);
else
    ip_col = opt.ip_keep;
end
% Max iter
if ~isfield(opt,'max_iter')
    opt.max_iter = prod(cellfun(@length,p_c))-1;
end
% Loewner as small or tall
if ~isfield(opt,'data_min')
    data_min = true;
else
    data_min = opt.data_min;
end

%%% Start
for i = 1:n 
    ip{i} = [p_c{i}(:); p_r{i}(:)];
end
max_samples         = max(abs(tab),[],'all');
err_mat             = abs(tab-mean(tab,'all'));
[max_err,max_idx]   = max(err_mat,[],'all');
fprintf('Initial rel max error %d\n',max_err(1)/max_samples)

% do this such that at least one iteration is done
max_err         = Inf;
max_err_best    = max_err;
max_Idx         = cell(1,n);
jj              = 0;
flop            = 0;
%
tab_sz          = size(tab);
tab_vec         = mlf.mat2vec(tab);
%abs_tab_vec     = abs(tab_vec);
%
comb    = mlf.combinations_dim(tab_sz);
N       = size(comb,1);
ipt     = zeros(N,n);
for ii = 1:n
    ipt(:,ii) = ip{ii}(comb(:,ii));
end
no_improvement = 0; CONTINUE = true;
while (max_err > max_samples * tol) && (jj < opt.max_iter) && CONTINUE
    %
    clear pc pr w c
    jj = jj + 1;
    % 
    [max_Idx{:}] = ind2sub(size(tab),max_idx);

    %%% Check if maximum order has been reached
    add_itpl = cellfun(@(ip,mi)length(ip)<mi,ip_col,num2cell(opt.max_itpl));
    if ~any(add_itpl)
        fprintf('Reached maximum number of interpolation points \n')
        break
    end

    %%% Add interpolation points
    for i = 1:n
        % make sure to keep at least one sample in LS partition
        if add_itpl(i)
            ip_col{i} = unique([ip_col{i},max_Idx{i}]);
        end
    end

    %%% Set column and row interpolation points
    for i = 1:n
        ip_row{i}   = setdiff(1:numel(ip{i}),ip_col{i});
        if data_min
            tmp         = numel(ip_col{i});
            ip_row{i}   = ip_row{i}(1:tmp);
        end
        pc{i}       = ip{i}(ip_col{i});
        pr{i}       = ip{i}(ip_row{i});
        ip_all{i,1} = [pc{i}(:); pr{i}(:)];
        pr_sz(i)    = numel(pr{i});
    end

    %%% Re-order "tab" and "W"
    idx_str     = '[';
    idx_strW    = '[';
    idx_strV    = '[';
    for ii = 1:n
        tmp         = [ip_col{ii}(:); ip_row{ii}(:)].';
        idx_str     = [idx_str num2str(tmp) '],['];
        tmp         = [ip_col{ii}(:)].';
        idx_strW    = [idx_strW num2str(tmp) '],['];
        tmp_compl   = setdiff(1:numel(ip{ii}),tmp);
        idx_strV    = [idx_strV num2str(tmp_compl) '],['];
    end
    idx_str     = idx_str(1:end-2); 
    idx_strW    = idx_strW(1:end-2);
    idx_strV    = idx_strV(1:end-2);
    eval(['tab_it = tab(' idx_str ');']);
    eval(['W_it   = tab(' idx_strW ');']);
    eval(['V_it   = tab(' idx_strV ');']);
    
    %%% Compute null space
    if strcmp(method,'full') && jj > 1
        LL      = mlf.loewnerMatrix(pc,pr,W_it,V_it);
        c       = mlf.null(LL,method_null);
        flop    = flop + size(LL,1)*size(LL,2)^2;
    else%if strcmp(method,'rec')
        [c,lag] = mlf.loewner_null_rec(pc,pr,tab_it,method_null,[]);
        flop    = flop + lag.nflop;
    % else
    %     error('Unknown method')
    end
    
    %%% Compute error
    w   = mlf.mat2vec(W_it);
    hr  = zeros(N,1);
    for i = 1:N
        hr(i) = mlf.eval_lagrangian(pc,w,c,ipt(i,:),false);
    end
    err_mat = abs(tab_vec-hr);%/max(abs(tab_vec));
    %err_mat = abs(tab_vec-hr)./abs_tab_vec;
    err_mat(isnan(err_mat)) = 0;
    err_mat(isinf(err_mat)) = 0;
    err_mat = mlf.vec2mat(err_mat,size(tab));
    
    %%% Model
    g = {pc w c};
    %[max_err,max_idx] = max(err_mat,[],'all');

    % set errors to zero to avoid no interpolation points to be added
    zero_idx = cell(1,n);
    for i = 1:n
        if length(ip_col{i}) >= opt.max_itpl(i)
            zero_idx{i} = 1:size(tab,i);
        else
            zero_idx{i} = ip_col{i};
        end
    end
    err_mat_greedy = err_mat;
    err_mat_greedy(zero_idx{:}) = 0;
    % maximum error for information
    [max_err,~] = max(err_mat,[],'all');
    % maximum error for greedy
    [~,max_idx] = max(err_mat_greedy,[],'all');

    %%% Plot
    fprintf('#%i rel max error %d \n \t interpolation points [ ',jj,max_err/max_samples)
    fprintf('%g ', cellfun(@length,ip_col));
    if (max_err < max_err_best) || (jj == 1)
        max_err_best    = max_err;
        g_best          = g;
        jj_best         = jj;
        for iii = 1:numel(pc)
            orders(iii) = length(pc{iii})-1;
        end
        fprintf(']*\n');
        no_improvement = 0;
    else
        fprintf(']\n');
        no_improvement = no_improvement + 1;
    end
    info.error(jj) = max_err;
    % if no_improvement > floor(opt.max_iter/5)
    %     CONTINUE = false;
    % end
end

%%% Output
r               = @(x) mlf.eval_lagrangian(g_best{1},g_best{2},g_best{3},x,false);
info.opt        = jj_best;
info.ord_est    = ip_est-1;
%
info.flop       = flop;
info.method     = method;
info.pr         = pr;
info.pc         = g_best{1};
info.w          = g_best{2};
info.c          = g_best{3};
info.ord        = orders;

