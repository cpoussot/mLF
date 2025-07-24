function [c,info_rec]  = loewner_null_rec_sym(p_c,p_r,tab,TAG)

n = length(p_c);
%
for ii = 1:numel(p_c)
    Nc(ii) = length(p_c{ii});
end
% Optimize flop
if nargin > 3
    if isempty(TAG)
        %
    elseif strcmp(TAG(1:5),'optim') && n > 1
        % Reorder 
        for ii = 1:length(p_c)
            ord(ii) = length(p_c{ii});
        end
        [ord_sort,idx] = sort(ord,'descend');
        % Update tab variables and p_c / p_r
        tab = permute(tab,idx);
        for ii = 1:numel(idx)
            tmp1{ii}    = p_c{idx(ii)};
            tmp2{ii}    = p_r{idx(ii)};
        end
        p_c = tmp1;
        p_r = tmp2;
    end
end
%
switch n
    case 1
        warning('No need of recursive computation')
        W               = tab(1:length(p_c{1}));
        V               = tab(1+length(p_c{1}):end);
        LL              = mlf.loewnerMatrix_sym(p_c,p_r,W,V);
        c               = null(LL);
        info_rec.LL     = LL;
        info_rec.nflop  = size(LL,1)*size(LL,2)^2;
        info_rec.c1     = c;
        info_rec.w1     = W(:);
        info_rec.d1     = 1;
    case 2
        [c,info_rec]    = mlf.loewner_null_rec2_sym(p_c,p_r,tab);
    otherwise
        [c,info_rec]    = mlf.loewner_null_recN_sym(p_c,p_r,tab,Nc);
end
% Optimize flop (put back in the right order / or not)
if nargin > 3
    if strcmp(TAG,'optim')
        combi1  = mlf.combinations_dim(ord);
        combi2  = mlf.combinations_dim(ord_sort);
        combi1  = combi1(:,idx);
        for ii = 1:size(combi2,1)
            jj          = 0;
            CONTINUE    = true;
            while CONTINUE
                jj = jj + 1;
                if combi1(jj,:) == combi2(ii,:)
                    idx_new(ii) = jj;
                    CONTINUE    = false;
                end
            end
        end
        for ii = 1:numel(idx_new)
            cnew(idx_new(ii),1) = c(ii);
        end
        c = cnew;
    elseif strcmp(TAG,'optim_fast') && n > 1
        info_rec.p_c    = p_c;
        info_rec.p_r    = p_r;
        info_rec.tab    = tab;%permute(tab,idx);
        str             = [];
        for ii = 1:numel(ord_sort)
            str = [str '1:' num2str(ord_sort(ii)) ','];
        end
        eval(['info_rec.W      = tab(' str(1:end-1) ');']);
        info_rec.idx    = idx;
    end
end