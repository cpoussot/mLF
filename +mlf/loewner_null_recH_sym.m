function [c,info_rec]  = loewner_null_recH_sym(p_c,p_r,H,TAG)

n = length(p_c);
for ii = 1:numel(p_c)
    Nc(ii) = length(p_c{ii});
end
% Optimize flop
if nargin > 3
    if isempty(TAG)
        %
    elseif strcmp(TAG,'optim')
        % Reorder 
        for ii = 1:length(p_c)
            ord(ii) = length(p_c{ii});
        end
        [ord_sort,idx] = sort(ord,'descend');
        % Update H variables and p_c / p_r
        for ii = 1:numel(idx)
            H           = subs(H,['s' num2str(idx(ii))],['t' num2str(ii)]);
            tmp1{ii}    = p_c{idx(ii)};
            tmp2{ii}    = p_r{idx(ii)};
        end
        for ii = 1:numel(ord_sort)
            H   = subs(H,['t' num2str(ii)],['s' num2str(ii)]);
        end
    end
    p_c = tmp1;
    p_r = tmp2;
end

switch n
    case 1
        warning('No need of recursive computation')
        W               = mlf.eval_sym(H,p_c,false);
        V               = mlf.eval_sym(H,p_r,false);
        LL              = mlf.loewnerMatrix_sym(p_c,p_r,W,V);
        c               = null(LL);
        info_rec.LL     = LL;
        info_rec.nflop  = size(LL,1)*size(LL,2)^2;
        info_rec.c1     = sym(c);
    case 2
        [c,info_rec]    = mlf.loewner_null_rec2H_sym(p_c,p_r,H);
    otherwise
        [c,info_rec]    = mlf.loewner_null_recNH_sym(p_c,p_r,H,Nc);        
end
% Optimize flop (put back in the right order)
if nargin > 3
    if strcmp(TAG,'optim')
        info_rec.p_c    = p_c;
        info_rec.p_r    = p_r;
        %
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
            % for jj = 1:size(combi1,1)
            %     if combi1(jj,:) == combi2(ii,:)
            %         idx_new(ii) = jj;
            %     end
            % end
        end
        info_rec.idx    = idx_new;
        for ii = 1:numel(idx_new)
            cnew(idx_new(ii),1) = c(ii);
        end
        c = cnew;
    end
end