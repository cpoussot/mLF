% Description
% Allow to generate a combination matrix replacing a sequence of FOR loop
% with a single one. More specifically, 
% for i1 = 1:k1
%    for i2 = 1:k2
%       ....
%       for in = 1:kn
%          a = tab(i1,i2,...,in)
% ...
% May be replaced by
% 
% 
% 
% 
% Input
%  * tab : is a n-dimensional tableau with dimensions k1 x k2 x ... x kn
% 
% Output
%  * comb : is the combination matrix of tab's entries to simulate an
%           imbricated for loop. 
% 

function comb = combinations(tab)

dim     = size(tab);
n       = length(dim);
allDim  = cell(n,1);
for ii = 1:n
    allDim{n-ii+1}  = 1:dim(ii);
end
comb        = cell(1,n);
[comb{:}]   = ndgrid(allDim{:});
comb        = reshape(cat(n+1,comb{:}),[],n);
comb        = fliplr(comb);
% mean(comb(:,1)) 
% if mean(comb(:,1)) == 1
%     comb = comb(:,2);
% end

% %%% Reshape
% avg         = mean(comb,1);
% idx_keep    = [];
% for ii = 1:n
%     if avg(ii) ~= 1
%         idx_keep = [idx_keep ii];
%     end
% end
% idx_rm  = setdiff(1:n,idx_keep);
% comb    = comb(:,idx_keep);
% N       = size(comb,1);
% comb_   = [];
% kk      = 1;
% for ii = 1:n
%     %[any(ii == idx_keep) any(ii == idx_rm)]
%     if any(ii == idx_keep)
%         comb_   = [comb_ comb(:,kk)];
%         kk      = kk + 1;
%     elseif any(ii == idx_rm)
%         comb_   = [comb_ ones(N,1)];
%     end
% end
% comb = comb_;
% % rm first
% if mean(comb(:,1)) == 1
%     comb = comb(:,2);
%     %comb = fliplr(comb);
% end