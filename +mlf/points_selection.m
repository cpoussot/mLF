function [P_c,P_r,W,V,tab] = points_selection(p_c,p_r,tab,ord,minimal,CONJUGATE)

if nargin < 6
    CONJUGATE = false;
end

n       = length(p_c);
ord     = ord+1;
tab_str = [];
W_str   = [];
V_str   = [];
for ii = 1:n
    ni(ii)      = length(p_c{ii});
    if CONJUGATE && (ii == 1)
        idxp_c{1}   = unique(floor(linspace(1,numel(p_c{1})-1,floor(ord(1))/2)));
        idxp_c{1}   = floor(idxp_c{1}/2)*2+1;
        idxp_c{1}   = sort([idxp_c{1} idxp_c{1}+1]);
    else
        idxp_c{ii}  = unique(floor(linspace(1,numel(p_c{ii}),ord(ii))));
    end
    %idxp_c{ii}  = unique(floor([1 linspace(3,numel(p_c{ii}),ord(ii)-1)]));
    jj          = 1;
    if ~CONJUGATE 
        while length(idxp_c{ii}) < ord(ii)
            idxp_c{ii}  = unique([idxp_c{ii} jj]);
            jj          = jj + 1;
        end
    end
    if minimal
        if CONJUGATE && (ii == 1)
            idxp_r{1}   = unique(floor(linspace(1,numel(p_r{1})-1,floor(ord(1))/2)));
            idxp_r{1}   = floor(idxp_r{1}/2)*2+1;
            idxp_r{1}   = sort([idxp_r{1} idxp_r{1}+1]);
        else
            idxp_r{ii}  = unique(floor(linspace(1,numel(p_r{ii}),ord(ii))));
        end
        jj = 1;
        if ~CONJUGATE 
            while length(idxp_r{ii}) < ord(ii)
                idxp_r{ii}  = unique([idxp_r{ii} jj]);
                jj          = jj + 1;
            end
        end
    else
        idxp_r{ii}  = 1:numel(p_r{ii});
    end
    p_c{ii} = p_c{ii}(idxp_c{ii});
    p_r{ii} = p_r{ii}(idxp_r{ii});
    %p_c{1}.'
    %p_r{1}.'
    %pause
    % 
    ii_str  = num2str(ii);
    tab_str = [tab_str '[idxp_c{' ii_str  '} ni(' ii_str ')+idxp_r{' ii_str '}],'];
    W_str   = [W_str 'idxp_c{' ii_str  '},'];
    V_str   = [V_str 'ni(' ii_str ')+idxp_r{' ii_str '},'];
end
eval(['W=tab(' W_str(1:end-1) ');']);
eval(['V=tab(' V_str(1:end-1) ');']);
eval(['tab=tab(' tab_str(1:end-1) ');']);
P_c     = p_c;
P_r     = p_r;
info    = [];

% % Tableau cell
% n1      = length(p_c{1});
% n2      = length(p_c{2});
% tabCell = {tab(1:n1,1:n2),     tab(1:n1,n2+1:end); ...
%            tab(1+n1:end,1:n2), tab(1+n1:end,n2+1:end)};
% % Get data to interpolate
% P_c     = [ptCell1(1), ptCell2(1), tabCell(1,1)];
% P_r     = [ptCell1(2), ptCell2(2), tabCell(2,2)];
% 
% 
% Check if first line is not constant
% ok  = true;
% if norm(abs(max(tab(1,:))-min(tab(1,:)))) < 1e-10
%     warning('''tableau(1,:)'' block not appropriate for CoD null-space computation')
%     ok = false;
% end
% info.tab_row1   = ok;
