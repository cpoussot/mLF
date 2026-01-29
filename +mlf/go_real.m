function [TL,TR] = go_real(p_c)

%%% Go real
J0  = sqrt(2)/2*[1 -1i; 1 1i];
n   = length(p_c);
n1  = length(p_c{1})/2;
for ii = 1:n
    if ii == 1
        Tc{ii} = [];
        for jj = 1:n1
            Tc{ii} = blkdiag(Tc{ii},J0);
        end
    end
    T{ii} = eye(length(p_c{ii}));    
end
for ii = 1:n
    if ii == 1
        TT{ii} = Tc{1};
        for jj = 2:n
            TT{ii} = kron(TT{ii},T{jj});
        end
    else
        TT{ii} = 1;
        for jj = 1:n
            TT{ii} = kron(TT{ii},T{jj});
        end
    end
end
TL = 1; TR = 1;
for ii = 1:n
    TL = TL'*TT{ii};
    TR = TR*TT{ii};
end