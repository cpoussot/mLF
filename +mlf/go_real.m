function [TL,TR] = go_real(p_c)

%%% Go real
J0  = sqrt(2)/2*[1 -1i; 1 1i];
n   = length(p_c);
n1r = sum(imag(p_c{1})==0);
n1c = (length(p_c{1})-n1r)/2;
for ii = 1:n
    if ii == 1
        Tc{ii} = eye(n1r);
        for jj = 1:n1c
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