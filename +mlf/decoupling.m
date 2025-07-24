function [Var,Lag,Bary] = decoupling(p_c,info_rec)

n = length(p_c);

%%% Symbolic
for ii = 1:n
    eval(['s' num2str(ii) ' = sym(''s' num2str(ii) ''');']);
end

%%% Barycentric coefficients
switch n
    case 1
        bary.s1  = info_rec.c1;
        baryW.w1 = info_rec.w1;
        baryD.d1 = info_rec.d1(:);
    case 2
        bary.s1  = info_rec.c1;
        baryW.w1 = info_rec.w1;
        baryD.d1 = info_rec.d1(:);
        bary.s2  = info_rec.c2;
        baryW.w2 = info_rec.w2;
        baryD.d2 = info_rec.d2(:);
    otherwise
        for ii = 1:n
            eval(['c.s' num2str(ii) '= [];']);
            eval(['cw.w' num2str(ii) '= [];']);
            eval(['cs.d' num2str(ii) '= [];']);
        end
        [bary,~]    = mlf.bary_data2cN(c,{info_rec});
        [baryW,~]   = mlf.bary_data2wN(cw,{info_rec});
        [baryD,~]   = mlf.bary_data2dN(cs,{info_rec});
end

% Force symbolic 
for ii = 1:n
    eval(['bary_.s' num2str(ii) '= sym(bary.s' num2str(ii) ');']);
    eval(['bary_.w' num2str(ii) '= sym(baryW.w' num2str(ii) ');']);
    eval(['bary_.d' num2str(ii) '= sym(baryD.d' num2str(ii) ');']);
end

%%% Kronecker groups
for ii = 1:length(p_c)
    kron_left{ii}   = 1;
    kron_right{ii}  = 1;
    for jj = 1:length(p_c)
        if ii > jj
            kron_left{ii}   = kron(kron_left{ii},ones(1,length(p_c{jj})));
        elseif ii < jj
            kron_right{ii}  = kron(kron_right{ii},ones(1,length(p_c{jj})));
        end
    end
    kron_left{ii}   = kron_left{ii}.';
    kron_right{ii}  = kron_right{ii}.';
end

%%% Lagrangian basis
for ii = 1:length(p_c)
    eval(['stmp = s' num2str(ii) ';'])
    Lag{ii} = 1./kron(kron(kron_left{ii},stmp-p_c{ii}(:)),kron_right{ii});
end

%%% Null variables toward KST
for ii = 1:length(p_c)
    eval(['tmp = bary_.s' num2str(n-ii+1) ';'])
    %eval(['tmpD = bary_.d' num2str(n-ii+1) ';'])
    %eval(['tmpW = bary_.w' num2str(n-ii+1) ';'])
    Var{ii}     = [];
    %Scal{ii}    = [];
    for jj = 1:size(tmp,2)
        Var{ii}     = [Var{ii};  kron(tmp(:,jj),kron_right{ii})];
        %Scal{ii}    = [Scal{ii};  kron(tmpD(jj),kron_right{ii})];
        %VarW{ii}    = [VarW{ii}; kron(tmpW(:,jj),kron_right{ii})];
    end
end

% Re-order bary
for ii = 1:n
    eval(['Bary.s' num2str(ii) '= sym(bary_.s' num2str(n+1-ii) ');']);
    eval(['Bary.w' num2str(ii) '= sym(bary_.w' num2str(n+1-ii) ');']);
    eval(['Bary.d' num2str(ii) '= sym(bary_.d' num2str(n+1-ii) ');']);
end
