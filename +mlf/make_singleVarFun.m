function [g,info] = make_singleVarFun(p_c,p_r,Bary)

SAVEIT = false;
jj_obj = 1;
info = [];
kk = 0;
n  = numel(p_c);
for ii = 1:n
    eval(['cc = Bary.s' num2str(ii) ';']);
    eval(['ww = Bary.w' num2str(ii) ';']);
    eval(['dd = Bary.d' num2str(ii) ';']);
    for jj = 1:size(cc,2)
        kk          = kk + 1;
        tmp         = mlf.tf_lagrangian({p_c{ii}},ww(:,jj),cc(:,jj),false);
        tmp         = subs(tmp,'s1',['s' num2str(ii)]);
        %tmp         = subs(tmp,'s1',['s' num2str(n+1-ii)]);
        g{kk,1}     = simplify(tmp);
        gs{kk,1}    = 1/dd(jj) * g{kk,1};
        % if (jj == jj_obj) && (g{kk,1} ~= 0)
        %     SAVEIT = true;
        % elseif (g{kk,1} ~= 0)
        %     jj_obj = jj_obj + 1;
        % end
        if (jj == 1) %SAVEIT 
            info.f{ii,1}    = g{kk,1};
            fun_tmp         = matlabFunction(info.f{ii});
            info.pc{ii}     = fun_tmp(p_c{ii});
            info.pr{ii}     = fun_tmp(p_r{ii});
            %SAVEIT          = false;
            %jj_obj          = 1;
        end
    end
end
info.gs = gs;