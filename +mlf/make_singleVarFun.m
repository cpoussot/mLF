function [g,gs] = make_singleVarFun(p_c,Bary)

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
    end
end