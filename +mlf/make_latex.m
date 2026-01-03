function latexList = make_latex(CAS,infoCas,tab,ip_data,opt)

FUN_SYM = @(x) vpa(x,3);
TAG     = false;
if nargin > 4
    TAG             = true;
    opt.ord_show    = false;
    if isfield(opt,'fun_sym')
        FUN_SYM = opt.fun_sym;
    end
end

%%% Generic stuff;
p_c     = ip_data.pc;
p_r     = ip_data.pr;
c       = ip_data.c;
w       = ip_data.w;
n       = numel(p_c);
%
vars    = '(';
vars_z  = '(';
vars_y  = '(';
for ii = 1:n
    vars    = [vars   's_{' num2str(ii) '},'];
    vars_z  = [vars_z 'z_{' num2str(ii) '},'];
    vars_y  = [vars_y 'y_{' num2str(ii) '},'];
    k_i(ii) = length(p_c{ii});
    q_i(ii) = length(p_r{ii});
end
N_i     = size(tab);
vars    = [vars(1:end-1)   ')'];
vars_z  = [vars_z(1:end-1) ')'];
vars_y  = [vars_y(1:end-1) ')'];
% Informations
infoData = [];
for ii = 1:numel(infoCas.tag)
    infoData = [infoData infoCas.tag{ii} ', '];
end
infoData = infoData(1:end-2);

% Intermediate 1-variables 
% vars_1  = 1; 
% ll      = 1;
% for ii = 1:numel(p_c)-1
%     ll      = ll*length(p_c{ii});
%     vars_1  = [vars_1 (ii+1)*ones(1,ll)];
% end

%%% Start
latexList   = [];

%%% Original function
latexList   = [latexList ['\noindent \textbf{Original function} ($l=1,\cdots,\ord=' num2str(n) '$) $$\bH' vars '=' infoCas.name(2:end-1) '$$ \\' ]];
latexList   = [latexList ['\noindent \textbf{Informations:} ' infoCas.ref ', ' infoCas.cite ' (' infoData ')\\~~\\' ]];
latexList   = [latexList ['\noindent \textbf{Original tensor} $$\tableau{' num2str(n) '}=' latex(sym(N_i)) '$$ \\' ]];

%%% Order detection
latexList   = [latexList ['\noindent \textbf{Order detection by single variable Loewner}: $$d_l=' latex(sym(ip_data.ord)) '$$' ]];
latexList   = [latexList ['\begin{figure}[H] \centering \includegraphics[width=\textwidth]{figures/svd_' num2str(CAS)  '} \caption{Normalized singular values in $' vars '$.} \end{figure}']];

%%% Interpolation points
latexList   = [latexList ['\noindent \textbf{Interpolation points}: $k_l=' latex(sym(k_i)) '$' ]];
if max(q_i) < 10
    latexList   = [latexList '$$ \begin{array}{rcl}'];
    for ii = 1:length(p_c)
        latexList = [latexList ['\lan{' num2str(ii) '} &=& ' latex(FUN_SYM(sym(p_c{ii}))) ' \\' ]];
        latexList = [latexList ['\mun{' num2str(ii) '} &=& ' latex(FUN_SYM(sym(p_r{ii}))) '\\' ]];
    end
    latexList   = [latexList '\end{array} $$' ];
else
    latexList   = [latexList '$$ \begin{array}{rclcrcl}'];
    for ii = 1:length(p_c)
        latexList = [latexList ['\lan{' num2str(ii) '} &=& ' latex(FUN_SYM(sym(p_c{ii}).')) ' &\text{ and }& ' ]];
        latexList = [latexList ['\mun{' num2str(ii) '} &=& ' latex(FUN_SYM(sym(p_r{ii}).')) '\\' ]];
    end
    latexList   = [latexList '\end{array} $$' ];
end

%%% Loewner c, w and c.w
latexList   = [latexList '\noindent \textbf{Lagrangian weights}:'];
top         = {'rc' 'rw' 'rcw'};
MAT         = sym([c w w.*c]);
for ii = 1:3
    tmp_data(:,ii) = [top{ii}; MAT(:,ii)];
end
latexList   = [latexList [ '$$' latex(FUN_SYM(tmp_data)) '$$' ]];

%%% Transfer function in Lagrangian
[glag,ilag] = mlf.tf_lagrangian(ip_data.pc,ip_data.w,ip_data.c,false);
[num,den]   = numden(glag);
latexList   = [latexList '\noindent \textbf{Lagrangian form} '];
latexList   = [latexList '(basis, numerator coefficients and denominator coefficients):\\~~\\'];
latexList   = [latexList '$$\textbf{Lag}=' latex(FUN_SYM([ilag.basis ilag.num_coeff ilag.den_coeff])) '$$ \\~~\\' ];
%
latexList   = [latexList '\noindent \textbf{Corresponding numerator and denominator}: \\~~\\'];
latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{lag}}' vars ' &=& \dfrac{\bn_{\textrm{lag}}' vars '}{\bd_{\textrm{lag}}' vars '}\\ &=& \dfrac{\sum \textbf{Lag}(:,1)^{-1}\odot \textbf{Lag}(:,2)}{\sum \textbf{Lag}(:,1)^{-1}\odot \textbf{Lag}(:,3)} \end{array}$$'];
latexList   = [latexList '$\bn_{\textrm{lag}}' vars ' = ' latex(FUN_SYM(num)) '$ \\~~\\'];
latexList   = [latexList '$\bd_{\textrm{lag}}' vars ' = ' latex(FUN_SYM(den)) '$ \\~~\\'];
%
% latexList   = [latexList ['\noindent \textbf{Equivalent simplified function $\bG_{\textrm{lag}}' vars '$}: \\~~\\ ']];
% latexList   = [latexList '$\bG_{\textrm{lag}}' vars ' = ' latex(FUN_SYM(simplify(num/den))) '$ \\~~\\'];

%%% Equivalent NN
if length(ip_data.c) < 20
    %
    latexList   = [latexList '\noindent \textbf{Connection with Neural Networks} (equivalent numerator $\bn_{\textrm{lag}}$ representation):\\ '];
    latexNN     = mlf.make_latex_NN_lag(p_c,double(w.*c)); 
    latexList   = [latexList ['\begin{figure}[H]\begin{center} \scalebox{.7}{' latexNN '} \caption{Equivalent NN representation of the numerator $\bn_{\textrm{lag}}$.}\end{center}\end{figure}']];
    %
    latexList   = [latexList '\noindent \textbf{Connection with Neural Networks} (equivalent denominator $\bd_{\textrm{lag}}$ representation):\\ '];
    latexNN     = mlf.make_latex_NN_lag(p_c,double(c)); 
    latexList   = [latexList ['\begin{figure}[H]\begin{center} \scalebox{.7}{' latexNN '} \caption{Equivalent NN representation of the denominator $\bd_{\textrm{lag}}$.}\end{center}\end{figure}']];
end

%%% Transfer function in Monomial
[gmon,imon] = mlf.tf_monomial(ip_data.pc,ip_data.w,ip_data.c,false);
%[num,den]   = numden(gmon);
num         = imon.num_coeff;
den         = imon.den_coeff;
% norm_max    = max(abs(den));
% num         = num/norm_max;
% den         = den/norm_max;
% den(abs(den)<1e-12) = 0;
% num(abs(num)<1e-12) = 0;
latexList   = [latexList '\noindent \textbf{Monomial form}'];
latexList   = [latexList ' (basis, numerator coefficients and denominator coefficients - small entries removed):\\~~\\'];
%latexList   = [latexList '$$' latex(fun([imon.basis imon.num_coeff imon.den_coeff ])) '$$ \\~~\\' ];
latexList   = [latexList '$$\textbf{Mon}=' latex(FUN_SYM([imon.basis num den])) '$$ \\~~\\' ];
latexList   = [latexList '\noindent \textbf{Corresponding numerator and denominator}:\\~~\\'];
num         = sum(num.*imon.basis);
den         = sum(den.*imon.basis);
latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{mon}}' vars ' &=& \dfrac{\bn_{\textrm{mon}}' vars '}{\bd_{\textrm{mon}}' vars '}\\ &=& \dfrac{\sum \textbf{Mon}(:,1)^{-1}\odot \textbf{Mon}(:,2)}{\sum \textbf{Mon}(:,1)^{-1}\odot \textbf{Mon}(:,3)} \end{array}$$'];
latexList   = [latexList '$\bn_{\textrm{mon}}' vars ' = ' latex(FUN_SYM(num)) '$ \\~~\\'];
latexList   = [latexList '$\bd_{\textrm{mon}}' vars ' = ' latex(FUN_SYM(den)) '$ \\~~\\'];
% latexList   = [latexList '\noindent \textbf{Equivalent simplified function $\bG_{\textrm{mon}}' vars '$}: \\~~\\ '];
% latexList   = [latexList '$\bG_{\textrm{mon}}' vars ' = ' latex(FUN_SYM(simplify(num/den))) '$ \\~~\\'];

%%% KST decoupling
[Bary,Lag,Cx]  = mlf.decoupling(ip_data.pc,ip_data.lag);
eval(['maxLength = length(Cx.d' num2str(n) ');'])
if maxLength > 10
    latexList = [latexList '\begin{landscape} '];
end
latexList = [latexList '\noindent \textbf{KST equivalent decoupling pattern}'];
% si's
latexList = [latexList '~\\ Barycentric weights: $$\begin{array}{rcll}'];
str_k     = [];
tmp_k     = [];
for ii = numel(p_c):-1:1
    eval(['tmp = Cx.s' num2str(ii) ';']);
    if ii == numel(p_c)
        latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \\' ]];
    else
        latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k '\\' ]];
    end
    tmp_k       = [tmp_k ['k_{' num2str(ii) '}']];
    str_k       = ['\bone_{' tmp_k '}'];
end
latexList = [latexList ['\end{array}$$' ]];
% w's
latexList = [latexList ['~\\ Data weights: $$\begin{array}{rcll}']];
str_k     = [];
tmp_k     = [];
for ii = numel(p_c):-1:1
    eval(['tmp = Cx.w' num2str(ii) ';']);
    if ii == numel(p_c)
        latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \\' ]];
    else
        latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k '\\' ]];
    end
    tmp_k       = [tmp_k ['k_{' num2str(ii) '}']];
    str_k       = ['\bone_{' tmp_k '}'];
end
latexList = [latexList '\end{array}$$' ];
% scalings
latexList = [latexList '~\\ Normalization scalings: $$\begin{array}{rcll}'];
str_k     = [];
tmp_k     = [];
for ii = numel(p_c):-1:1
    eval(['tmp = Cx.d' num2str(ii) ';']);
    if ii == numel(p_c)
        latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \\' ]];
    else
        latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k '\\' ]];
    end
    tmp_k       = [tmp_k ['k_{' num2str(ii) '}']];
    str_k       = ['\bone_{' tmp_k '}'];
end
latexList = [latexList '\end{array}$$' ];
if maxLength > 10
    latexList = [latexList '\end{landscape} '];
end

% KST den in matrix form
latexList = [latexList '\noindent Barycentric univariate vector-wise functions:'];
for ii = 1:length(p_c)
    tmp_den(:,ii)   = [['bary' num2str(ii)]; simplify((Lag{ii}.*Bary{ii}))];
end
latexList   = [latexList [ '$$\textbf{Bary}=' latex(FUN_SYM(tmp_den)) '$$' ]];
tmp_den     = tmp_den(2:end,:);
% Equivalent num & den
Den = 1;
for ii = 1:n
    Den = Den.*Lag{ii}.*Bary{ii}; % den
end
[tmp_num,~] = numden(simplifyFraction(sum(w.*Den)));
[tmp_den,~] = numden(simplifyFraction(sum(Den)));
%
latexList   = [latexList '\noindent \textbf{Corresponding numerator and denominator} (without support points): \\~\\ '];
latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{kst}}' vars ' &=& \dfrac{\bn_{\textrm{kst}}' vars '}{\bd_{\textrm{kst}}' vars '}\\ &=& \dfrac{\sum_i \bw\cdot \prod_j \textbf{Bary}_{i,j}}{\sum_i \prod_j \textbf{Bary}_{i,j}} \end{array}$$'];
latexList   = [latexList ['$\bn_{\textrm{kst}}' vars ' = ' latex(FUN_SYM(tmp_den)) '$ \\~\\']];
latexList   = [latexList ['$\bd_{\textrm{kst}}' vars ' =' latex(FUN_SYM(tmp_num)) '$ \\']];
% % KST num./den
% latexList   = [latexList '~\\ \noindent \textbf{Equivalent simplified function $\bG_{\textrm{kst}}' vars '$}: \\~\\ '];
% latexList   = [latexList ['$\bG_{\textrm{kst}}' vars '=' latex(FUN_SYM(simplify(tmp_num./tmp_den))) '$ \\']];

%%% Eval figure plot
latexList   = [latexList ['\begin{figure}[H] \centering \includegraphics[width=\textwidth]{figures/eval_' num2str(CAS)  '} \caption{Evaluation of $\bG_{\textrm{lag}}$ / $\bG_{\textrm{mon}}$ / $\bG_{\textrm{kst}}$ vs. original function $\bH$.} \end{figure}']];

%%% KST univariate functions 
[gkst1,kst] = mlf.make_singleVarFun(ip_data.pc,ip_data.pr,Cx);
gkst1       = kst.f;
latexList   = [latexList '~\\ \noindent \textbf{KST-like univariate functions} \\'];
latexList   = [latexList '~\\ \noindent Equivalent ' num2str(length(gkst1)) ' scaled univariate functions $\bphi_{1,\cdots,' num2str(length(gkst1)) '}$: \\ '];
latexListF  = '$$\left\{\begin{array}{rcrcl}';
idx         = [];
%inc         = 1;
for ii = 1:length(gkst1)
    eq_ii       = latex(FUN_SYM(gkst1{ii}));
    % ZVAR        = ['& &'];
    % if vars_1(ii) == inc
    %     ZVAR    = ['z_{' num2str(vars_1(ii)) '}&=&'];
    %     inc     = inc + 1;
    % end
    ZVAR  = ['z_{' num2str(ii) '} &=&'];
    if length(eq_ii) < 100
        %latexListF  = [latexListF [ZVAR '\bphi_{' num2str(ii) '}(s_{' num2str(vars_1(ii)) '}) &=& ' eq_ii '\\' ]];
        latexListF  = [latexListF [ZVAR '\bphi_{' num2str(ii) '}(s_{' num2str(ii) '}) &=& ' eq_ii '\\' ]];
    else
        idx         = [idx ii];
        [Num,Den]   = numden(gkst1{ii});
        nn{ii}      = latex(FUN_SYM(Num));
        dd{ii}      = latex(FUN_SYM(Den));
        %latexListF  = [latexListF [ZVAR '\bphi_{' num2str(ii) '}(s_{' num2str(vars_1(ii)) '}) &=& \frac{\bn_{' num2str(ii) '}}{\bd_{' num2str(ii) '}} \\' ]];
        latexListF  = [latexListF [ZVAR '\bphi_{' num2str(ii) '}(s_{' num2str(ii) '}) &=& \frac{\bn_{' num2str(ii) '}}{\bd_{' num2str(ii) '}} \\' ]];
    end
end
latexListF  = [latexListF '\end{array} \right.$$' ];
latexListF  = strrep(latexListF,'\frac','\cfrac');
latexList   = [latexList latexListF];
if ~isempty(idx)
    latexList   = [latexList '\noindent where, \\ '];
    for ii = 1:length(idx)
        latexList   = [latexList '\noindent $\bn_{' num2str(idx(ii)) '}=' nn{idx(ii)} '$ and \\ '];
        latexList   = [latexList '\noindent $\bd_{' num2str(idx(ii)) '}=' dd{idx(ii)} '$, \\ '];
    end
end

%%% Lifted interpolation points
if TAG
    [Var,Lag,Cx]    = mlf.decoupling(ip_data.pc,ip_data.lag);
    [gkst1,ikst]    = mlf.make_singleVarFun(ip_data.pc,ip_data.pr,Cx);
    ok              = mlf.check(ikst.pc,ikst.pr);
    [psi,iloe2]     = mlf.alg1(ip_data.tabr,ikst.pc,ikst.pr,opt);
    psi             = sym(simplify(mlf.tf_lagrangian(iloe2.pc,iloe2.w,iloe2.c,false)));
    % if isnan(psi)
    %     psi = sym(NaN);%'\textrm{NaN}';
    % end
    % KST like connecting function
    latexList   = [latexList ['~\\ \noindent \textbf{Connecting the scaled univariate function with rational $\bPsi' vars_y '$}: \\~\\ ']];
    for ii = 1:n
        psi = subs(psi,['s' num2str(ii)],['y' num2str(ii)]);
    end
    latexList   = [latexList ['$\bPsi' vars_y '= ' latex(FUN_SYM(psi)) ' $ \\']];
    latexList   = [latexList ['~\\ \noindent \textbf{Leading to}:\\~\\ ']];
    gkst        = psi;
    for ii = 1:n
        gkst = subs(gkst,['y' num2str(ii)],ikst.f{ii});
    end
    gkst        = simplify(gkst);
    vars_Phi    = '(';
    for ii = 1:n
        vars_Phi = [vars_Phi '\bphi(\var{' num2str(ii) '})'];
    end
    vars_Phi    = [vars_Phi ')'];
    latexList   = [latexList ['$\bPsi' vars_Phi '=\bG_{\textrm{kst}}' vars '=' latex(FUN_SYM(gkst)) '$']];
    %latexList   = [latexList ['$\bG_{\textrm{kst}}' vars '= \bPsi' vars_z '=' latex(FUN_SYM(gkst)) '$']];
end

%%% Some minor text-style
latexList = strrep(latexList,'\mathrm{rc}','\bc');
latexList = strrep(latexList,'\mathrm{rw}','\bw');
latexList = strrep(latexList,'\mathrm{rcw}','\bc\cdot\bw');
% Replace bary with \bc 
for ii = 1:numel(p_c)
    latexList = strrep(latexList,['\mathrm{bary}_{' num2str(ii) '}'],['\bc^{\var{' num2str(ii) '}} \mathbf{Lag}(\var{' num2str(ii) '})']);
end
% Replace s with \var 
for ii = 1:numel(p_c)
    latexList = strrep(latexList,['s_{' num2str(ii) '}'],['\var{' num2str(ii) '}']);
end
% \frac
% latexList = strrep(latexList,'\frac','\cfrac');