function latexList = make_latex(CAS,infoCas,ip_data,opt)

FUN_SYM = @(x) vpa(x,3);
TAG     = false;
if nargin > 4
    TAG             = true;
    opt.ord_show    = false;
    if isfield(opt,'fun_sym')
        FUN_SYM = opt.fun_sym;
    end
end
Nmax = 75;

%%% Generic stuff;
p_c     = ip_data.pc;
p_r     = ip_data.pr;
c       = ip_data.c;
w       = ip_data.w;
lag     = ip_data.lag;
n       = numel(p_c);
N       = numel(c);
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
for ii = 1:n
    N_i(ii) = length(infoCas.ip{ii});
end
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
latexList   = [latexList ['\begin{figure}[H] \centering \includegraphics[width=\textwidth]{case_' num2str(CAS) '/svd} \caption{Normalized singular values in $' vars '$.} \end{figure}']];

if N < Nmax
    %%% Interpolation points
    latexList   = [latexList '\noindent \textbf{Interpolation points} ($k_l=' latex(sym(k_i)) '$, $q_l=' latex(sym(q_i)) '$):'];
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
    [glag,ilag] = mlf.tf_lagrangian(p_c,w,c,false);
    [num,den]   = numden(glag);
    latexList   = [latexList '\noindent \textbf{Lagrangian form} '];
    latexList   = [latexList '(basis, numerator and denominator coefficients):'];
    latexList   = [latexList '$$\left(\begin{array}{ccc}\mathcal{B}_\textrm{lag}' vars ' & \bN_\textrm{lag} &\bD_\textrm{lag}\end{array}\right) =$$ $$' latex(FUN_SYM([ilag.basis ilag.num_coeff ilag.den_coeff])) '.$$' ];
    latexList   = [latexList '\noindent The corresponding function is:'];
    latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{lag}}' vars ' &=& \dfrac{\bn_{\textrm{lag}}' vars '}{\bd_{\textrm{lag}}' vars '}\\ && \\ &=& \dfrac{\sum_{\textrm{row}} \bN_\textrm{lag} \odot\mathcal{B}^{-1}_\textrm{lag}' vars '}{\sum_{\textrm{row}} \bD_\textrm{lag} \odot\mathcal{B}^{-1}_\textrm{lag}' vars '}, \end{array}$$'];
    latexList   = [latexList '\noindent where,\\'];
    latexList   = [latexList '$\bn_{\textrm{lag}}' vars ' = ' latex(FUN_SYM(num)) '$ \\~~\\'];
    latexList   = [latexList '$\bd_{\textrm{lag}}' vars ' = ' latex(FUN_SYM(den)) '$ \\~~\\'];
    
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
    [gmon,imon] = mlf.tf_monomial(p_c,w,c,false);
    num         = imon.num_coeff;
    den         = imon.den_coeff;
    latexList   = [latexList '\noindent \textbf{Monomial form}'];
    latexList   = [latexList ' (basis, numerator and denominator coefficients - entries $<10^{-12}$ removed):'];
    latexList   = [latexList '$$\left(\begin{array}{ccc}\mathcal{B}_\textrm{mon}' vars ' & \bN_\textrm{mon} &\bD_\textrm{mon}\end{array}\right) =$$ $$' latex(FUN_SYM([imon.basis imon.num_coeff imon.den_coeff])) '$$' ];
    latexList   = [latexList '\noindent The corresponding function is:'];
    num         = sum(num.*imon.basis);
    den         = sum(den.*imon.basis);
    latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{mon}}' vars ' &=& \dfrac{\bn_{\textrm{mon}}' vars '}{\bd_{\textrm{mon}}' vars '}\\ && \\&=& \dfrac{\sum_{\textrm{row}} \bN_\textrm{mon} \odot \mathcal{B}_\textrm{mon}' vars '}{\sum_{\textrm{row}} \bD_\textrm{mon} \odot\mathcal{B}_\textrm{mon}' vars '},  \end{array}$$'];
    latexList   = [latexList '\noindent where,\\'];
    latexList   = [latexList '$\bn_{\textrm{mon}}' vars ' = ' latex(FUN_SYM(num)) '$ \\~~\\'];
    latexList   = [latexList '$\bd_{\textrm{mon}}' vars ' = ' latex(FUN_SYM(den)) '$ \\~~\\'];
    
    %%% KST decoupling
    [Bary,Lag,Cx]  = mlf.decoupling(p_c,lag);
    eval(['maxLength = length(Cx.d' num2str(n) ');'])
    if maxLength > 10
        latexList = [latexList '\begin{landscape} '];
    end
    latexList = [latexList '\noindent \textbf{KST equivalent decoupling pattern}'];
    % si's
    latexList = [latexList ' (Barycentric weights $\bc^{\var{l}}$): $$\begin{array}{rclll}'];
    str_k     = [];
    tmp_k     = [];
    for ii = numel(p_c):-1:1
        eval(['tmp = Cx.s' num2str(ii) ';']);
        if ii == numel(p_c)
            latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) & := \textbf{Bary}(\var{' num2str(ii) '}) \\' ]];
        else
            latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k ' & := \textbf{Bary}(\var{' num2str(ii) '}) \\' ]];
        end
        tmp_k   = [tmp_k ['k_{' num2str(ii) '}']];
        str_k   = ['\bone_{' tmp_k '}'];
    end
    latexList = [latexList ['\end{array}$$' ]];
    % % w's
    % latexList = [latexList ['~\\ Data weights $\bw^{\var{l}}$: $$\begin{array}{rcll}']];
    % str_k     = [];
    % tmp_k     = [];
    % for ii = numel(p_c):-1:1
    %     eval(['tmp = Cx.w' num2str(ii) ';']);
    %     if ii == numel(p_c)
    %         latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \\' ]];
    %     else
    %         latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k '\\' ]];
    %     end
    %     tmp_k   = [tmp_k ['k_{' num2str(ii) '}']];
    %     str_k   = ['\bone_{' tmp_k '}'];
    % end
    % latexList = [latexList '\end{array}$$' ];
    % % scalings
    % latexList = [latexList '~\\ Normalization scalings $d^{\var{l}}$: $$\begin{array}{rcll}'];
    % str_k     = [];
    % tmp_k     = [];
    % for ii = numel(p_c):-1:1
    %     eval(['tmp = Cx.d' num2str(ii) ';']);
    %     if ii == numel(p_c)
    %         latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \\' ]];
    %     else
    %         latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k '\\' ]];
    %     end
    %     tmp_k   = [tmp_k ['k_{' num2str(ii) '}']];
    %     str_k   = ['\bone_{' tmp_k '}'];
    % end
    % latexList = [latexList '\end{array}$$' ];
    if maxLength > 10
        latexList = [latexList '\end{landscape} '];
    end
    % 
    latexList = [latexList '~\\ Then, with the above notations, one defines the following univariate vector functions:  $$ \left\{ \begin{array}{rcl}'];
    for ii = 1:n
        latexList = [latexList ['\bPhi_{' num2str(ii) '}(\var{' num2str(ii) '}) &:=& \textbf{Bary}(\var{' num2str(ii) '}) \odot \mathcal{B}_{\textrm{lag}}^{-1}(\var{' num2str(ii) '}) \\']];
    end
    latexList = [latexList '\end{array} \right. $$'];

    %%% Resulting KST
    latexList   = [latexList '\noindent The corresponding function is:'];
    latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{kst}}' vars ' &=& \dfrac{\bn_{\textrm{kst}}' vars '}{\bd_{\textrm{kst}}' vars '}\\ &=& \dfrac{\sum_{\text{rows}} \bw \odot \bPhi_{1}(\var{1}) \odot \cdots \odot\bPhi_{' num2str(n) '}(\var{' num2str(n) '})}{\sum_{\text{rows}} \bPhi_{1}(\var{1}) \odot \cdots \odot\bPhi_{' num2str(n) '}(\var{' num2str(n) '})} . \end{array}$$'];

    %%% Eval figure plot
    latexList   = [latexList ['\begin{figure}[H] \centering \includegraphics[width=\textwidth]{case_' num2str(CAS) '/eval} \caption{Evaluation of $\bG_{\textrm{lag}}$ / $\bG_{\textrm{mon}}$ / $\bG_{\textrm{kst}}$ vs. original function $\bH$.} \end{figure}']];

    %%% KST univariate functions 
    [gkst1,kst] = mlf.make_singleVarFun(p_c,p_r,Cx);
    gkst1       = kst.f;
    latexList   = [latexList '~\\ \noindent \textbf{KST-like univariate functions} '];
    latexList   = [latexList '(equivalent scaled univariate functions $\bphi_{1,\cdots,' num2str(length(gkst1)) '}$): '];
    latexListF  = '$$\left\{\begin{array}{rcrcl}';
    idx         = [];
    %inc         = 1;
    for ii = 1:length(gkst1)
        eq_ii   = latex(FUN_SYM(gkst1{ii}));
        ZVAR    = ['z_{' num2str(ii) '} &=&'];
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
    latexListF  = [latexListF '\end{array} \right. .$$' ];
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
else
    %%% Eval figure plot
    latexList   = [latexList ['\begin{figure}[H] \centering \includegraphics[width=\textwidth]{case_' num2str(CAS) '/eval} \caption{Evaluation of $\bG_{\textrm{lag}}$ / $\bG_{\textrm{mon}}$ / $\bG_{\textrm{kst}}$ vs. original function $\bH$.} \end{figure}']];
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