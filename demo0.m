clearvars, close all, clc
set(groot,'DefaultFigurePosition', [300 300 1000 400]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
%addpath("location_of_mlf") % Add the location of the +mlf package

%%% Define a multivariate handle function 
n       = 3; % number of variables
%H       = @(s1,s2,s3) (s3/100-1)*(s2-pi/2)*(s1+atan(2*s2)*tanh(5*(s2-pi)))/(s1^2+s3/10*cos(3*s1)+3)/(s2+10);
H       = @(s1,s2,s3) s1*(s2/10)^2+s3;
%H       = @(s1,s2,s3) s2^3/(s1^2+2*s1-s3/2+10);
% /!\ The tolerence is an important parameter when the data are generated from an irrational function
tol_ord = 1e-7; 
% Interpolation points (IP) - separate columns and rows (as Section 3, eq. 13-15)
for ii = 1:n
    p_c{ii} = linspace(-6,6,30);
    dx      = abs(p_c{ii}(2)-p_c{ii}(1))/2;
    p_r{ii} = p_c{ii}+dx;
end

%%% Elementary ingredients
% Check that column/row IP are disjoint (Section 3)
ok                  = mlf.check(p_c,p_r);
% Construct tab_n (Section 3)
tab                 = mlf.make_tab(H,p_c,p_r,true);

%%% Estimate order along each variables and select a subset of IP
ord                 = mlf.compute_order(p_c,p_r,tab,tol_ord,[],5,true);
[pc,pr,W,V,tab_red] = mlf.points_selection(p_c,p_r,tab,ord,false);
w                   = mlf.mat2vec(W);

%%% Recursive null-space computation (Section 5, Theorem 5.3)
[c,info]            = mlf.loewner_null_rec(pc,pr,tab_red,'svd0');

%%% Curse of dimensionality: flops & memory estimation (Section 5, Theorem 5.5 & Theorem 5.6)
fprintf('FLOPS\n')
fprintf(' * recursive: %d\n',info.nflop)
fprintf(' * full     : %d\n',length(c)^3)
fprintf('MEMORY\n')
fprintf(' * recursive: %d MB\n',max(ord+1)^2/2^20)
fprintf(' * full     : %d MB\n',prod(ord+1)^2/2^20)

%%% 3D plot 
% Along first and second variables 
% Other variables are randomly chosen between bounds
x1      = linspace(min(p_r{1}),max(p_r{1}),40)+rand(1)/10;
x2      = linspace(min(p_r{1}),max(p_r{1}),41)+rand(1)/10;
[X,Y]   = meshgrid(x1,x2);
rnd_p   = [];
if n > 2; rnd_p = mlf.rand(n-2,p_r(3:end),false); end
for ii = 1:numel(x1)
    for jj = 1:numel(x2)
        param                   = [x1(ii) x2(jj) rnd_p];
        paramStr                = regexprep(num2str(param,36),'\s*',',');
        eval(['tab_ref(jj,ii)   = H(' paramStr ');'] )
        tab_app(jj,ii)          = mlf.eval_lagrangian(pc,w,c,param,false);
    end
end
%
figure
subplot(1,2,1); hold on, grid on
surf(X,Y,tab_app,'EdgeColor','none'), hold on
surf(X,Y,tab_ref,'EdgeColor','k','FaceColor','none')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('Original vs. Approximation','Interpreter','latex')
axis tight, zlim([min(tab_ref(:)) max(tab_ref(:))]), view(-30,40)
subplot(1,2,2); hold on, grid on, axis tight
imagesc(log10(abs(tab_ref-tab_app)/max(abs(tab_ref(:)))),'XData',x1,'YData',x2)
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('{\bf log}(abs. err.)/max.','Interpreter','latex')
colorbar,

%%% Decoupling features for KST
[Var,Lag,Bary]      = mlf.decoupling(pc,info);
%[Lag{1} Lag{2} Lag{3}]