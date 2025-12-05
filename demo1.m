clearvars, close all, clc
set(groot,'DefaultFigurePosition', [300 300 1000 400]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory')); index_interpreter = find(contains(list_factory,'Interpreter')); for i = 1:length(index_interpreter); set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex'); end
%addpath("location_of_mlf") % Add the location of the +mlf package
%addpath('/Users/charles/Mon Drive/Research/Codes/')

%%% Pick an example from 1 to 40
CAS         = 1% 33

%%% mlf parameters
alg1_tol    = 1e-9; 
alg2_tol    = 1e-9;

%%% Chose model
[H,infoCas] = mlf.examples(CAS)
n           = infoCas.n;
p_c         = infoCas.p_c;
p_r         = infoCas.p_r;

%%% Data tensor/rand
[y,x,dim]   = mlf.make_tab_vec(H,p_c,p_r);
tab         = mlf.vec2mat(y,dim);

%%% Alg. 1: direct pLoe [A/G/P-V, 2025]
opt = [];
tic;
opt.ord_tol     = alg1_tol;
opt.method_null = 'svd0';
opt.method      = 'rec';
% opt.ord_obj     = [];
% opt.ord_N       = 10;
% opt.ord_show    = false;
% opt.data_min    = true;
[r_loe1r,i1_r]  = mlf.alg1(tab,p_c,p_r,opt);
toc

%%% Alg. 2: iterative pLoe [A/G/P-V, 2025]
opt = [];
tic
opt.tol         = alg2_tol;
opt.method_null = 'svd0';
opt.max_iter    = 25;
opt.method      = 'rec';
[r_loe2r,i2_r]  = mlf.alg2(tab,p_c,p_r,opt);
toc

%%% 
% Along first and second variables 
% Other variables are randomly chosen between bounds
x1      = linspace(min(infoCas.bound{1}),max(infoCas.bound{1}),40)+rand(1)/10;
x2      = linspace(min(infoCas.bound{2}),max(infoCas.bound{2}),41)+rand(1)/10;
[X,Y]   = meshgrid(x1,x2);
rnd_p   = [];
if n > 2; rnd_p = mlf.rand(n-2,p_r(3:end),false); end
for ii = 1:numel(x1)
    for jj = 1:numel(x2)
        param           = [x1(ii) x2(jj) rnd_p];
        tab_ref(jj,ii)  = H(param);
        tab_app1(jj,ii) = r_loe1r(param);
        tab_app2(jj,ii) = r_loe2r(param);
    end
end
%
figure
subplot(2,2,1); hold on, grid on
surf(X,Y,tab_app1,'EdgeColor','none'), hold on
surf(X,Y,tab_ref,'EdgeColor','k','FaceColor','none')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('mLF alg. 1 (direct)','Interpreter','latex')
axis tight, zlim([min(tab_ref(:)) max(tab_ref(:))]), view(-30,40)
subplot(2,2,2); hold on, grid on, axis tight
imagesc(log10(abs(tab_ref-tab_app1)/max(abs(tab_ref(:)))),'XData',x1,'YData',x2)
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('{\bf log}(abs. err.)/max.','Interpreter','latex')
colorbar,
%
subplot(2,2,3); hold on, grid on
surf(X,Y,tab_app2,'EdgeColor','none'), hold on
surf(X,Y,tab_ref,'EdgeColor','k','FaceColor','none')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('mLF alg. 2 (iterative)','Interpreter','latex')
axis tight, zlim([min(tab_ref(:)) max(tab_ref(:))]), view(-30,40)
subplot(2,2,4); hold on, grid on, axis tight
imagesc(log10(abs(tab_ref-tab_app2)/max(abs(tab_ref(:)))),'XData',x1,'YData',x2)
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('{\bf log}(abs. err.)/max.','Interpreter','latex')
colorbar,
%
for i = 1:infoCas.n;  infoCas.name = strrep(infoCas.name,['\var{' num2str(i) '}'],['x_{' num2str(i) '}']); end
sgtitle(infoCas.name,'Interpreter','latex')
