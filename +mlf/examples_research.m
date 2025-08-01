% param.x_bnd
% param.x_grid
% param.x_Npts


function [tab,xpt,info] = examples_research(CAS,param)

if nargin < 2
    param = [];
end
%%% Default "function"
p       = 1;
gamma   = 1;
if isa(param,'struct')
    if isfield(param,'p')
        p = param.p;
    end
    if isfield(param,'gamma')
        gamma = param.gamma;
    end
end
%
switch CAS
    case '1-1'
        n       = 1;
        H       = @(x) 1+0.*x;
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$1$';
        dom     = 'R';
        tag     = {'constant' '$C^\infty$'};
    case '1-2'
        n       = 1;
        H       = @(x) 1+gamma.*x(:,1).^p;
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$1+\gamma x_1^p$';
        dom     = 'R';
        tag     = {'polynomial' '$C^\infty$'};
    case '1-3'
        n       = 1;
        H       = @(x) 1./(1+gamma.*x(:,1).^p);
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{1}{1+\gamma x_1^p}$';
        dom     = 'R';
        tag     = {'rational' '$C^\infty$'};
    case '1-4'
        n       = 1;
        H       = @(x) (1-x(:,1).^p)./(1+gamma.*x(:,1).^p);
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{1-\gamma x_1^p}{1+\gamma x_1^p}$';
        dom     = 'R';
        tag     = {'rational' '$C^\infty$'};
    case '1-5'
        n       = 1;
        H       = @(x) x(:,1).^2 .* atan(1./x(:,1));
        ref     = 'Salazar Celis PhD thesis';
        cite    = '\cite{SalazarCelis:2008}';
        nam     = '$x_1^2 \mathbf{atan}\frac{1}{x_1}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
    case '2-1'
        n       = 2;
        H       = @(x) 1e-1*x(:,1)+gamma.*x(:,2).^p;
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$10^{-1}x_1+\gamma x_2^p$';
        dom     = 'R';
        tag     = {'polynomial' '$C^\infty$'};
    case '2-2'
        n       = 2;
        H       = @(x) x(:,1)+x(:,2);
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$x_1+x_2$';
        dom     = 'R';
        tag     = {'polynomial' '$C^\infty$'};
end

%%% Default "x"
x_grid  = @linspace;
for i = 1:n
    x_bnd{i}    = [-1 1];
    x_Npts(i)   = 10;
end
if isa(param,'struct')
    if isfield(param,'x_bnd')
        x_bnd = param.x_bnd;
    end
    if isfield(param,'x_grid')
        x_grid = str2func(param.x_grid);
    end
    if isfield(param,'x_Npts')
        x_Npts = param.x_Npts;
    end
end

%%% x-data
ip = cell(1,n);
for i = 1:n
    xmin    = min(x_bnd{i});
    xmax    = max(x_bnd{i});
    tmp     = x_grid(xmin,xmax,x_Npts(i));
    ip{i}   = tmp(:);
end

%%% Evaluation
[tab,xpt,dim] = mlf.build_tab(H,ip);

%%% 
info.n      = n;
info.H      = H;
info.p      = p;
info.gamma  = gamma;
%
info.bnd    = x_bnd;
info.ip     = ip;
info.dim    = dim;
%
info.tex    = nam;
info.tag    = tab;
info.ref    = ref;
info.cite   = cite;