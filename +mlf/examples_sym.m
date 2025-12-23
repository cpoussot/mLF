function [H,info] = examples_sym(CAS)

% default
Nip     = 20;
xbnd    = [-1 1];
%
tag = {};
%
syms s1 s2 s3 s4 s5
switch CAS
    case 1
        n       = 2;
        H       = s1*s2^2;
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\var{1} \var{2}^2$';
        dom     = 'R';
        tag     = {'polynomial'};
        %
        Nip     = 10;
        xbnd    = [-1 1];
    case 2
        n       = 3;
        H       = s1*s2+s1*s3+s2*s3;
        %
        ref     = 'G. P\''olya and G.Szeg\"o';
        cite    = '\cite{Polya:1925}';
        nam     = '$\var{1}\var{2}+\var{1}\var{3}+\var{2}\var{3}$';
        dom     = 'R';
        tag     = {'polynomial'};
        %
        xbnd    = [-1/2,1]; 
        Nip     = 6;
    case 3
        n       = 2;
        H       = (s1^3*s2^2-s1^2*s2^2+s1^2+1)/(s1+s2+s1^2*s2+s1^2+2*s1^3);
        %
        ref     = 'A.C. Antoulas';
        cite    = '[Personnal communication]';
        nam     = '$\var{1}$';
        dom     = 'R';
        tag     = {'polynomial'};
        %
        xbnd    = [-1/2,1]; 
        Nip     = 8;
    % case  
    %     n       = 5;
    %     H       = 12*s3 + 24.*s1.*s3 - 6.*s2.*s3 + 4.*s3.*s4 + 2.*s2.*s4 - 18.*s3.*s5 ...
    %               - 12.*s1.*s2.*s3 + 8.*s1.*s3.*s4 - 2.*s1.*s2.*s4 - 36.*s1.*s3.*s5 + 4.*s2.*s3.*s4 + 9.*s2.*s3.*s5 - 2.*s2.*s4.*s5 - 6.*s3.*s4.*s5 ...
    %               - 10.*s1.*s2.*s3.*s4 + 18.*s1.*s2.*s3.*s5 + 2.*s1.*s2.*s4.*s5 - 12.*s1.*s3.*s4.*s5 ... 
    %               + 12.*s1.*s2.*s3.*s4.*s5;
    %     %H       = @(x) (1 + 2*s1).*(-2 + s2).*(-s3).*(3 + s4).*(2 - 3*s5) + ... 
    %     %               (-1 + s1).*(2*s2).*(1 + 3*s3).*(-s4).*(1 - s5);
    %     %
    %     ref     = 'G/al. 2025 (Ex 3.1)';
    %     cite    = '\cite{GHK:2025}';
    %     nam     = '$\begin{array}{c}(1 + 2\var{1})(-2 + \var{2})(-\var{3})(3 + \var{4})(2- 3\var{5}) \\ + (-1 + \var{1})(2\var{2})(1 + 3\var{3})(-\var{4})(1 -\var{5})\end{array}$';
    %     dom     = 'R';
    %     tag     = {'polynomial'};
    %     %
    %     xbnd    = [-2,2]; 
    %     Nip     = 6;
end
%
if ~isa(xbnd,'cell')
    for ii = 1:n
        xlim{ii} = xbnd;
    end
else
    xlim = xbnd;
end
%
p_c     = cell(1,n);
p_r     = cell(1,n);
ip      = cell(1,n);
for ii = 1:n
    ip{ii}  = sym(linspace(xlim{ii}(1),xlim{ii}(2),Nip));
    p_c{ii} = ip{ii}(2:2:end);
    p_r{ii} = ip{ii}(1:2:end);
end

%%% 
info.n      = n;
info.ref    = ref;
info.cite   = cite;
info.name   = nam;
info.domain = dom;
info.tag    = tag;
info.nip    = Nip;
info.bound  = xlim;
info.ip     = ip;
info.p_c    = p_c;
info.p_r    = p_r;
info.tab_MB = mlf.file_size(Nip*ones(1,n));