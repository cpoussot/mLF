function [H,info] = examples(CAS)

% default
rng(1)
rnd     = 0;
Nip     = 20;
xbnd    = [-1 1];
mpts    = 1; % plot factor multiplier
%
tag = {};
%
switch CAS
    case 1 % ReLU
        n       = 2;
        H       = @(x) 1/2*(x(:,1) + abs(x(:,1))) + 1/10*x(:,2);
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\mathrm{ReLU}(\var{1})+\frac{1}{100}\var{2}$';
        dom     = 'R';
        tag     = {'irrational' '$C^0$'};
        %
        xbnd    = {[-1 1],[-1 -1e-10]};
    case 2
        n       = 2;
        H       = @(x) exp( sin(x(:,1)) + x(:,2).^2 );
        ref     = 'L/al. 2024';
        cite    = '\cite{Liu:2025}';
        nam     = '$\mathrm{exp}\left(\sin(\var{1}) + \var{2}^2\right)$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
    case 3
        n       = 2;
        H       = @(x) x(:,1).*x(:,2);
        ref     = 'L/al. 2024';
        cite    = '\cite{Liu:2025}';
        nam     = '$\var{1}\cdot \var{2}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
    case 4
        n       = 3;
        H       = @(x) 1/n*sum(sin(pi*x/2).^2,2);
        ref     = 'L/al. 2024';
        cite    = '\cite{Liu:2025}';
        nam     = '$\frac{1}{3} \sum_{i=1}^3 \sin(\pi x_i/2)^2$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
    case 5
        n       = 4;
        H       = @(x) exp(.5*(sin(pi*(x(:,1).^2+x(:,2).^2)) + ...
                               sin(pi*(x(:,3).^2+x(:,4).^2)) ) );
        ref     = 'L/al. 2024';
        cite    = '\cite{Liu:2025}';
        nam     = '$\mathrm{exp}\left(1/2 \left( \sin(\pi(\var{1}^2+\var{2}^2) + \sin(\pi(\var{3}^2+\var{4}^2) \right) \right)$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        Nip     = 20;%10;
    case 6
        n       = 2;
        H       = @(x) exp(x(:,1).*x(:,2))./((x(:,1).^2-1.44).*(x(:,2).^2-1.44));
        ref     = 'A/al. 2021 (A.5.1)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\mathrm{exp}\left(\var{1} \var{2}\right)}{(\var{1}^2-1.44)(\var{2}^2-1.44)}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
    case 7
        n       = 2;
        H       = @(x) log2(2.25-x(:,1).^2-x(:,2).^2);
        ref     = 'A/al. 2021 (A.5.2)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\mathrm{log}(2.25-\var{1}^2-\var{2}^2)$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
    case 8
        n       = 2;
        H       = @(x) tanh(4*(x(:,1)-x(:,2)));
        ref     = 'A/al. 2021 (A.5.3)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\mathrm{tanh}(4(\var{1}-\var{2}))$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        Nip     = 37;
    case 9
        n       = 2;
        H       = @(x) exp(-(x(:,1).^2+x(:,2).^2)/1000);
        ref     = 'A/al. 2021 (A.5.4)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\mathrm{exp}(\frac{-(\var{1}^2+\var{2}^2)}{1000})$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
    case 10
        n       = 2;
        H       = @(x) (abs(x(:,1)-x(:,2))).^3;
        ref     = 'A/al. 2021 (A.5.5)';
        cite    = '\cite{Austin:2021}';
        nam     = '$|\var{1}-\var{2}|^3$';
        dom     = 'R';
        tag     = {'irrational' '$C^0$'};
        %
        Nip     = 41;
    case 11
        n       = 2;
        H       = @(x) (x(:,1)+x(:,2).^3)./(x(:,1).*x(:,2).^2+2);
        ref     = 'A/al. 2021 (A.5.6)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}+\var{2}^3}{\var{1}\var{2}^2+2}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [1e-10 1];
    case 12
        n       = 2;
        H       = @(x) (x(:,1).^2+x(:,2).^2+x(:,1)-x(:,2)-1)./((x(:,1)-1.1).*(x(:,2)-1.1));
        ref     = 'A/al. 2021 (A.5.7)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}^2+\var{2}^2+\var{1}-\var{2}-1}{(\var{1}-1.1)(\var{2}-1.1)}$';
        dom     = 'R';
        tag     = {'rational'};
    case 13
        n       = 2;
        H       = @(x) (x(:,1).^4+x(:,2).^4+x(:,1).^2.*x(:,2).^2+x(:,1).*x(:,2))./((x(:,1)-1.1).*(x(:,2)-1.1));
        ref     = 'A/al. 2021 (A.5.8)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}^4+\var{2}^4+\var{1}^2\var{2}^2+\var{1}\var{2}}{(\var{1}-1.1)(y_2-1.1)}$';
        dom     = 'R';
        tag     = {'rational'};
    case 14
        n       = 4;
        H       = @(x) (x(:,1).^2+x(:,2).^2+x(:,1)-x(:,2)+1)./((x(:,3)-1.5).*(x(:,4)-1.5));
        ref     = 'A/al. 2021 (A.5.9)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}^2+\var{2}^2+\var{1}-\var{2}+1}{(\var{3}-1.5)(\var{4}-1.5)}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        Nip     = 10;
    case 15 
        n       = 2;
        H       = @(x) (x(:,1).^2+x(:,2).^2+x(:,1)-x(:,2)-1)./(x(:,1).^3+x(:,2).^3+4);
        ref     = 'A/al. 2021 (A.5.10)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}^2+\var{2}^2+\var{1}-\var{2}-1}{\var{1}^3+\var{2}^3+4}$';
        dom     = 'R';
        tag     = {'rational'};
    case 16
        n       = 2;
        H       = @(x) (x(:,1).^3+x(:,2).^3)./(x(:,1).^2+x(:,2).^2+3);
        ref     = 'A/al. 2021 (A.5.11)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}^3+\var{2}^3}{\var{1}^2+\var{2}^2+3}$';
        dom     = 'R';
        tag     = {'rational'};
    case 17
        n       = 2;
        H       = @(x) (x(:,1).^4+x(:,2).^4+x(:,1).^2.*x(:,2).^2+x(:,1).*x(:,2))./(x(:,1).^2.*x(:,2).^2-2.*x(:,1).^2-2.*x(:,2).^2+4);
        ref     = 'A/al. 2021 (A.5.12)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}^4+\var{2}^4+\var{1}^2\var{2}^2+\var{1}\var{2}}{\var{1}^2\var{2}^2-2\var{1}^2-2\var{2}^2+4}$';
        dom     = 'R';
        tag     = {'rational'};
    case 18
        n       = 2;
        H       = @(x) (x(:,1).^3+x(:,2).^3)./(x(:,1).^2.*x(:,2).^2-2.*x(:,1).^2-2.*x(:,2).^2+4);
        ref     = 'A/al. 2021 (A.5.13)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}^3+\var{2}^3}{\var{1}^2\var{2}^2-2\var{1}^2-2\var{2}^2+4}$';
        dom     = 'R';
        tag     = {'rational'};
    case 19
        n       = 2;
        H       = @(x) (x(:,1).^4+x(:,2).^4+x(:,1).^2.*x(:,2).^2+x(:,1).*x(:,2))./(x(:,1).^3+x(:,2).^3+4);
        ref     = 'A/al. 2021 (A.5.14)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\var{1}^4+\var{2}^4+\var{1}^2\var{2}^2+\var{1}\var{2}}{\var{1}^3+\var{2}^3+4}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        %Nip     = 51;
    case 20 
        n       = 3;
        H0      = @(E,G,M) ( 2.*sqrt(2).*M.*G.*sqrt(M.^2.*(M.^2+G.^2)) ) ./ ... 
                           ( (pi.*sqrt(M.^2+sqrt(M.^2.*(M.^2+G.^2)))).*((E.^2-M.^2).^2+M.^2.*G.^2) );
        H       = @(x) H0(x(:,1),x(:,2),x(:,3));
        ref     = 'A/al. 2021 (A.5.15)';
        cite    = '\cite{Austin:2021}';
        nam     = 'Breit Wigner function';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        e_val   = [80 100];
        g_val   = [5 10];
        m_val   = [90 93];
        xbnd    = {e_val g_val m_val};
    case 21
        n       = 4;
        H       = @(x) (atan(x(:,1))+atan(x(:,2))+atan(x(:,3))+atan(x(:,4)))./(x(:,1).^2.*x(:,2).^2-x(:,1).^2-x(:,2).^2+1);
        ref     = 'A/al. 2021 (A.5.16)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\sum_{i=1}^4\mathrm{atan}(x_i)}{\var{1}^2\var{2}^2-\var{1}^2-\var{2}^2+1}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        Nip     = 10;
        xbnd    = [-1 1]*.95;
    case 22
        n       = 4;
        H       = @(x) exp(x(:,1).*x(:,2).*x(:,3).*x(:,4))./(x(:,1).^2+x(:,2).^2-x(:,3).*x(:,4)+3);
        ref     = 'A/al. 2021 (A.5.17)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\frac{\mathrm{exp}(\var{1}\var{2}\var{3}\var{4})}{\var{1}^2+\var{2}^2-\var{3}\var{4}+3}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        Nip     = 10;
    case 23
        n       = 4;
        H       = @(x) 10*(sin(x(:,1)).*sin(x(:,2)).*sin(x(:,3)).*sin(x(:,4)))./(x(:,1).*x(:,2).*x(:,3).*x(:,4));
        %
        ref     = 'A/al. 2021 (A.5.18)';
        cite    = '\cite{Austin:2021}';
        nam     = '$10\prod_{i=1}^4\mathrm{sinc}(x_i)$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        Nip     = 11;
        xbnd    = [1e-6 4*pi];
    case 24
        n       = 2;
        H       = @(x) 10*(sin(x(:,1)).*sin(x(:,2)))./(x(:,1).*x(:,2));
        %
        ref     = 'A/al. 2021 (A.5.19)';
        cite    = '\cite{Austin:2021}';
        nam     = '$10\mathrm{sinc}(\var{1})\mathrm{sinc}(\var{2})$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        Nip     = 21;
        xbnd    = [1e-6 4*pi];
    case 25
        n       = 2;
        H       = @(x) x(:,1).^2+x(:,2).^2+x(:,1).*x(:,2)-x(:,2)+1;
        %
        ref     = 'A/al. 2021 (A.5.20)';
        cite    = '\cite{Austin:2021}';
        nam     = '$\var{1}^2+\var{2}^2+\var{1}\var{2}-\var{2}+1$';
        dom     = 'R';
        tag     = {'polynomial'};
    case 26
        n       = 3;
        H       = @(x) (x(:,1)+x(:,2)+x(:,3)) ./ ( 2*3+cos(x(:,1))+cos(x(:,2))+cos(x(:,3)) );
        %
        ref     = 'B/G 2025';
        cite    = '\cite{Balicki:2025}';
        nam     = '$\frac{\var{1}+\var{2}+\var{3}}{6+\cos(\var{1})+\cos(\var{2})+\cos(\var{3})}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        Nip     = 30;
        xbnd    = [-1 1]*10;
    case 27
        n       = 5;
        H       = @(x) (x(:,1)+x(:,2)+x(:,3)+x(:,4)+x(:,5)) ./ ( 2*5+cos(x(:,1))+cos(x(:,2))+cos(x(:,3))+cos(x(:,4))+cos(x(:,5)) );
        %
        ref     = 'B/G 2025';
        cite    = '\cite{Balicki:2025}';
        nam     = '$\frac{\var{1}+\var{2}+\var{3}+\var{4}+\var{5}}{10+\cos(\var{1})+\cos(\var{2})+\cos(\var{3})+\cos(\var{4})+\cos(\var{5})}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        Nip     = 13;
        xbnd    = [-1 1]*4;
    case 28
        n       = 2;
        H       = @(x) (x(:,1)./(x(:,1)+1)).^4 .* (1+exp(-x(:,2).^2)) .* (1+x(:,2).*cos(x(:,2)).*exp(-x(:,1).*x(:,2)./(x(:,1)+1)));
        %
        ref     = 'J/al. 2024 (Toy function)';
        cite    = '[none]';
        nam     = '$\left(\frac{\var{1}}{\var{1}+1}\right)^4 (1+\mathrm{exp}(-\var{2}^2)) \left(1+\var{2} \cos(\var{2}) \mathrm{exp}\frac{(-\var{1}\var{2})}{\var{1}+1}\right)$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        xbnd    = [1e-10 10];
        Nip     = 31;
    case 29
        n       = 2;
        H       = @(x) min(abs(x(:,1)).*10,1).*sign(x(:,1)) + (x(:,1).*x(:,2).^3)./10;
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\min(10|\var{1}|,1)\mathrm{sign}(\var{1}) + \frac{\var{1}\var{2}^3}{10}$';
        dom     = 'R';
        tag     = {'irrational' '$C^0$'};
    case 30
        n       = 8;
        frac1   = @(Tu,Hu,Hl) 2 * pi .* Tu .* (Hu-Hl);
        frac2a  = @(rw,r,Tu,L,Kw) 2.*L.*Tu ./ (log(r./rw).*rw.^2.*Kw);
        frac2b  = @(Tu,Tl) Tu ./ Tl;
        frac2   = @(rw,r,Tu,Tl,L,Kw) log(r./rw) .* (1+frac2a(rw,r,Tu,L,Kw)+frac2b(Tu,Tl));
        H0      = @(rw,r,Tu,Hu,Tl,Hl,L,Kw) frac1(Tu,Hu,Hl) ./ frac2(rw,r,Tu,Tl,L,Kw);
        H       = @(x) H0(x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6),x(:,7),x(:,8));
        %
        ref     = 'Borehole function';
        cite    = '\url{sfu.ca/~ssurjano}';
        nam     = '$\begin{array}{c}\mathrm{f}(r_w,r,T_u,H_u,T_l,H_l,L,K_w) =\\ \frac{2\pi T_u(H_u-H_l)}{\ln\left(\frac{r}{r_w}\right) \left(1+\frac{2LT_u}{\ln(r/r_w)r_w^2K_w}\right) + \frac{T_u}{T_l}} \end{array}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        rw      = [0.05, 0.15]; %	radius of borehole (m)
        r       = [100, 50000]; %	radius of influence (m)
        Tu      = [63070, 115600]; %   	transmissivity of upper aquifer (m2/yr)
        Hu      = [990, 1110]; % potentiometric head of upper aquifer (m)
        Tl      = [63.1, 116]; % transmissivity of lower aquifer (m2/yr)
        Hl      = [700, 820]; % potentiometric head of lower aquifer (m)
        L       = [1120, 1680]; % length of borehole (m)
        Kw      = [9855, 12045]; % hydraulic conductivity of borehole (m/yr)
        %
        xbnd    = {rw, r, Tu, Hu, Tl, Hl, L, Kw}; 
        Nip     = 4;
    case 31 
        n       = 6;
        H       = @(x) x(:,1).^2.*x(:,2).^3. + x(:,3).*x(:,4) - x(:,5).^2 + x(:,6);% + x(:,7);
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\var{1}^2 \var{2}^3 \var{3} \var{4} - \var{5}^2 + \var{6}$';
        dom     = 'R';
        tag     = {'polynomial'};
        %
        Nip     = 8;
        xbnd    = [-2 2];
    case 32
        n       = 2;
        H       = @(x) atan(x(:,1)) + x(:,2).^3;
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\texttt{atan}(\var{1}) + \var{2}^3$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        xbnd    = [-2 2];
    case 33
        n       = 2;
        H       = @(x) (x(:,1)+x(:,2))./( cos(x(:,1)).^2+cos(x(:,2))+3 );
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{\var{1}+\var{2}}{\cos(\var{1})^2+\cos(\var{2}) + 3}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        Nip     = 30;
        xbnd    = [-10 10];
    case 34
        n       = 2;
        N       = 1e3;
        % % Brute force Trefethen's way
        % zetF    = @(z_val) sum((N:-1:1).^(-z_val(:)),2); 
        % Continuation: 1/(1-2^(1-s))*(1-1/2^s+1/3^s-1/4^s+1/5^s ...)
        zetF    = @(z_val) (1./(1-2.^(1-z_val(:)))).*(sum((-1).^((N:-1:1)-1)./(N:-1:1).^(z_val(:)),2));
        H       = @(x) real(zetF(x(:,1)+1i*x(:,2)));
        %
        ref     = 'Riemann $\zeta$ function (real part)';
        cite    = '[none]';
        nam     = '$Re(\mathbf{\zeta}(\var{1}+\imath \var{2}))$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        xbnd    = {[.45 .55] [1 50]}; 
        Nip     = 200;
        mpts    = 4;
    case 35
        n       = 2;
        N       = 1e3;
        % % Brute force Trefethen's way
        % zetF    = @(z_val) sum((N:-1:1).^(-z_val(:)),2); 
        % Continuation: 1/(1-2^(1-s))*(1-1/2^s+1/3^s-1/4^s+1/5^s ...)
        zetF    = @(z_val) (1./(1-2.^(1-z_val(:)))).*(sum((-1).^((N:-1:1)-1)./(N:-1:1).^(z_val(:)),2));
        H       = @(x) imag(zetF(x(:,1)+1i*x(:,2)));
        %
        ref     = 'Riemann $\zeta$ function (imaginary part)';
        cite    = '[none]';
        nam     = '$Im(\mathbf{\zeta}(\var{1}+\imath \var{2}))$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        xbnd    = {[.45 .55] [1 50]}; 
        Nip     = 200;
        mpts    = 4;
    case 36
        n       = 3;
        H       = @(x) x(:,2)./(3+1/3*x(:,2).*x(:,1)-x(:,3).^2);
        %H       = @(x) x(:,2).*x(:,1)-x(:,3);
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{\var{2}}{3+1/3 \var{2}\var{1}-\var{3}^2}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [1/10 1];
        Nip     = 10;
    case 37
        n       = 4;
        H       = @(x) x(:,1).*x(:,4).^3+sin(2.*x(:,2)).*x(:,3);
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\var{1}\var{4}^3+\sin(2\var{2})\var{3}$';
        dom     = 'R';
        tag     = {'irrational' '$C^\infty$'};
        %
        xbnd    = [1e-3 1];
        Nip     = 10;%7;
        % n       = 5;
        % H       = @(x) (2+2.*x(:,3)/(2*pi)) .* (atan( 20*(x(:,1)-.5+x(:,2)/6).*exp(x(:,5))) + pi/2) + ...
        %                (2+2.*x(:,4)/(2*pi)) .* (atan( 20*(x(:,1)-.5-x(:,2)/6).*exp(x(:,5))) + pi/2);
        % %
        % ref     = 'P/P 2025';
        % cite    = '\cite{Poluektov:2025}';
        % nam     = '$\frac{2+2\var{3}}{2\pi} \mathbf{atan}( 20(\var{1}-.5+\var{2}/6)\mathbf{exp}(\var{5})) + \pi/2) + \frac{2+2\var{4}}{2\pi} \mathbf{atan}( 20(\var{1}-.5-\var{2}/6)\mathbf{exp}(\var{5})) + \pi/2) $';
        % dom     = 'R';
        % %
        % xbnd    = [0 1];
        % Nip     = 10;%7;
    case 38
        n       = 3;
        H       = @(x) (x(:,1).^9.*x(:,2).^7 + x(:,1).^3 + 5*x(:,3).^2 ) ./ ...
                       (5*x(:,1).^4 + 4*x(:,1).^2 + x(:,3).*x(:,2).^3 + 1);
        %
        ref     = 'A.C. Antoulas presentation';
        cite    = '[none]';
        nam     = '$\frac{\var{1}^9 \var{2}^7 + \var{1}^3 + 5 \var{3}^2}{5 \var{1}^4 + 4 \var{1}^2 + \var{3}\var{2}^3 + 1}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [-1 1]*1.1;
        Nip     = 30;
    case 39 
        n       = 3;
        H       = @(x) (x(:,3)+x(:,1).^4)./(x(:,1).^3+x(:,2).^2+1);
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{\var{3}+\var{1}^4}{\var{1}^3+\var{2}^2+1}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [.1 10];
        Nip     = 20;
    case 40
        n       = 4;
        H       = @(x) x(:,3).*x(:,1)./(x(:,1).^2+x(:,2)+x(:,3).^2+1)+x(:,4).^3;
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{\var{3}\var{1}}{\var{1}^2+\var{2}+\var{3}^2+1}+\var{4}^3$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [1 4];
        Nip     = 20;
    case 41
        n       = 5;
        H       = @(x) (x(:,5).^3.*x(:,3).*x(:,1)+x(:,3).^2)./(x(:,1).^3+x(:,2).*x(:,3)+x(:,4));
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{\var{5}^3\var{3}\var{1}+\var{3}^2}{\var{1}^3+\var{2}\var{3}+\var{4}}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [.1 1];
        Nip     = 5;
    case 42
        n       = 6;
        H       = @(x) (x(:,1)+x(:,3)-sqrt(2).*x(:,6).^2)./(x(:,1).^4+x(:,2).*x(:,3)+x(:,4).^3+x(:,5).^2+x(:,6));
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{\var{1}+\var{3}-\sqrt{2}\var{6}^2}{\var{1}^4+\var{2}\var{3}+\var{4}^3+\var{5}^2+\var{6}}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [.1 1];
        Nip     = 5;
    case 43
        n       = 7;
        H       = @(x) (x(:,3).*x(:,2).^3+1)./(x(:,1).^4+x(:,2).^2.*x(:,3)+x(:,4).^2+x(:,5)+x(:,6).^3+x(:,7));
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{\var{3}\var{2}^3+1}{\var{1}^4+\var{2}^2\var{3}+\var{4}^2+\var{5}+\var{6}^3+\var{7}}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [1 10];
        Nip     = 5;
    case 44
        n       = 8;
        H       = @(x) 1./(x(:,1).^4+x(:,2).^2.*x(:,3)+x(:,4).^2+x(:,5)+x(:,6)+x(:,7)+x(:,8) );
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{1}{\var{1}^4+\var{2}^2\var{3}+\var{4}^2+\var{5}+\var{6}+\var{7}+\var{8}}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [.1 20];
        Nip     = 5;
    case 45
        n       = 9;
        H       = @(x) 1./(x(:,1).^2+x(:,2).^2.*x(:,3)+x(:,4).^2+x(:,5)+x(:,6)+x(:,7)+x(:,8)+x(:,9));
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{1}{\var{1}^2+\var{2}^2\var{3}+\var{4}^2+\var{5}+\var{6}+\var{7}+\var{8}+\var{9}}$';
        dom     = 'R';
        tag     = {'rational'};
        % 
        xbnd    = [1 5];
        Nip     = 3;
    case 46
        n       = 10;
        H       = @(x) 1./(x(:,1)+x(:,1).^2.*x(:,2).*x(:,3)+x(:,4)+x(:,5)+x(:,6)+x(:,7).*x(:,8)+x(:,9).^2+x(:,10));
        %
        ref     = 'Personal communication';
        cite    = '[none]';
        nam     = '$\frac{1}{\var{1}+\var{1}^2\var{2}\var{3}+\var{4}+\var{5}+\var{6}+\var{7}\var{8}+\var{9}^2+\var{10}}$';
        dom     = 'R';
        tag     = {'rational'};
        %
        xbnd    = [1 5];
        Nip     = 3;
    % case 11
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11) s1/(s1^6+s2*s3+s4+s5+s6+s7*s8+s9+s10+s11);
    %     ord         = [6 1 1, 1 1 1, 1 1 1, 1 1];
    %     lam         = [1 5];
    %     N           = 2;
    %     %
    %     name_n      = 11;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex';
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 12
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12) 1/(s1^6+s1^2*s2*s3+s4+s5+s6+s7*s8+s9+s10+s11+s12);
    %     ord         = [6 1 1, 1 1 1, 1 1 1, 1 1 1];
    %     lam         = [1 5];
    %     N           = 2;
    %     %
    %     name_n      = 12;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex';
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 13
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13) (3*s1+4*s8+s12)/(s1^7+s2^2*s3+s4+s5^3+s6+s7*s8+s9*s10*s11+s13^2); 
    %     ord         = [7 2 1, 1 3 1, 1 1 1, 1 1 1, 2];
    %     lam         = [.1 10];
    %     %
    %     name_n      = 13;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 14
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14) (3*s1+4*s8+s12+3*s14^3)/(s1^7+s2^2*s3+s4+s5^3+s6+s7*s8+s9*s10*s11+s13^2+2); 
    %     ord         = [7 2 1, 1 3 1, 1 1 1, 1 1 1, 2 3];
    %     lam         = [.1 10];
    %     %
    %     name_n      = 14;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 15
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15) (3*s1+4*s8+s12+3*s14^3)/(s1^8+s2^2*s3+s4+s5^3+s6+s7*s8+s9*s10*s11+s13^2+2+s15^5); 
    %     ord         = [8 2 1, 1 3 1, 1 1 1, 1 1 1, 2 3 5];
    %     lam         = [.1 10];
    %     %
    %     name_n      = 15;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 16
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16) (3*s1^8+4*s8+s12+s13*s14+s15)/(s1+s2^2*s3+s4+s5+s6+s7*s8+s9*s10*s11+s13+s16^3*pi); 
    %     ord         = [8 2 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1, 3];
    %     lam         = [.1 2];
    %     %
    %     name_n      = 16;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 17
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17) (3*s1^9+4*s8+s12+s13*s14+s15)/(s1+s2^2*s3+s4+s5+s6+s7*s8+s9*s10*s11+s13+s16^3*pi+s17^2);
    %     ord         = [9 2 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1, 3 2];
    %     lam         = [.1 2];
    %     %
    %     name_n      = 17;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 18
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18) (3*s1^9+4*s8+s12+s13*s14+s15+s18^2)/(s1+s2^2*s3+s4+s5+s6+s7*s8+s9*s10*s11+s13+s16^3*pi+s17^2);
    %     ord         = [9 2 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1, 3 2 2];
    %     lam         = [.1 2];
    %     %
    %     name_n      = 18;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 19
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19) (3*s1^10+4*s8+s12+s13*s14+s15+s18^2)/(s1+s2^2*s3+s4+s5+s6+s7*s8+s9*s10*s11+s13+s16^3*pi+s17^2+s19);
    %     ord         = [10 2 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1, 3 2 2, 1];
    %     lam         = [.1 2];
    %     %
    %     name_n      = 19;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 20
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20) (3*s1^3+4*s8+s12+s13*s14+s15)/(s1^10+s2^2*s3+s4+s5+s6+s7*s8+s9*s10*s11+s13+s16^3*pi+s17+s18*s19-s20); 
    %     ord         = [10 2 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1, 3 1 1, 1 1];
    %     lam         = [.1 2];
    %     %
    %     name_n      = 20;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
    % case 30
    %     H           = @(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24,s25,s26,s27,s28,s29,s30) (3*s1^3+4*s8+s12+s13*s14+s15+s21^4-s30*s2)/(s1+s2^2*s3+s4+s5+s6+s7*s8+s9*s10*s11+s13+s16^3*pi+s17+s18*s19-s20+s22*s23-s24/2+s25*s26-sqrt(2)*s27-s28+s29+s30^3); 
    %     ord         = [3 2 1, 1 1 1, 1 1 1, 1 1 1, 1 1 1, 3 1 1, 1 1 4, 1 1 1, 1 1 1, 1 1 3];
    %     lam         = [.1 5];
    %     %
    %     name_n      = 30;
    %     name_type   = 'Rational function';
    %     name_domain = 'Complex'; 
    %     name_rmk    = '';
    %     name_ref    = '';
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
ord     = ones(1,n)*Nip;
p_c     = cell(1,n);
p_r     = cell(1,n);
ip      = cell(1,n);
xmin    = inf;
xmax    = -inf;
for ii = 1:n
    xmin    = min(xmin,xlim{ii}(1));
    xmax    = max(xmax,xlim{ii}(2));
    p_c{ii} = linspace(xlim{ii}(1),xlim{ii}(2),ord(ii));
    dx      = abs(p_c{ii}(end)-p_c{ii}(end-1))/2;
    p_r{ii} = (1+rnd)*p_c{ii}+dx;
    ip{ii}  = [p_c{ii}(:); p_r{ii}(:)];
end

%%% 
info.n      = n;
info.ref    = ref;
info.cite   = cite;
info.name   = nam;
info.domain = dom;
info.tag    = tag;
info.nip    = Nip;
info.mpts   = mpts;
info.bound  = xlim;
info.xmin   = xmin;
info.xmax   = xmax;
info.ip     = ip;
info.p_c    = p_c;
info.p_r    = p_r;
info.tab_MB = mlf.file_size(2*ord);