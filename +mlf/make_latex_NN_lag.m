function latexList = make_latex_NN_lag(p_c,c)

n = length(p_c);
for ii = 1:n
    k(ii) = length(p_c{ii});
end
sumK    = sum(k);
prodK   = prod(k);
distY   = 2;
distX   = 4;

%
latexList   = [];
latexList   = [latexList '\begin{tikzpicture}[line width=0.4mm]'];
latexList   = [latexList '\tikzstyle{place}=[circle, draw=black, minimum size = 8mm]'];
latexList   = [latexList '\tikzstyle{placeInOut}=[circle, draw=orange, minimum size = 8mm]'];
% Input
%kk = 0;
for ii = 1:n
    %kk          = k(ii) + kk;
    latexList   = [latexList ['\node at (0,' num2str(-distY*ii) ') [placeInOut] (first_' num2str(ii) '){$\var{' num2str(ii) '}$};']];
end
% Hidden 1
kk = 0;
for ii = 1:n
    for jj = 1:k(ii)
        kk = kk + 1;
        latexList   = [latexList ['\node at (' num2str(distX) ',' num2str(-kk*distY) ') [place] (secondL' num2str(ii) '_' num2str(jj) '){$\frac{1}{\var{' num2str(ii) '}-\lani{' num2str(ii) '}{' num2str(jj) '}}$};']];
    end
end
% Hidden 2
for ii = 1:prodK
    latexList   = [latexList ['\node at (' num2str(2*distX) ',' num2str(-ii*distY) ') [place] (third_' num2str(ii) '){$\prod$};']];
end
% Output
latexList   = [latexList ['\node at (' num2str(3*distX) ',' num2str(-prodK-distY/2) ') [placeInOut] (output){$\bSigma$};']];

% Input -> Hidden 1
kk = 0;
for ii = 1:n
    for jj = 1:k(ii)
        kk = kk + 1;
        latexList   = [latexList ['\draw[->] (first_' num2str(ii) ')--(secondL' num2str(ii) '_' num2str(jj) ') node[above,sloped,pos=0.75] { };']];
    end
end
% Hidden 1 -> Hidden 2
combi = mlf.combinations_dim(k);
for ii = 1:n
    kk = 0;
    for jj = 1:prodK 
        kk = kk + 1;
        latexList   = [latexList ['\draw[->] (secondL' num2str(ii) '_' num2str(combi(jj,ii)) ')--(third_' num2str(kk) ') node[above,sloped,pos=0.25] {};']];
    end
end
% Hidden 2 -> Out
for ii = 1:prodK
    latexList   = [latexList ['\draw[->] (third_' num2str(ii) ')--(output) node[above,sloped,pos=0.25] {' num2str(c(ii)) '};']];
end
%
latexList   = [latexList '\end{tikzpicture}'];