# Multivariate Loewner Framework

### Overview

The Multivariate Loewner Framework is introduced  by A.C. Antoulas, I.V. Gosea and C. Poussot-Vassal in ["On the Loewner framework, the Kolmogorov superposition theorem, and the curse of dimensionality"](https://arxiv.org/abs/2405.00495). It allows constructing a $n$-variate rational function approximating a $n$-dimensional tensor (either real of complex valued). It is suited to approximate, from any $n$-dimensional tensor 
- $n$-variables static functions
- $n$-variables (parametric) dynamical systems

### Contributions claim

We claim the following innovations:

- We propose a ***generalized realization*** form for rational functions in n-variables (for any $n$), which are described in the Lagrange basis;
- We show that the $n$-dimensional Loewner matrix can be written as the solution of a ***series of cascaded Sylvester equations***;
- We demonstrate that the required variables to be determined, i.e. the barycentric coefficients, can be computed using a sequence of small-scale 1-dimensional Loewner matrices instead of the large-scale ($N\times N$) $n$-dimensional one, therefore drastically ***taming the curse of dimensionality***, i.e. reducing the both computational effort and memory needs, and improving accuracy;
- We show that this decomposition achieves variables decoupling; thus connecting the Loewner framework for rational interpolation of multivariate functions and the ***Kolmogorov Superposition Theorem (KST)***, restricted to rational functions. The result is the formulation of KST for the special case of rational functions;
- Connections with KAN neural nets follows (detailed in future work).


### Main reference

```
@article{AGPV:2025,
	Author = {A.C. Antoulas and I-V. Gosea and C. Poussot-Vassal},
	Doi = {},
	Journal = {SIAM Review},
	Number = {},
	Pages = {},
	Title = {On the Loewner framework, the Kolmogorov superposition theorem, and the curse of dimensionality},
	Volume = {},
	Year = {},
    note ={\url{https://arxiv.org/abs/2405.00495}}, 
}
```

### Additional support material (presentations, reports, ...)

- A.C. Antoulas presentation at [BANFF](https://www.birs.ca/events/2025/5-day-workshops/25w5376/schedule)
- C. Poussot-Vassal presentation at [GT Identification](https://drive.google.com/file/d/1qEirwD7c5h56h1gRTPJmmyJNQY2qa4B1/view?usp=sharing) with [YouTube video](https://youtu.be/M2SX3C4VCt8)
- Benchmark result 


# The "mLF" MATLAB package 

The code embedded in this GitHub page is given for open science perspectives. Please cite the reference above and do not hesitate to contact us in case of spotted bug.

## Functions description

### mlf.check

```Matlab
ok = mlf.check(p_c,p_r)
```

### mlf.make_tab

```Matlab
tab = mlf.make_tab(H,p_c,p_r,true);
```

### mlf.compute_order

```Matlab
ord = mlf.compute_order(p_c,p_r,tab,1e-10,[],5,true);
```
### mlf.points_selection

```Matlab
[p_c,p_r,W,V,tab_red] = mlf.points_selection(p_c,p_r,tab,ord,true);
```

## A simple MATLAB code example

Here is a simple code that describes how to deploy the cascaded 1-D Loewner null space construction.
Please refer to https://arxiv.org/abs/2405.00495  for notations and related equations.

```Matlab
%addpath("location_of_mlf") % Add the location of the +mlf package
%%% Define a multivariate handle function (here toy case)
n       = 3; % number of variables
H       = @(s1,s2,s3) (s3/100-1)*(s2-pi/2)*(s1+atan(2*s2)*tanh(5*(s2-pi)))/(s1^2+s3/10*cos(3*s1)+3)/(s2+10);
% /!\ important tuning variable when data generated from irrational function 
% used to set the approximation order 
tol_ord = 1e-7; 
% Interpolation points (IP) - separate columns and rows (as Section 3, eq. 13-15)
for ii = 1:n
    p_c{ii} = linspace(-10,10,30);
    dx      = abs(p_c{ii}(2)-p_c{ii}(1))/2;
    p_r{ii} = p_c{ii}+dx;
end

%%% Elementary ingredients
% Check that column/row IP are disjoint (Section 3)
ok                  = mlf.check(p_c,p_r);
% Construct tab_n (Section 3)
tab                 = mlf.make_tab(H,p_c,p_r,true);
% Estimate order along each variables and select a subset of IP
ord                 = mlf.compute_order(p_c,p_r,tab,tol_ord,[],5,true);
[pc,pr,W,V,tab_red] = mlf.points_selection(p_c,p_r,tab,ord,false);
w                   = mlf.mat2vec(W);

%%% Recursive null-space computation (Section 5, Theorem 5.3)
[c_rec1,info_rec1]  = mlf.loewner_null_rec(pc,pr,tab_red,'svd0');

%%% Curse of dimensionality: flops & memory estimation (Section 5, Theorem 5.5 & Theorem 5.6)
fprintf('FLOPS\n')
fprintf(' * recursive: %d\n',info_rec1.nflop)
fprintf(' * full: %d\n',length(c_rec1)^3)
fprintf('MEMORY\n')
fprintf(' * recursive: %d MB\n',max(ord+1)^2/2^20)
fprintf(' * full: %d MB\n',prod(ord+1)^2/2^20)
```

Evaluate simplified function and display results

```Matlab
%%% 3D plot 
% Along first and second variables 
% Other variables are randomly chosen betwwen bounds
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
        tab_app(jj,ii)          = mlf.eval_lagrangian(pc,w,c_rec1,param,false);
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
%
subplot(1,2,2); hold on, grid on, axis tight
imagesc(log10(abs(tab_ref-tab_app)/max(abs(tab_ref(:)))),'XData',x1,'YData',x2)
xlim([min(x1) max(x1)])
ylim([min(x2) max(x2)])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('{\bf log}(abs. err.)/max.','Interpreter','latex')
colorbar,
```

## Feedbacks

Please send any comment to [C. Poussot-Vassal](charles.poussot-vassal@onera.fr) if you report a bug or user experience issues.

## Claim

This deposit consitutes a research code and is not aimed to be included in any third party software without the consent of the authors. Authors decline responsabilities in case of problem when applying the code.

