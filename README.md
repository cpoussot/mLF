# Multivariate Loewner Framework

### Overview

The Multivariate Loewner Framework is introduced  by A.C. Antoulas, I-V. Gosea and C. Poussot-Vassal in ["On the Loewner framework, the Kolmogorov superposition theorem, and the curse of dimensionality"](https://arxiv.org/abs/2405.00495), to appear in SIAM Review. It allows constructing a $n$-variate rational function approximating a $n$-dimensional tensor (either real or complex valued). It is suited to approximate, from any $n$-dimensional tensor 
- $n$-variable static functions
- $n$-variable (parametric) dynamical systems

### Contributions claim


- We propose a ***generalized realization*** form for rational functions in n-variables (for any $n$), which are described in the Lagrange basis;
- We show that the $n$-dimensional Loewner matrix can be written as the solution of a ***series of cascaded Sylvester equations***;
- We demonstrate that the required variables to be determined, i.e. the barycentric coefficients, can be computed using a sequence of small-scale 1-dimensional Loewner matrices instead of the large-scale ($Q\times K$, $Q\geq K$)  $n$-dimensional one, therefore drastically ***taming the curse of dimensionality***, i.e. reducing both the computational effort and the memory needs, and, in addition improving accuracy;
- We show that this decomposition achieves variables decoupling; thus connecting the Loewner framework for rational interpolation of multivariate functions and the ***Kolmogorov Superposition Theorem (KST)***, restricted to rational functions. The result is the formulation of KST for the special case of rational functions;
- Connections with KAN neural nets follows (detailed in future work).


### Main reference

```
@article{AGPV:2025,
	Author	= {A.C. Antoulas and I-V. Gosea and C. Poussot-Vassal},
	Doi 	= {},
	Journal = {to appear in SIAM Review},
	Number 	= {},
	Pages 	= {},
	Title 	= {On the Loewner framework, the Kolmogorov superposition theorem, and the curse of dimensionality},
	Volume 	= {},
	Year 	= {},
    	Note 	= {\url{https://arxiv.org/abs/2405.00495}}, 
}
```

### Additional support material (presentations, reports, ...)

- A.C. Antoulas presentation [BANFF Video]([https://www.birs.ca/events/2025/5-day-workshops/25w5376/schedule](https://www.birs.ca/events/2025/5-day-workshops/25w5376/videos/watch/202504090859-Antoulas.html))
- C. Poussot-Vassal presentation [GT Identification video](https://youtu.be/M2SX3C4VCt8), [slides](https://drive.google.com/file/d/1qEirwD7c5h56h1gRTPJmmyJNQY2qa4B1/view?usp=sharing)
- Benchmark results and comparison (to come)


# The "mLF" MATLAB package 

The code (`+mlf` folder)  provided in this GitHub page is given for open science purpose. Its principal objective is to accompany the [paper](https://arxiv.org/abs/2405.00495) by the authors, thus aims at being educative rather than industry-oriented. Evolutions (numerical improvements) may come with time. Please, cite the reference above if used in your work and do not hesitate to contact us in case of bug of problem when using it. Below we present an example of use, then functions list are given.

## A simple MATLAB code example

Here is a simple code that describes how to deploy the cascaded 1-D Loewner null space construction. Refer to https://arxiv.org/abs/2405.00495 for notations and related equations. Code below is `demo.m`.

First add the path where the `+mlf` package is.

```Matlab
%addpath("location_of_mlf") % Add the location of the +mlf package
```

Then compute the barycentic coefficients for a given function `H`.

```Matlab
%%% Define a multivariate handle function 
n       = 3; % number of variables
H       = @(s1,s2,s3) (s3/100-1)*(s2-pi/2)*(s1+atan(2*s2)*tanh(5*(s2-pi)))/(s1^2+s3/10*cos(3*s1)+3)/(s2+10);
%H       = @(s1,s2,s3) s2/(s1-s3^2/2);
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
% Estimate order along each variables and select a subset of IP
ord                 = mlf.compute_order(p_c,p_r,tab,tol_ord,[],5,true);
[pc,pr,W,V,tab_red] = mlf.points_selection(p_c,p_r,tab,ord,false);
w                   = mlf.mat2vec(W);

%%% Recursive null space computation (Section 5, Theorem 5.3)
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
	% Evaluate the model in Lagrangian form at values = param
        tab_app(jj,ii)          = mlf.eval_lagrangian(pc,w,c_rec1,param,false);
    end
end
% Evaluation of the functions (original and reconstructed) 
figure
subplot(1,2,1); hold on, grid on
surf(X,Y,tab_app,'EdgeColor','none'), hold on
surf(X,Y,tab_ref,'EdgeColor','k','FaceColor','none')
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('Original vs. Approximation','Interpreter','latex')
axis tight, zlim([min(tab_ref(:)) max(tab_ref(:))]), view(-30,40)
% Evaluation of the error
subplot(1,2,2); hold on, grid on, axis tight
imagesc(log10(abs(tab_ref-tab_app)/max(abs(tab_ref(:)))),'XData',x1,'YData',x2)
xlim([min(x1) max(x1)])
ylim([min(x2) max(x2)])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
title('{\bf log}(abs. err.)/max.','Interpreter','latex')
colorbar,
```

## Functions description

### `mlf.check`

Function which checks that column `p_c` ($P_c^{n}$) and row `p_r` ($P_r^{n}$) interpolation points sets are disjoints. Note that both `p_c` and `p_r` are $n$-dimensional cells, where each `p_c{i}` (`i=1...n`) gathers the interpolatoin points along each variables.

#### Input

- `p_c`: column interpolation points $P_c^{n}$ ($n$-dimneiosnl cell with double)
- `p_r`: row interpolation points $P_r^{n}$ ($n$-dimneiosnl cell with double)

#### Output

- `ok`: tag assessing that interpolation points are disjoints (boolean)

#### Syntax 

```Matlab
ok = mlf.check(p_c,p_r)
```

### `mlf.make_tab`

Function which constructs the $n$-dimensional tensor $\mathbf{tab}_n$, from the column `p_c` ($P_c^{n}$) and row `p_r` ($P_r^{n}$) interpolation point sets. Note that both `p_c` and `p_r` are $n$-dimensional cells, where each `p_c{i}` (`i=1...n`) gathers the interpolation points along each variables.

#### Input

- `p_c`: column interpolation points $P_c^{n}$ ($n$-dimensional cell with double)
- `p_r`: row interpolation points $P_r^{n}$ ($n$-dimensional cell with double)
- `show`: tag used to show the advancement of the tableau construction (boolean)

#### Output

- `tab`: $n$-dimensional tensor $\mathbf{tab}_n$ ($n$-dimensional double)

#### Syntax 

```Matlab
tab = mlf.make_tab(H,p_c,p_r,show);
```

### `mlf.compute_order`

Function which estimates the orders along each variables.

#### Input

- `p_c`: column interpolation points $P_c^{n}$ ($n$-dimensional cell with double)
- `p_r`: row interpolation points $P_r^{n}$ ($n$-dimensional cell with double)
- `tab`: $n$-dimensional tensor $\mathbf{tab}_n$ ($n$-dimensional double)
- `ord_tol`: normalized singular value threshild for order selection (positive real scalar)
- `ord_obj`: maximal orders tolerated along each variables ($n\times 1$ integer)
- `ord_n`: maximal single variable Loewner sincular value computation per variables (integer scalar)
- `show`: tag used to show the advancement of the tableau construction (boolean)

#### Output

- `ord`: $n$-dimensional vector with order along each variables ($n\times 1$ double)

#### Syntax 

```Matlab
ord = mlf.compute_order(p_c,p_r,tab,ord_tol,ord_obj,ord_n,show);
```
### `mlf.points_selection`

Function which selects a subset of interpolation points, accordingly to the order along each variables.

#### Input

- `p_c`: column interpolation points $P_c^{n}$ ($n$-dimensional cell with double)
- `p_r`: row interpolation points $P_r^{n}$ ($n$-dimensional cell with double)
- `tab`: $n$-dimensional tensor $\mathbf{tab}_n$ ($n$-dimensional double)
- `ord`: $n$-dimensional vector with order along each variables ($n\times 1$ double)
- `square`: tag used to set wheather row interpolation points have the same dimension as columns or may be greather;  `true` leads to same number while `false` allows more row than columns (boolean)

#### Output

- `pc`: reduced column interpolation points $P_c^{n}$ ($n$-dimnensional cell with double)
- `pr`: reduced row interpolation points $P_r^{n}$ ($n$-dimensional cell with double)
- `W`: tensor corresponding to $P_r^{n}$ evaluation ($n$-dimensional cell with double)
- `V`: tensor corresponding to $P_c^{n}$ evaluation ($n$-dimensional cell with double)
- `tab_red`: reduced $n$-dimensional tensor $\mathbf{tab}_n$ ($n$-dimensional double)

#### Syntax 

```Matlab
[pc,pr,W,V,tab_red] = mlf.points_selection(p_c,p_r,tab,ord,square);
```

### `mlf.loewner_null_rec`

Function which construct the Loewner matrix null space (barycentric coefficients) using the recursive approach (without constructing the $N\times N$ Loewner matrix).

#### Input

- `pc`: reduced column interpolation points $P_c^{n}$ ($n$-dimensional cell with double)
- `pr`: reduced row interpolation points $P_r^{n}$ ($n$-dimensional cell with double)
- `tab_red`: reduced $n$-dimensional tensor $\mathbf{tab}_n$ ($n$-dimensional double)
- `method`: method used for null space computation (string); for the time being, we suggest to set this variable to `svd0`.

#### Output

- `c`: Loewner matrix null space, being the barycentric coefficients ($N \times 1$ double)
- `info`: information related to the process (e.g. \texttt{flop} count)

#### Syntax 

```Matlab
[c,info]  = mlf.loewner_null_rec(pc,pr,tab_red,method);
```

### `mlf.eval_lagrangian`

Function which evaluates the model in Lagrangian form.

#### Input

- `pc`: reduced column interpolation points $P_c^{n}$ ($n$-dimensional cell with double)
- `w`: original vectorized data ealuated at $P_c^{n}$ ($n$-dimensional cell with double). In practice we have `w = mlf.mat2vec(W);`.
- `c`: Loewner matrix null space, being the barycentric coefficients ($N \times 1$ double)
- `param`: $n$-variable vector where the the model has to be evaluated ($n\times 1$ double)
- `show`: tag used to show the advancement of the model evaluation (boolean)

#### Output

- `val`: value of the Lagrangian model (double)

#### Syntax 
```Matlab
val = mlf.eval_lagrangian(pc,w,c,param,false);
```

## Feedbacks

Please send any comment to C. Poussot-Vassal (charles.poussot-vassal@onera.fr) if you want to report any bug or user experience issues.

## Disclaimer

Once again, this deposit consitutes a research code that accompany the paper mentionned above. It is not aimed to be included in any third party software without the consent of the authors. Authors decline responsabilities in case of problem when applying the code.

Notice also that pathological cases may appear. A more advanced code, to deal with practical and theoretical issues/limitations is currently under developpement by the authors.
