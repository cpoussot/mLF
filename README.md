# Multivariate Loewner Framework - "mLF" MATLAB package 

## Overview

The Multivariate Loewner Framework is introduced  by A.C. Antoulas, I.V. Gosea and C. Poussot-Vassal in https://arxiv.org/abs/2405.00495. It allows constructing a $n$-variate rational function approximating a $n$-dimensional tensor (either real of complex valued). It is suited to approximate, from any $n$-dimensional tensor 
- $n$-variables static functions
- $n$-variables (parametric) dynamical systems

## Contributions claim

We claim the following innovations:

- We propose a ***generalized realization*** form for rational functions in n-variables (for any $n$), which are described in the Lagrange basis;
- We show that the $n$-dimensional Loewner matrix can be written as the solution of a ***series of cascaded Sylvester equations***;
- We demonstrate that the required variables to be determined, i.e. the barycentric coefficients, can be computed using a sequence of small-scale 1-dimensional Loewner matrices instead of the large-scale ($N\times N$) $n$-dimensional one, therefore drastically ***taming the curse of dimensionality***, i.e. reducing the both computational effort and memory needs, and improving accuracy;
- We show that this decomposition achieves variables decoupling; thus connecting the Loewner framework for rational interpolation of multivariate functions and the ***Kolmogorov Superposition Theorem (KST)***, restricted to rational functions. The result is the formulation of KST for the special case of rational functions;
- Connections with KAN neural nets follows (detailed in future work).

## Reference and additional support material (presentations, reports, ...)

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
	Bdsk-Url-1 = {https://arxiv.org/abs/2405.00495},
}
```

### Some presentations and reports

- A.C. Antoulas presentation at BANFF https://www.birs.ca/events/2025/5-day-workshops/25w5376/schedule
- C. Poussot-Vassal presentation 
- Benchmark result 


# The mLF MATLAB package

The code embedded in this GitHub page is given for open science perspectives. Please cite the reference above and do not hesitate to contact us in case of spotted bug.

## mlf.check

```Matlab
ok = mlf.check(p_c,p_r)
```

## mlf.make_tab

```Matlab
tab = mlf.make_tab(H,p_c,p_r,true);
```

## mlf.compute_order

```Matlab
ord = mlf.compute_order(p_c,p_r,tab,1e-10,[],5,true);
```
## mlf.points_selection

```Matlab
[p_c,p_r,W,V,tab_red] = mlf.points_selection(p_c,p_r,tab,ord,true);
```

## A simple MATLAB code example

Here is a simple code that describes how to deploy the cascaded 1-D Loewner null space construction.
Please refer to https://arxiv.org/abs/2405.00495  for notations and related equations.

```Matlab
%%% Define a multivariate handle function 
n = 2; % number of variables
H = @(s1,s2) (s1+s2)/(cos(s1)^2+cos(s2)+3);
%%% Interpolation points - separate columns and rows (as Section 3, eq. 13-15)
for ii = 1:n
    p_c{ii} = linspace(-1,1,10);
    dx      = abs(p_c{ii}(2)-p_c{ii}(1))/2;
    p_r{ii} = p_c{ii}+dx;
end
%%% Elementary ingredients
ok                      = mlf.check(p_c,p_r);
tab                     = mlf.make_tab(H,p_c,p_r,true);
ord                     = mlf.compute_order(p_c,p_r,tab,1e-10,[],5,true);
[p_c,p_r,W,V,tab_red]   = mlf.points_selection(p_c,p_r,tab,ord,true);

%%% Recursive null-space computation (using tab_n)
tic
[c_rec,info_rec]    = mlf.loewner_null_rec(p_c,p_r,tab_red,'svd0');
toc
%%% Recursive null-space computation (using tab_n)
%[c_rec,info_rec]    = mlf.loewner_null_recH(p_c,p_r,H,'svd0');
info_rec

%%% Evaluate the interpolated model vs. orginial

%%% Decoupling features, toward rational KST
[Var,Lag,Bary]  = mlf.decoupling(p_c,info_rec);

```

## Feedbacks

Please send any comment to [C. Poussot-Vassal](charles.poussot-vassal@onera.fr) if you report a bug or user experience issues.

## Claim

This deposit consitutes a research code and is not aimed to be included in any third party software without the consent of the authors. Authors decline responsabilities in case of problem when applying the code.

