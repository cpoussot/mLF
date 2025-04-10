# Multivariate Loewner Framework - mLF 
Multivariate Loewner Framework 


# Citation and feedbacks

The code below is given for open science perspectives. Please cite the reference below and do not hesitate to contact us in case of spotted bug.

## Reference

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

## Feedbacks

Please send any comment to charles.poussot-vassal@onera.fr if you report a bug or user experience issue.


# A simple MATLAB code

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
ok                      = data.check(p_c,p_r);
tab                     = data.make_tab(H,p_c,p_r,true);
ord                     = data.compute_order(p_c,p_r,tab,1e-10,[],5,true);
[p_c,p_r,W,V,tab_red]   = data.points_selection(p_c,p_r,tab,ord,true);

%%% Full null-space span computation (tableau)
tic
LL                  = data.loewnerMatrix(p_c,p_r,W,V);
c_full              = data.null(LL,'svd0');
toc 

%%% Recursive null-space computation (tableau)
tic
[c_rec,info_rec]    = data.loewner_null_rec(p_c,p_r,tab_red,'svd0');
%[c_rec,info_rec]    = data.loewner_null_recH(p_c,p_r,H,'svd0','optim');
info_rec
toc

[Var,Lag,Bary]  = data.decoupling(p_c,info_rec);

```

# Claim

This deposit consitutes a research code and is not aimed at 


