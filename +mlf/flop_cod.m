function nflop = flop_cod(nbr_pts)

nflop       = 0;
kappa(1)    = 1;
for ii = 1:numel(nbr_pts)
    nflop       = nflop + kappa(ii)*(nbr_pts(ii))^3;
    kappa(ii+1) = kappa(ii)*(nbr_pts(ii)); 
end