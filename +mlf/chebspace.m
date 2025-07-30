% Chebychev point repartition
% Author: C. Poussot-Vassal [MOR Digital Systems & ONERA]
%

function x = chebspace(a,b,n)

x = zeros(1,n);
for k = 0:n-1
    x(1,k+1) = cos((k+.5)*pi/n);
end

x = .5*(a+b)+.5*(b-a)*x;
x = fliplr(x);