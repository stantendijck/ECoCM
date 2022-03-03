function f = GPD_PDF(x, mu, sigma, xi)

x(x == mu) = mu + 1e-8;
f = 1/sigma * (1 + xi * (x - mu)/sigma) .^(-1/xi-1);
I = x < mu;
f(I) = 0;
if xi < 0
    I = x > mu - sigma/xi;
    f(I) = 0;
end
    

end