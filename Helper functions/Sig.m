function z = Sig(x,prm)
% returns Sigma(x) = prm(1) + prm(2) * exp(prm(3) * x) for x > 0; 1
% otherwise

I = x <= 0;
z = nan(size(x));

z(I) = 1;
z(~I) = sqrt(prm(1) + prm(2) * exp(prm(3) * x(~I)));


end