function z = Mu(x,prm)
% returns mu(x) = prm(1) + prm(2) * x ^ (prm(3)) for x >= 0, zero otherwise

I = x <= 0;
z = nan(size(x));

z(I) = 0;
z(~I) = prm(1) + prm(2) * x(~I).^(prm(3));


end