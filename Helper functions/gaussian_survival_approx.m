function y = gaussian_survival_approx(z)

logmillsRatioBnd2 = @(z)(-1/2*log(2*pi) -z.^2/2 + log(z) - log(z.^2+1));

y = zeros(size(z));

I1 = z > 10;
y(I1) = logmillsRatioBnd2(z(I1));

I2 = z >= -5 & z <= 10;
y(I2) = log(normcdf(z(I2),'upper'));

I3 = z < -5;
y(I3) = -exp(logmillsRatioBnd2(-z(I3)));


end
