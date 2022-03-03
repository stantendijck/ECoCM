function x = Weibull_iCDF(u, lambda, k)
%Compute the inverse of the cumulative distribution function of a Weibull
%random variable with parameters (lambda, k) -> default = (lambda = 1, k = 1)

if nargin < 2
    lambda = 1;
end
if nargin < 3
    lambda = 1;
end

x = lambda * (-log(1-u)).^(1/k);

end