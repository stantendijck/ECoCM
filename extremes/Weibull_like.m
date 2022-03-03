function NLOGL = Weibull_like(x, lambda, k)

if lambda < 0 || k < 0 || any(x < 0)
    NLOGL = inf;
    return
end

n = length(x);
NLOGL = n*k*log(lambda) - n*log(k) - (k-1) * sum(log(x)) + sum((x/lambda).^k);

end