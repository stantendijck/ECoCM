function NLOGL = LogNormal_like(x, mu, sig)

if sig < 0 || any(x < 0)
    NLOGL = inf;
    return
end

n = length(x);
NLOGL = n/2 * log(2*pi) + n*log(sig) + sum(log(x)) + sum((log(x) - mu).^2)/(2*sig^2);

end