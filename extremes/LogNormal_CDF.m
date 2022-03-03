function F = LogNormal_CDF(x, mu, sig)

if sig < 0
    F = nan;
    return
end

F = zeros(size(x));
F(x>0) = normcdf((log(x) - mu)/sig);

end
