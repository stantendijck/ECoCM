function F = LogNormal_iCDF(p, mu, sig)

if sig < 0
    F = nan;
    return
end

F = exp(mu + norminv(p) * sig);

end
