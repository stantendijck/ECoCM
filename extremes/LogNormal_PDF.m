function f = LogNormal_PDF(x, mu, sig)

if sig < 0
    f = nan;
    return
end

f = zeros(size(x));
f(x>0) = 1./(sqrt(2*pi) * sig * x) .* exp( -(log(x) - mu).^2/(2*sig^2) );

end
