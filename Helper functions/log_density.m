function z = log_density(x,prm,distribution_index)
%distribution_index: 
% 1     Weibull(prm(1),prm(2))
% 2     logNormal(prm(1),prm(2)) for x < prm(3); Weibull(prm(4),prm(5)) else
%output: log density

switch distribution_index
    case 1
        lambda = prm(1);
        k = prm(2);
        
        I = (x < 0 | (x == 0 & k > 1));
        
        z = nan(size(x));
        z(I) = -1e5;
        z(~I) = log(k) - k*log(lambda) + (k-1) * log(x(~I)) -(x(~I)/lambda).^k;
    case 2
        logMu = prm(1);
        logSigma = prm(2);
        thr = prm(3);
        lambda = prm(4);
        k = prm(5);

        I = (x < 0 | (x == 0 & k > 1));
        Inonexc = ~I & (x <= thr);
        Iexc = ~I & (x > thr);

        z = nan(size(x));
        z(I) = -1e5;
        z(Inonexc) = -1/2*log(2*pi) - log(logSigma) - log(x(Inonexc)) - (log(x(Inonexc)) - logMu).^2/(2*logSigma^2);
        z(Iexc) = log(k) - k*log(lambda) + (k-1) * log(x(Iexc)) -(x(Iexc)/lambda).^k;
    otherwise
        z = nan(size(x));
end



end