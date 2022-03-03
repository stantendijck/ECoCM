function [p,u,pconf] = chi_bar(x)
%We need x to be on standard uniform margins, then this function 
%calculates chi(u) = P(X_1 > u | X_2 > u) for different values of u

if nargout == 3
    flag_conf = true;
else
    flag_conf = false;
end
n = size(x,1);

poss_val = sort(min(x,[],2));
max_u = poss_val(end-1);
min_u = median(x(:,1));

u = linspace(min_u,max_u,1000);
    
p = zeros(length(u),1);
pconf = zeros(length(u),2);

for i = 1:length(u)
%     p(i) = mean(x(x(:,1)>u(i),2) > u(i));
%     p(i) = 2*log(1-u(i))./log(1-2*u(i) + mean(x(:,1)<u(i) & x(:,2)<u(i))) - 1;
%     p(i) = 2*log(mean(x(:,1)>u(i)))./log(mean(x(:,1)>u(i) & x(:,2)>u(i))) - 1;
    chi_i = mean(x(:,1)>u(i) & x(:,2)>u(i));
    p(i) = 2*log(1-u(i))./log(chi_i) - 1;
    if isnan(p(i))
        p(i) = 0;
    end
    if flag_conf
        if n*chi_i < 30
            binomconf = binoinv([0.025,0.975],n,chi_i);
        else
            binomconf = norminv([0.025,0.975],n*chi_i,sqrt(n*chi_i*(1-chi_i)));
        end
        pconf(i,:) = 2*log(1-u(i))./log(1/n*binomconf) - 1;
        if isnan(pconf(i,1))
            pconf(i,1) = 0;
        end
        if isnan(pconf(i,2))
            pconf(i,2) = 0;
        end
    end
end


end