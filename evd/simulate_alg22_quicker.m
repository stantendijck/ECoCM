function [x_laplace, x_gumbel, x_frechet] = simulate_alg22_quicker(n, d, alpha, theta)
% simulate from asymmetric logistic model using Algorithm 2.2 in Stephenson
% (2003)
%
% INPUT:
% n         number of samples
% d         dimension
% alpha     (2^d - d - 1) x 1 dependence parameters
% theta     d*2^(d-1) x 1 dependence parameters
%
% OUTPUT:
% x         (n x d) sample on Gumbel scale
% y         (n x d) sample on Frechet scale

% Check input
if length(alpha) ~= 2^d - d - 1
    error('In simulate_alg22, length(alpha) needs to be equal to 2^d - d - 1');
end
if length(theta) ~= d*2^(d-1)
    error('In simulate_alg22, length(theta) needs to be equal to d*2^(d-1)');
end
if any(alpha) <=0  || any(alpha > 1)
    error('In simulate_alg22, 0 < alpha(i) <= 1 for all i');
end

B = cell(1,d);
for iB = 1:d
    B{iB} = combnk(1:d,iB);
end

for i = 1:d
    s = 0;
    counter = 1;
    for j = 1:d
        for k = 1:size(B{j},1)
            for m = 1:size(B{j},2)
                if B{j}(k,m) == i
                    s = s + theta(counter);
                end
                counter = counter + 1;
            end
        end
    end
    if s ~= 1
        error('In simulate_alg22, sum of theta_%d over all sets b does not equal 1',i);
    end
end

flag = false;
counter_alpha = 1;
counter_theta = d + 1;
for i = 2:d
    for j = 1:size(B{i},1)
        if alpha(counter_alpha) == 1
            for k = 1:size(B{i},2)
                if theta(counter_theta) ~= 0
                    flag = true;
                end
                counter_theta = counter_theta + 1;
            end
        end
        counter_alpha = counter_alpha + 1;
    end
end
if flag
    error('In simulate_alg22, if alpha(b) = 1, then theta(i,b) = 0 for all i must hold');
end

flag = false;
counter_theta = d + 1;
for i = 2:d
    for j = 1:size(B{i},1)
        s = 0;
        for k = 1:size(B{i},2)
            s = s + (theta(counter_theta)==0);
            counter_theta = counter_theta + 1;
        end
        if s == 1
            flag = true;
        end
    end
end
if flag
    error('In simulate_alg22, if theta(j,b)=0 for all j in b_{-i}, then theta(i,b) needs to be 0');
end

w = cell(1,d);

for i = 1:d
    sz_B = size(B{i}); sz_B(1) = sz_B(1)*n;
    w{i} = -log(rand(sz_B));
end

x_b = cell(1,d);

sb = simulate_PS(n,[ones(1,d),alpha]); 

counter_alpha = 1;
counter_theta = 1;
counter_sb = 1;
t_alpha = 1;
for i = 1:d
    sz_B = size(B{i}); sz_B(1) = sz_B(1)*n;
    x_b{i} = zeros(sz_B);
    for j = 1:size(B{i},1)
        if i >= 2
            t_alpha = alpha(counter_alpha);
            counter_alpha = counter_alpha + 1;
        end
        for k = 1:size(B{i},2)
            x_b{i}(j:size(B{i},1):end,k) = theta(counter_theta)*(sb(:,counter_sb)./w{i}(j:size(B{i},1):end,k)).^(t_alpha);
            counter_theta = counter_theta + 1;
        end
        counter_sb = counter_sb + 1;
    end
end

x_frechet = zeros(n,d) - Inf;
for i = 1:d
    for j = 1:d
        for k = 1:size(B{j},1)
            for m = 1:size(B{j},2)
                if B{j}(k,m) == i
                    x_frechet(:,i) = max(x_frechet(:,i),x_b{j}(k:size(B{j},1):end,m));
                end
            end
        end
    end
end

x_gumbel = log(x_frechet);
x_unif = exp(-x_frechet.^(-1));
x_laplace = -sign(x_unif - 1/2) .* log(1 - 2*abs(x_unif-1/2));


end