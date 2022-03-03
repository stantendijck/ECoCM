addpath 'evd'
addpath 'standard functions'
addpath 'extremes'
addpath 'Helper Functions'

%% Generate Figure 1 from the Paper

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.15 0.05], [0.11 0.03]);
if ~make_it_tight,  clear subplot;  end

muprm = [1.134, 0.892, 0.225];
sigprm = [0.005,0.120,-0.455];
densprm = [0.893, 0.573, 3.803, 2.908, 1.550];
dist_ind = 2;

xl = linspace(0.01,100,1000);
yl = [10, 20, 30, 50, 100, 200];
col = jet(length(yl));
leg = cell(length(yl),1);
figure(1); clf;
subplot(1,1,1)
for iy = 1:length(yl)
    z = gaussian_survival_approx((log(yl(iy))- Mu(xl,muprm))./Sig(xl,sigprm).^(2/dist_ind)) + log_density(xl,densprm,dist_ind);
    plot(xl,z,'Color',col(iy,:)); hold on;
    leg{iy} = sprintf('y=%d',yl(iy));
end
set(gca,'YScale','log');
legend(leg);
xlabel('$x$','Interpreter','LaTeX'); ylabel('$\log g_y(x)$','Interpreter','LaTeX');
xlim([-3 90]);
savePics('figures/logDensity.pdf',1,'paper',5,5*1/2)


%% Generate Figure 2 from the Paper

make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.1 0.05], [0.1 0.05]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.18 0.1], [0.05 0.03]);
if ~make_it_tight,  clear subplot;  end

nl = 100;
alpha = linspace(0.01,0.99,nl); beta = linspace(0.01,0.99,nl)';
delta = (1-beta).^(-1); gamma = [1, 1.5, 2, 5];
figure(2); clf; 
for igamma = 1:length(gamma)
    subplot(1,4,igamma);
    Z = gamma(igamma)*(1-alpha).^(delta-1).*(delta-1+alpha);
    Z(Z>1) = 1.0001;
    
    C = zeros(length(alpha),length(beta),3);
    for ialpha = 1:length(alpha)
        for ibeta = 1:length(beta)
            C(ialpha,ibeta,:) = [0,1,0];
        end
    end
    surf(alpha,beta,log(1+Z),C);
    shading interp
    xlabel('\alpha'); ylabel('\beta'); zlabel('FUN'); title(sprintf('\\gamma = %.1f',gamma(igamma)));
    C = zeros(2,2,3);
    for ialpha = 1:2
        for ibeta = 1:2
            C(ialpha,ibeta,:) = [1,0,0];
        end
    end
    hold on; surf([0.01,0.99],[0.01;0.99],[0.01,0.99].*[0.01;0.99].*0+log(2),C);
    view(2);
end
savePics('figures/cEqualToOne.pdf',1,'paper',10,10*1/5);


%% Generate Figure 3 from the Paper
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.15 0.1], [0.05 0.03]);
if ~make_it_tight,  clear subplot;  end


u = 100;
nl = 100;
alpha = linspace(0.01,0.99,nl); beta = linspace(0.01,0.99,nl);
gamma = [1, 1.5, 2, 5];
figure(3); clf; 
for igamma = 1:length(gamma)
    M = zeros(length(alpha),length(beta));
    eta = zeros(length(alpha),length(beta));
    for ialpha = 1:length(alpha)
        for ibeta = 1:length(beta)
            delta = (1-beta(ibeta))^(-1);
            h = @(x,u)(-gamma(igamma)*(u-alpha(ialpha)*x).^delta./x.^(beta(ibeta)*delta) - x);

            x0 = fminbnd(@(x)(-h(x,u)),u,u/alpha(ialpha));
            if abs(x0 - u) < 1e-4
                x0 = u;
            elseif abs(x0 - u/alpha(ialpha)) < 1e-4
                x0 = u/alpha(ialpha);
            end
            c = x0/u;

            M(ialpha,ibeta) = c;

            eta(ialpha,ibeta) = (gamma(igamma)*(1-alpha(ialpha)*c)^delta/c^(delta-1) + c)^(-1);
        end
    end

    subplot(1,4,igamma);
    surf(alpha,beta,real(eta)','edgecolor','none');
    xlabel('\alpha'); ylabel('\beta'); zlabel('\eta'); title(sprintf('\\gamma = %.1f',gamma(igamma)));
    zlim([0 1]);
end
savePics('figures/htModelPrmComb.pdf',1,'paper',10,10*1/5);


%% Generate Figure 4 from the Paper

% Initialisation
n =  1e4;
rng(12345);

xi = 0.35;
alpha = 0; beta = 1 - xi;
C = 1; gamma = xi; delta = 1/xi;

% Simulate from inverted logistic:
[A_laplace, ~, ~] = simulate_alg21(n,2,xi);
Data = -A_laplace;
writematrix(Data, 'simulatedData.txt')

% True eta
true_eta = 2^(-xi);
true_chi_bar = 2*true_eta - 1;

% HT eta
HT_eta = 1/(1+gamma);
HT_chi_bar = 2*HT_eta - 1;

% Empirical eta
[pOrg,uOrg,pconfOrg] = chi_bar(Laplace_CDF(Data));

%% Numerical integration of eta(p) of the HT model

H_overline = @(x)(C * exp(-gamma*x.^(delta)));
integrand = @(x,u)(H_overline(u.*x.^(-beta)) .* exp(-x));
NumInt = @(u)(integral(@(x)(integrand(x+u,u)),0,inf));

% const = C * sqrt(2 * pi)/(2 * sqrt(gamma*delta*(delta - 1)));
% exp_factor = @(u)(-(1+gamma)*u + (gamma * (delta - 1) - 1)/(2*gamma*delta*(delta - 1)) * 1/u);
% est = @(u)(const * u^(-1/2) * exp(exp_factor(u)));

% Using the log method:
logH_overline = @(x)(log(C) - gamma*x.^(delta));
log_integrand = @(x,u)(log(C) - gamma * u^(delta) * x.^(-beta*delta) - x);
log_NumInt = @(u)(log_integrand(u,u) + log(integral(@(x)(exp(log_integrand(x+u,u) - log_integrand(u,u))),0,inf)));

B = 1e3;
ul = sort([linspace(1,300,B-100)';log(1./linspace(0.1,0.5)')]);
p1 = zeros(B,1);
% p2 = zeros(B,1);
log_p1 = zeros(B,1);
% log_p2 = zeros(B,1);
for iu = 1:B
    u = ul(iu);
    p1(iu) = NumInt(u);    
    log_p1(iu) = log_NumInt(u);
end

%% Plot Figure 4
make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.1 0.05], [0.1 0.05]);
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.1 0.05], [0.06 0.03]);
if ~make_it_tight,  clear subplot;  end

figure(4);
clf;
subplot(1,2,1);
plot(uOrg,(pOrg+1)/2,'b-');
hold on;
% plot(uOrg,(pconfOrg(:,1)+1)/2,'b--');
% plot(uOrg,(pconfOrg(:,2)+1)/2,'b--','HandleVisibility','off');
plot(1-1/2*exp(-ul),-ul./log_p1,'k-');
plot(xlim,true_eta*ones(2,1),'r-');
plot(xlim,HT_eta*ones(2,1),'g-');
xlim([0.7 1.01]); ylim([0.7 0.9]);
patch([uOrg fliplr(uOrg)],([pconfOrg(:,1)' fliplr(pconfOrg(:,2)')] + 1)/2,'b','FaceAlpha',0.2,'EdgeAlpha',0);
xlabel('$p$','Interpreter','latex'); ylabel('$\eta$','Interpreter','latex')
legend({'$\hat{\eta}(p)$','$\eta_{HT}(p)$','$\eta$','$\eta_{HT}$'},'Interpreter','latex');

subplot(1,2,2);
plot(ul,-ul./log_p1,'k-'); hold on;
plot(xlim,true_eta*ones(2,1),'r-');
plot(xlim,HT_eta*ones(2,1),'g-');
% plot(1-1/2*exp(-ul),2*(-ul./log_p2) - 1,'k-.'); hold on;
xlabel('$-\log(2-2p)$','Interpreter','latex');
ylabel('$\eta$','Interpreter','latex');
ylim([0.7 0.9]);
legend({'$\eta_{HT}(p)$','$\eta$','$\eta_{HT}$'},'Interpreter','latex');
% savePics('figures/logisticEta.pdf',1,'paper',10,10*1/2)
% Note: eta = eta(p) for Inv-Logistic model


