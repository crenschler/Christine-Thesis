% model parameters
Sc=3661.7;          % hepatocyte recruitment
dc=1.4e-3;          % hepatocyte death rate
beta_c=7.4e-8;      % rate of infection
r1=0.31;            % uninfected hepatocyte proliferation rate
r2=4.4;             % infected hepatocyte proliferation rate
c=11.5;             % rate of virion clearance
p=40.5;             % rate of virion production 
Tc_max=6.33e6;      % total maximum hepatocyte count
delta=2.7266;       % virion death rate
alpha=2.5e-5;       % dependence of clearance on CD4+
Sh=9e-7;            % CD4+ recruitment
dh=0.1;             % CD4+ death rate
gamma=0.10;         % dependence of recruitment on viral load

% define time
tspan = [0,365*30];

% define initial conditions
Tc_init= 2.4e6;
Vc_init= 1.0e6;
Th_init= 0.3068;

% run model
[T_out, N_out]=ode45(@TcVcTh, tspan, [Tc_init, Vc_init, Th_init], []);

% rename variables
Tc = N_out(:,1);
Vc = N_out(:,2);
Th = N_out(:,3);

csvwrite('T_out.csv', T_out)
csvwrite('N_out.csv', N_out)



% plot
figure(1)
clf
semilogy(T_out,Tc,'b-','LineWidth',2);
hold on;
semilogy(T_out,Vc,'r-','LineWidth',2);
hold on;
semilogy(T_out,Th,'g-','LineWidth',2);
xlabel('Time (days)')
ylabel('HCV RNA & cells/ml')
legend('Tc','Vc','Th')
title('ODE model','Fontsize',12)
ylim([0 10000000])
