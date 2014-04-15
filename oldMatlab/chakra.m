% model parameters
k=7.09e-7;      %infection rate per virion per day
delta=0.2279;     %death rate per effector cell
N=2.9;          %virions produced per day per infected cell
c=6;            %virion cleraance rate per day
alpha=9e-7;     %source term for effector cells per day
beta=6e-8;      %epansion rate of effector cells per infected cell per day
mew=0.1;        %death rate of effector cells per day
Hmax=2.4e6;     %maximum number of hepatocytes scaled to per ml of serum

% define time
tspan = [0,500];

% define initial conditions
H_init= 2.4e6;
I_init= 1;
V_init= 1;
E_init= 9e-8;

% run model
[T_out, N_out]=ode45(@IVE, tspan, [H_init, I_init, V_init, E_init], []);

% rename variables
H = N_out(:,1);
I = N_out(:,2);
V = N_out(:,3);
E = N_out(:,4);

% plot
figure(3)
clf
semilogy(T_out,H,'b-','LineWidth',2);
hold on;
semilogy(T_out,I,'r-','LineWidth',2);
hold on;
semilogy(T_out,V,'g-','LineWidth',2);
hold on;
semilogy(T_out,E,'y-','LineWidth',2);
xlabel('Time (days)')
ylabel('HCV RNA & cells/ml')
legend('H','I','V','E')
title('ODE model','Fontsize',12)
ylim([0 10000000])