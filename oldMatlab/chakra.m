% model parameters
k=1.05e-7;      %infection rate per virion per day
delta=0.1632;     %death rate per effector cell
N=2.9;          %virions produced per day per infected cell
c=6;            %virion cleraance rate per day
alpha=9e-7;     %source term for effector cells per day
beta=6e-8;      %epansion rate of effector cells per infected cell per day
mew=0.1;        %death rate of effector cells per day
Hmax=2.4e6;     %maximum number of hepatocytes scaled to per ml of serum

% define time
tspan = [0,500];

% define initial conditions
I_init= 10;
V_init= 10;
E_init= 9e-8;

% run model
[T_out, N_out]=ode45(@IVE, tspan, [I_init, V_init, E_init], []);

% rename variables
I = N_out(:,1);
V = N_out(:,2);
E = N_out(:,3);

% plot
figure(3)
clf
semilogy(T_out,I,'b-','LineWidth',2);
hold on;
semilogy(T_out,V,'r-','LineWidth',2);
hold on;
semilogy(T_out,E,'y-','LineWidth',2);
xlabel('Time (days)')
ylabel('HCV RNA & cells/ml')
legend('I','V','E')
title('ODE model','Fontsize',12)
ylim([0 100000000000])