% define time-span
tspan = [0,365*40]

% define initial conditions
Tc_init = 2.4e6;
Vc_init = 1.0e6;
Th_init = 0.3068;
Vh_init = 10;
init = [Tc_init, Vc_init, Th_init, Vh_init]

% run deterministic model
[T, N] = ode45(@deterministic, tspan, init, [])

% rename variables
Tc = N(:,1);
Vc = N(:,2);
Th = N(:,3);
Vh = N(:,4);

% plot results from deterministic model
figure(1)
clf
semilogy(T_out,Tc,'b-','LineWidth',2);
hold on;
semilogy(T_out,Vc,'r-','LineWidth',2);
hold on;
semilogy(T_out,Th,'g-','LineWidth',2);
hold on;
semilogy(T_out,Vh,'y-','LineWidth',2);
xlabel('Time (days)')
ylabel('RNA & cells/ml')
legend('Tc','Vc','Th','Vh')
title('Deterministic Model','Fontsize',12)
ylim([0 100000000])

% call stochastic function
HCC_cells = stochastic(T, N)

% interpret results 

% check if cancer evolved
threshold = 10e6; % threshold HCC cell levels to be considered cancer

% if greater than threshold, then that individual developed cancer
if HCC_cells > threshold
    
