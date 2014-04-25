% define time
tspan = [0,500];

% define initial conditions
Tc_init= 5.31e6;
Vc_init= 1000000;
Th_init= 1000;

% run model
[T_out, N_out]=ode45(@TcVcTh, tspan, [Tc_init, Vc_init, Th_init], []);

% rename variables
Tc = N_out(:,1);
Vc = N_out(:,2);
Th = N_out(:,3);

csvwrite('T_out.csv', T_out)
csvwrite('N_out.csv', N_out)


% plot
figure(2)
clf
semilogy(T_out,Tc,'b-','LineWidth',2);
hold on;
semilogy(T_out,Vc,'r-','LineWidth',2);
hold on;
semilogy(T_out,Th,'g-','LineWidth',2);
xlabel('Time (days)')
ylabel('HCV RNA & cells/ml')
legend('Tc','Vc','Th')
title('Monoinfection ODE','Fontsize',12)
ylim([0 10000000])
