 tspan = 1:1:500;
  
% define initial conditions
  Tc_init = 5.31e6;
  Vc_init = 1.0e6;
  Th_init = 1000;
  Vh_init = 0.3;
  init = [Tc_init, Vc_init, Th_init, Vh_init];

  % run deterministic model using inputs for this simulation
  [T, N] = ode45(@deterministic_p, tspan, init);

  % extract variables
  Tc = N(:,1);
  Vc = N(:,2);
  Th = N(:,3);
  Vh = N(:,4);


  %% plot results from deterministic model
  figure(1)
  clf
  semilogy(T,Tc,'b-','LineWidth',2);
  hold on;
  semilogy(T,Vc,'r-','LineWidth',2);
  hold on;
  semilogy(T,Th,'g-','LineWidth',2);
  hold on;
  semilogy(T,Vh,'y-','LineWidth',2);
  xlabel('Time (days)')
  ylabel('RNA & cells/ml')
  legend('Tc','Vc','Th','Vh')
  title('Deterministic Model','Fontsize',12)
  ylim([0 100000000])