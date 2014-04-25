function dndt = deterministic(t,n,beta_c,delta_c,alpha)
   % calculate instantaneous rate of change for uninfected
   % hepatocytes (Tc), HCV virons (Vc), activated CD4+ cells (Th)
   % and HIV virons (Vh)

   % initialize population vector
   dndt = zeros(size(n));

   % fixed model parameters
   Sc=3661.7;          % hepatocyte recruitment
   dc=1.4e-3;          % uninfected hepatocyte death rate
   
   % COMMENTED OUT
   %beta_c=2.25e-7;      % rate of HCV infection
   

   r1=0.31;            % uninfected hepatocyte proliferation rate
   r2=4.4;             % infected hepatocyte proliferation rate  
   
   % COMMENTED OUT
   %delta_c=0.26;       % infected hepatocyte death rate   
   

   Tc_max=6.33e6;      % total maximum hepatocyte count
   c=11.5;             % rate of virion clearance
   p=15;               % rate of virion production 
   
   % COMMENTED OUT
   %alpha=0.003;       % dependence of virion clearance on CD4+
   
   Sh=9;               % CD4+ recruitment
   dh=9e-3;            % CD4+ death rate
   gamma=1e-8;         % dependence of CD4+ recruitment on viral load
   beta_h=4.1e-6;      % rate of HIV infection 
   delta_h=0.6;        % HIV virion death rate
   k=75;               % HIV burst size 
   e=0.1;              % HIV viral clearance

   
   % current population sizes
   Tc = n(1); % uninfected hepatocytes
   Vc = n(2); % HCV virons
   Th = n(3); % activated CD4+ cells
   Vh = n(4); % HIV virons

   % dTc/dt
   dndt(1) = Sc+r1*Tc*(1-((Tc+(c/p*Vc))/Tc_max))-dc*Tc-beta_c*Tc*Vc;
   
   % dVc/dt
   dndt(2) = (p/c)*beta_c*Tc*Vc+r2*Vc*(1-((Tc+(c/p*Vc))/Tc_max))-delta_c*(1+alpha*Th)*Vc;

   % dTh/dt
   dndt(3) = Sh*(1+gamma*Vc)-dh*Th-beta_h*Th*Vh;

   % dVh/dt
   dndt(4) = k/e*beta_h*Th*Vh-delta_h*Vh;

   % return dndt

end
