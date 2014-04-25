function dndt = TcVcTh(t,n) 
% calculate instantaneous rate of change for uninfected hepatocytes (Tc),
% virions(Vc), and activated CD4+ (Th)

dndt=zeros(size(n));

% model parameters
   Sc=3661.7;          % hepatocyte recruitment
   dc=1.4e-3;          % uninfected hepatocyte death rate
   beta_c=2.25e-7;     % rate of HCV infection
   r1=0.31;            % uninfected hepatocyte proliferation rate
   r2=4.4;             % infected hepatocyte proliferation rate  
   delta=0.26;         % infected hepatocyte death rate   
   Tc_max=6.33e6;      % total maximum hepatocyte count
   c=11.5;             % rate of virion clearance
   p=15;              % rate of virion production 
   alpha=0.003;     % dependence of virion clearance on CD4+
   Sh=9;               % CD4+ recruitment
   dh=.009;            % CD4+ death rate
   gamma=1e-8;         % dependence of CD4+ recruitment on viral load


Tc = n(1);           % uninfected hepatocyes
Vc = n(2);           % virions/infected hepatocytes
Th = n(3);           % activated CD4+ cells

dndt(1) = Sc+r1*Tc*(1-((Tc+(c/p*Vc))/Tc_max))-dc*Tc-beta_c*Tc*Vc;
dndt(2) = (p/c)*beta_c*Tc*Vc+r2*Vc*(1-((Tc+(c/p*Vc))/Tc_max))-delta*(1+alpha*Th)*Vc;
dndt(3) = Sh*(1+gamma*Vc)-dh*Th;

end
