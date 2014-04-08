function dndt = TcVcThVh(t,n) 
% calculate instantaneous rate of change for uninfected hepatocytes (Tc),
% virions(Vc), and activated CD4+ (Th)

dndt=zeros(size(n));

% model parameters
Sc=3661.7;          % hepatocyte recruitment
dc=1.4e-3;          % hepatocyte death rate
beta_c=7.4e-8;      % rate of HCV infection
r1=0.31;            % uninfected hepatocyte proliferation rate
r2=4.4;             % infected hepatocyte proliferation rate
c=11.5;             % rate of virion clearance
p=40.5;             % rate of virion production 
Tc_max=6.33e6;      % total maximum hepatocyte count
delta=2.7266;       % HCV virion death rate
alpha=2.5e-5;       % dependence of clearance on CD4+
Sh=9e-7;            % CD4+ recruitment
dh=0.1;             % CD4+ death rate
gamma=0.10;         % dependence of recruitment on viral load
beta_h=4.1e-6;      % rate of HIV infection 
delta_h=0.3;        % HIV virion death rate
k=75;               % HIV burst size 
e=0.6;              % HIV viral clearance

Tc = n(1);          % uninfected hepatocyes
Vc = n(2);          % HCV free virus/infected hepatocytes
Th = n(3);          % activated CD4+ cells
Vh = n(4);          % HIV free virus

dndt(1) = Sc+r1*Tc*(1-((Tc+(c/p*Vc))/Tc_max))-dc*Tc-beta_c*Tc*Vc;
dndt(2) = (p/c)*beta_c*Tc*Vc+r2*Vc*(1-((Tc+(c/p*Vc))/Tc_max))-delta*(1+alpha*Th)*Vc;
dndt(3) = Sh*(1+gamma*Vc)-dh*Th-beta_h*Th*Vh;
dndt(4) = (k/e)*beta_h*Th*Vh-delta_h*Vh;

end