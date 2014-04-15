function dndt = IVE(t,n) 
% calculate instantaneous rate of change for uninfected hepatocytes (Tc),
% virions(Vc), and activated CD4+ (Th)

dndt=zeros(size(n));

% model parameters
k=2.25e-7;      %infection rate per virion per day
delta=0.26;     %death rate per effector cell
N=2.9;          %virions produced per day per infected cell
c=6;            %virion cleraance rate per day
alpha=9e-7;     %source term for effector cells per day
beta=6e-8;      %epansion rate of effector cells per infected cell per day
mew=0.1;        %death rate of effector cells per day
Hmax=2.4e6;     %maximum number of hepatocytes scaled to per ml of serum

I = n(1);           % uninfected hepatocyes
V = n(2);           % virions/infected hepatocytes
E = n(3);           % activated CD4+ cells
H = n(4);

dndt(1) = (k*H*V)-(delta*E*I);
dndt(2) = (N*I)-(c*V);
dndt(3) = alpha+(beta*E*I)-(mew*E);
dndt(4) = Hmax-I;

end