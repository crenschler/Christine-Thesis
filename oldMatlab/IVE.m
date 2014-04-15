function dndt = IVE(t,n) 
% calculate instantaneous rate of change for uninfected hepatocytes (Tc),
% virions(Vc), and activated CD4+ (Th)

dndt=zeros(size(n));

% model parameters
k=7.09e-7;      %infection rate per virion per day
delta=0.2279;     %death rate per effector cell
N=2.9;          %virions produced per day per infected cell
c=6;            %virion cleraance rate per day
alpha=9e-7;     %source term for effector cells per day
beta=6e-8;      %epansion rate of effector cells per infected cell per day
mew=0.1;        %death rate of effector cells per day
Hmax=2.4e6;     %maximum number of hepatocytes scaled to per ml of serum

H = n(1);
I = n(2);           % uninfected hepatocyes
V = n(3);           % virions/infected hepatocytes
E = n(4);           % activated CD4+ cells

dndt(1) = Hmax-I;
dndt(2) = (k*H*V)-(delta*E*I);
dndt(3) = (N*I)-(c*V);
dndt(4) = alpha+(beta*E*I)-(mew*E);

end