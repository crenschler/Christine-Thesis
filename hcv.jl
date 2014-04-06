# load libraries
using ODE
using DataFrames
using Gadfly

# define the model
function TcVcTh(t,x,param) # calculate instantaneous rate of change for
                           # uninfected hepatocytes (Tc),virions(Vc),
                           # and activated CD4+ (Th)

    Tc = x[1]           # uninfected hepatocyes
    Vc = x[2]           # virions/infected hepatocytes
    Th = x[3]           # activated CD4+ cells

    Sc = param[1]       # hepatocyte recruitment
    dc = param[2]       # hepatocyte death rate
    beta_c = param[3]   # rate of infection
    r1 = param[4]       # uninfected hepatocyte proliferation rate
    r2 = param[5]       # infected hepatocyte proliferation rate
    c = param[6]        # rate of virion clearance
    p = param[7]        # rate of virion production
    Tc_max = param[8]   # total maximum hepatocyte count
    delta = param[9]    # virion death rate
    alpha = param[10]   # dependence of clearance on CD4+
    Sh = param[11]      # CD4+ recruitment
    dh = param[12]      # CD4+ death rate
    gamma = param[13]   # dependence of recruitment on viral load

    Tc2 = Sc+r1*Tc*(1-((Tc+(c/p*Vc))/Tc_max))-dc*Tc-beta_c*Tc*Vc
    Vc2 = (p/c)*beta_c*Tc*Vc+r2*Vc*(1-((Tc+(c/p*Vc))/Tc_max))-delta*(1+alpha*Th)*Vc
    Th2 = Sh*(1+gamma*Vc)-dh*Th

    [Tc2,Vc2,Th2]

end

# initialize model
t = linspace (0,10000,1);
inits = [2.4e6,1.0e6,0.3068];
param = [3661.7, 1.4e-3, 7.4e-8, 0.31, 4.4, 11.5, 40.5, 6.33e6, 2.7266, 2.5e-5, 9e-7, 0.1, 0.25];
      # [s, dc, beta_c, r1, r2, c, p, Tc_max, delta, alpha, Sh, dh, gamma]

# run model
result = ode45((t,x)->TcVcTh(t,x,param),t,inits);

# collate results in DataFrame
df=DataFrame();
df["t"]=result[1];
df["Tc"]=result[2][:,1];
df["Vc"]=result[2][:,2];
df["Th"]=result[2][:,3];

# Plot using Gadfly
p = Gadfly.plot(df,x="t",y="Vc",Geom.line)


