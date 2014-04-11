clc
syms Tc Vc Th Vh Sc r1 c p Tc_max dc beta_c r2 delta alpha Sh gamma dh beta_h delta_h k e 'real'
e1 = eval('Sc+r1*Tc*(1-((Tc+(c/p*Vc))/Tc_max))-dc*Tc-beta_c*Tc*Vc');
e2 = eval('(p/c)*beta_c*Tc*Vc+r2*Vc*(1-((Tc+(c/p*Vc))/Tc_max))-delta*(1+alpha*Th)*Vc');
e3 = eval('Sh*(1+gamma*Vc)-dh*Th-beta_h*Th*Vh');
e4 = eval('(k/e)*beta_h*Th*Vh-delta_h*Vh');
sol=solve(e1,Tc,e2,Vc,e3,Th,e4,Vh);
sol.Tc=simplify(sol.Tc);
sol.Vc=simplify(sol.Vc);
sol.Th=simplify(sol.Th);
sol.Vh=simplify(sol.Vh);
disp([sol.Tc sol.Vc sol.Th sol.Vh])

