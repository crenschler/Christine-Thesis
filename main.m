% define time-span
tspan = [0,365*40]

% define initial conditions
Tc = 2.4e6;
Vc = 1.0e6;
Th = 0.3068;
Vh = 10;
init = [Tc, Vc, Th, Vh]

% run deterministic model
[T, N] = ode45(@deterministic, tspan, initial_conditions, [])

% plot results from deterministic model


% call stochastic function
HCC_cells = stochastic(T, N)

%%%% interpret results ?

% check if cancer evolved
threshold = ; % threshold HCC cell levels to be considered cancer

% if the sum below is greater than 1, then that person got cancer
sum(HCC_cells > threshold)
