% define number of simulations
numSims = 10;

% define time-span
tspan = 1:1:365*40;

% define cancer threshold
threshold = 10e6;

% vector for storing simulation cancer results
hasCancer = zeros(numSims);

% vector for storing day cancer appeared
dayOfCancer = zeros(numSims);

% vector for storing beta_c
beta_c = zeros(numSims);

% vector for storing delta_c
delta_c = zeros(numSims);

% vector for storing alpha;
alpha = zeros(numSims);

% generate hypercube: numSims by 3 parameters
probs = lhsdesign(numSims, 3);

%%%%%%%%%%%%%%% Start simulations
for sim = 1:numSims
	     
  % draw random number between min and max for simulation
  %beta_c(sim) = unifrnd(1e-8, 1e-6);
  beta_c(sim) = logninv(probs(sim,1), 2.25e-7, 4.5e-7);

  % draw random number between min and max for simulation
  %delta(sim) = unifrnd(1e-3, 1);
  delta_c(sim) = logninv(probs(sim,2), 0.26, 0.16);
  
  alpha(sim) = logninv(probs(sim,3), 0.04, 0.02);

  % define initial conditions
  Tc_init = 5.31e6;
  Vc_init = 1.0e6;
  Th_init = 1000;
  Vh_init = 300;
  init = [Tc_init, Vc_init, Th_init, Vh_init];

  % run deterministic model using inputs for this simulation
  [T, N] = ode45(@deterministic, tspan, init, [], beta_c(sim), delta_c(sim), alpha(sim));

  % extract variables
  Tc = N(:,1);
  Vc = N(:,2);
  Th = N(:,3);
  Vh = N(:,4);


  %% plot results from deterministic model
  %figure(1)
  %clf
  %semilogy(T,Tc,'b-','LineWidth',2);
  %hold on;
  %semilogy(T,Vc,'r-','LineWidth',2);
  %hold on;
  %semilogy(T,Th,'g-','LineWidth',2);
  %hold on;
  %semilogy(T,Vh,'y-','LineWidth',2);
  %xlabel('Time (days)')
  %ylabel('RNA & cells/ml')
  %legend('Tc','Vc','Th','Vh')
  %title('Deterministic Model','Fontsize',12)
  %ylim([0 100000000])

  % call stochastic function
  HCC_cells = stochastic(T, N, beta_c(sim), delta_c(sim), alpha(sim));

  % above threshold vector of 0 (FALSE) and 1 (TRUE)
  aboveThreshold = HCC_cells > threshold;

  % check if patient got cancer
  if sum(aboveThreshold) > 1
    % Has CANCER, so store results and figure out
    % what day cancer appeared
    hasCancer(sim) = 1; % set hasCancer to TRUE for this simulation
    
    for t = 1:length(aboveThreshold)
      if aboveThreshold(t) == 1
	dayOfCancer(sim) = t; % store first day cancer appeared for this simulation
	break; % break out of for loop bc we only care about first day cancer appeared
      end 
    end
  end
end

% hasCancer has record of who got cancer
hasCancer

% dayOfCancer has record of what day they got cancer on
dayOfCancer

% beta_c has record of the value of param1 for each person
% delta has a record of the value of param2 for each person

