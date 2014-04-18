% define number of simulations
numSims = 1000;

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

% vector for storing delta
delta = zeros(numSims);

% generate hypercube: numSims by 2 parameters
probs = lhsdesign(numSims, 2);

%%%%%%%%%%%%%%% Start simulations
for sim = 1:numSims
	     
  % draw random number between min and max for simulation
  %beta_c(sim) = unifrnd(1e-8, 1e-6);
  beta_c(sim) = logninv(probs(sim,1), MEAN, STD)

  % draw random number between min and max for simulation
  %delta(sim) = unifrnd(1e-3, 1);
  delta(sim) = logninv(probs(sim,2), MEAN, STD)

  % define initial conditions
  Tc_init = 2.4e6;
  Vc_init = 1.0e6;
  Th_init = 0.3068;
  Vh_init = 300;
  init = [Tc_init, Vc_init, Th_init, Vh_init];

  % run deterministic model using inputs for this simulation
  [T, N] = ode45(@deterministic, tspan, init, [], beta_c(sim), delta(sim));

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
  HCC_cells = stochastic(T, N);

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
% dayOfCancer has record of what day they got cancer on
% beta_c has record of the value of param1 for each person
% delta has a record of the value of param2 for each person

