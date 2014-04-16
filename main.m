%P define number of simulations
numSims = 1000;

% define time-span
tspan = 1:1:365*40;

% define cancer threshold
threshold = 10e6;

% vector for storing simulation cancer results
hasCancer = zeros(numSims)

% vector for storing day cancer appeared
dayOfCancer = zeros(numSims)

% vector for storing BETA_C
BETA_C = zeros(numSims)

% vector for storing DELTA
DELTA = zeros(numSims)

%%%%%%%%%%%%%%% START SIMULATIONS
for sim = 1:numSims
  % sim is the counter that keeps track of the current simulation
	     
  % draw random number between min and max for this simulation
  BETA_C(sim) = unifrnd(ENTER MIN, ENTER MAX);

  % draw random number between min and max for this simulation
  DELTA(sim) = unifrnd(ENTER MIN, ENTER MAX);

  % define initial conditions
  Tc_init = 2.4e6;
  Vc_init = 1.0e6;
  Th_init = 0.3068;
  Vh_init = 300;
  init = [Tc_init, Vc_init, Th_init, Vh_init];

  % run deterministic model using inputs for this simulation
  [T, N] = ode45(@deterministic, tspan, init, [], BETA_C(sim), DELTA(sim));

  % extract variables
  Tc = N(:,1);
  Vc = N(:,2);
  Th = N(:,3);
  Vh = N(:,4);


  %%%% COMMENT OUT PLOTTING FOR NOW
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
  %%% CHECK IF THIS WORKS, might need to change the notation
  %%% if this doesnt produce a vector but i think it does
  aboveThreshold = HCC_cells > threshold

  % check if patient got cancer
  if sum(aboveThreshold) > 1
    % Has CANCER, so store results and figure out
    % what day cancer appeared
    hasCancer(sim) = 1 % set hasCancer to TRUE for this simulation
    
    for t = 1:length(aboveThreshold)
      if aboveThreshold(t) == 1
	dayOfCancer(sim) = t; % store first day cancer appeared for this simulation
	break; % break out of for loop bc we only care about first day cancer appeared
      end 
    end
  end
end
%%%%%%%%%%%%%%%%% END SIMULATIONS

% hasCancer should have record of who got cancer
% dayOfCancer should have a record of what day they got cancer on
% BETA_C has record of the value of param1 for each person
% DELTA has a record of the value of param2 for each person    
