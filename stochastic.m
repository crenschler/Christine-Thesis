function y = stochastic(T, N)
  % Runs the stochastic model for HCC development beginning with
  % the start of HCV infection and lasting 40 years using the time
  % courses from the deterministic model (deterministic.m)

  % T is the time vector produced by the ode solver of the deterministic
  % equation

  % N is the population vector produced by the ode solver of the
  % deterministic equation

  % set time span
  days = 365*50; % fifty years

  % y will store HCC cells over time
  y = zeros(days+1);

  % set parameters for stochastic model
  b = ;         % birth rate of HCC cells
  p = ;         % death rate of HCC cells
  alpha = ;     % mutation rate of hepatocytes
  delta = ;     % 

  % loop over time and track population of HCC cells
  for i = 1:(365*50)
    % get I pop for time i
    I = SOMETHING;

    % get E pop for time i
    E = SOMETHING;

    % Birth
    if (b*y(i)) < 5
      dyb = poissrnd(b*y(i));
    else
      dyb = b*y(i) + (normrnd(0,1)*sqrt(b*y(i)));

    % Death
    if SOMETHING2 < 5
      dyd = ;
    else
      dyd = ;

    % Mutation
    if SOMETHING3 < 5
      dym = ;
    else
      dym = ;

    % Update Total
    y(i+1) = max(y(i) + dyb - dyd + dym, 0);

  end

  % return y

end
