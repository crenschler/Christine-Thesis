function y = stochastic(T, N)
  % Runs the stochastic model for HCC development beginning with
  % the start of HCV infection and lasting 40 years using the time
  % courses from the deterministic model (deterministic.m)

  % T is the time vector produced by the ode solver of the deterministic
  % equation

  % N is the population vector produced by the ode solver of the
  % deterministic equation

  % set time span
  days = 365*40; % forty years

  % y will store HCC cells over time
  y = zeros(days+1);

  % set parameters for stochastic model
  b = 0.26;         % birth rate of HCC cells
  p = 0.01;         % death rate of HCC cells
  alpha = 0.4e-10;  % mutation rate of hepatocytes
  delta = 0.26;     % 

  % loop over time and track population of HCC cells
  for i = 1:days
    % get Vc pop for time i
    Vc = N(i,2);

    % get Th pop for time i
    Th = N(i,3);

    % Birth
    if (b*y(i)) < 5
      dyb = poissrnd(b*y(i));
    else
      dyb = b*y(i) + (normrnd(0,1)*sqrt(b*y(i)));
    end
    
    % Death
    if (p*y(i)) < 5
      dyd = poissrnd(p*y(i));
    else
      dyd = p*y(i) + (normrnd(0,1)*sqrt(p*y(i)));
    end
    
    % Mutation
    if (alpha*delta*Th*Vc) < 5
      dym = poissrnd(alpha*delta*Th*Vc);
    else
      dym = alpha*delta*Th*Vc + (normrnd(0,1)*sqrt(alpha*delta*Th*Vc));
    end
    
    % Update Total
    y(i+1) = max(y(i) + dyb - dyd + dym, 0);

  end

  % return y

end
