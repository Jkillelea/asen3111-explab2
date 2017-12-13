function [e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N)
  % RETURNS:
  % e    -> span efficiency factor
  % c_L  -> coefficient of lift
  % c_Di -> coefficient of induced drag
  % PARAMS:
  % b      -> span (feet)
  % a0_t   -> cross sectional lift slope at tip
  % a0_r   -> cross sectional lift slope at root
  % c_t    -> chord line at tip (feet)
  % c_r    -> chord line at root (feet)
  % aero_t -> zero-lift angle of attack at tip
  % aero_r -> zero-lift angle of attack at root
  % geo_t  -> geometric angle of attack at tip
  % geo_r  -> geometric angle of attack at root
  % N      -> number of odd terms to include in the series expansion for circulation

  S  = (c_r + c_t)*b/2;
  AR = b^2/S;

  % lift slope and chord as functions of theta
  a0     = @(t) a0_r   + (a0_t - a0_r)*cos(t);     % lift slope
  c      = @(t) c_r    + (c_t - c_r)*cos(t);       % chord length
  alpha0 = @(t) aero_r + (aero_t - aero_r)*cos(t); % zero lift angle of attack
  alpha  = @(t) geo_r   + (geo_t - geo_r)*cos(t);  % geometric angle of attack
                                                   % (what we're flying at)
  % the theta range
  thetas = (1:N)*pi/(2*N); % half the wing

  % A1, A3, A5 ...
  nvals = 1:2*N;
  nvals = nvals(mod(nvals, 2) == 1); % filter even terms out

  % each An term is multiplied by a constant
  C = @(n, t) 4*b*sin(n*t)/(a0(t)*c(t)) + n*(sin(n*t)/sin(t));

  % build coefficient matrix for the An's
  coeffs = zeros(length(thetas), length(nvals));
  for i = 1:length(thetas) % row
    t = thetas(i);
    for j = 1:length(nvals) % col
      n = nvals(j);
      coeffs(i, j) = C(n, t);
    end
  end

  alphas = alpha(thetas)' - alpha0(thetas)';
  alphas = deg2rad(alphas);

  % Forier An terms
  As = coeffs\alphas;

  A1 = As(1);
  delta = 0;
  for i = 2:length(nvals)
    n = nvals(i);
    delta = delta + n*(As(i)/A1)^2;
  end
  e = (1 + delta)^-1;

  c_L  = A1*pi*AR;
  c_Di = c_L^2/(pi*e*AR);
end
