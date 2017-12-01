clear;
clc;
close all;

data_dir = 'data/';
alphas   = -11:20;

[short, long, elliptical] = getfiles(data_dir);
data = [short; long; elliptical];

moment_arm = (6.67 + 72.75)/1000; % m

wing_areas = [square_mm_to_m(14400); % Long
              square_mm_to_m(10800); % Short
              square_mm_to_m(2700*pi)]; % Elliptical
chord = 60/1000; % meters

for i = 1:length(data)
  fname = data(i).name;
  file_info = groupdata(fname);

  if strcmp(file_info.type, 'long')
    S = wing_areas(1);
  elseif strcmp(file_info.type, 'short')
    S = wing_areas(2);
  else % elliptical
    S = wing_areas(3);
  end

  [airspeed, rho, alpha, N, A, M] = load_csv([data_dir, fname]);

  % make seleections based on airspeed
  speed0  = round(airspeed) < 5;
  speed15 = round(airspeed) < 17 & round(airspeed) >= 5;
  speed25 = round(airspeed) >= 17;

  % dynamic pressure (Pa)
  q15 = mean(0.5 * rho(speed15) .* airspeed(speed15).^2);
  q25 = mean(0.5 * rho(speed25) .* airspeed(speed25).^2);

  % remove zero-airspeed forces
  N15 = N(speed15) - N(speed0);
  A15 = A(speed15) - A(speed0);
  M15 = M(speed15) - M(speed0);
  N25 = N(speed25) - N(speed0);
  A25 = A(speed25) - A(speed0);
  M25 = M(speed25) - M(speed0);

  % lift and drag from normal and axial force (Newtons)
  L15 = N15.*cosd(alpha(speed15)) - A15.*sind(alpha(speed15));
  D15 = N15.*sind(alpha(speed15)) + A15.*cosd(alpha(speed15));
  L25 = N25.*cosd(alpha(speed25)) - A25.*sind(alpha(speed25));
  D25 = N25.*sind(alpha(speed25)) + A25.*cosd(alpha(speed25));

  MLE15 = M15 - N15*moment_arm; % N*m
  MLE25 = M25 - N25*moment_arm; % N*m

  CL15   = L15/(q15*S);
  CD15   = D15/(q15*S);
  CMLE15 = M15/(q15*S*chord);
  CL25   = L25/(q25*S);
  CD25   = D25/(q25*S);
  CMLE25 = M25/(q25*S*chord);

  figure; hold on; grid on; % 25 m/s plot
  l = fit(alpha(speed25), CL25,   'smoothingspline');
  d = fit(alpha(speed25), CD25,   'smoothingspline');
  m = fit(alpha(speed25), CMLE25, 'smoothingspline');
  plot(alpha(speed25), feval(l, alpha(speed25)), 'displayname', 'C_L');
  plot(alpha(speed25), feval(d, alpha(speed25)), 'displayname', 'C_D');
  plot(alpha(speed25), feval(m, alpha(speed25)), 'displayname', 'C_{MLE}');
  title(sprintf('%s wing, %s angles, 25 m/s', ...
                file_info.type, file_info.angles));
  xlabel('Angle of Attack (degrees)');
  ylabel('Force (N)');
  legend('show', 'location', 'southeast');
  print(['graphs/', gen_filename(file_info, 25)], '-dpng');

  figure; hold on; grid on; % 15 m/s plot
  l = fit(alpha(speed15), L15,   'smoothingspline');
  d = fit(alpha(speed15), D15,   'smoothingspline');
  m = fit(alpha(speed25), CMLE15, 'smoothingspline');
  plot(alpha(speed15), feval(l, alpha(speed15)), 'displayname', 'C_L');
  plot(alpha(speed15), feval(d, alpha(speed15)), 'displayname', 'C_D');
  plot(alpha(speed15), feval(m, alpha(speed15)), 'displayname', 'C_{MLE}');
  title(sprintf('%s wing, %s angles, 15 m/s', ...
                file_info.type, file_info.angles));
  xlabel('Angle of Attack (degrees)');
  ylabel('Force (N)');
  legend('show', 'location', 'southeast');
  print(['graphs/', gen_filename(file_info, 15)], '-dpng');

  close all;
end
