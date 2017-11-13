clear;
clc;
close all;

data_dir = 'data/';
% alphas = -11:20;

[short, long] = getfiles(data_dir);
data = [short; long];

for i = 1:length(data)
  fname = data(i).name;
  file_info = groupdata(fname);

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

  figure; hold on; grid on; % 25 m/s plot
  l = fit(alpha(speed25), L25, 'smoothingspline');
  d = fit(alpha(speed25), D25, 'smoothingspline');
  plot(alpha(speed25), feval(l, alpha(speed25)), 'displayname', 'Lift');
  plot(alpha(speed25), feval(d, alpha(speed25)), 'displayname', 'Drag');
  title(sprintf('%s wing, %s angles, 25 m/s', file_info.type, file_info.angles));
  xlabel('Angle of Attack (degrees)');
  ylabel('Force (N)');
  legend('show', 'location', 'southeast');

  filename = sprintf('section%d_group%d_%s_%s_25ms', ...
      file_info.section, file_info.group, file_info.angles, file_info.type);
  print(['graphs/', filename], '-dpng');


  figure; hold on; grid on; % 15 m/s plot
  l = fit(alpha(speed15), L15, 'smoothingspline');
  d = fit(alpha(speed15), D15, 'smoothingspline');
  plot(alpha(speed25), feval(l, alpha(speed25)), 'displayname', 'Lift');
  plot(alpha(speed25), feval(d, alpha(speed25)), 'displayname', 'Drag');
  title(sprintf('%s wing, %s angles, 15 m/s', file_info.type, file_info.angles));
  xlabel('Angle of Attack (degrees)');
  ylabel('Force (N)');
  legend('show', 'location', 'southeast');

  filename = sprintf('section%d_group%d_%s_%s_15ms', ...
      file_info.section, file_info.group, file_info.angles, file_info.type);
  print(['graphs/', filename], '-dpng');

  close all;
end
