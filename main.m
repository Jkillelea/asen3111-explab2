clear;
clc;
close all;

outfile = fopen('results.txt', 'w');
if outfile < 0
  error('Unable to create an output file');
end

data_dir             = 'data/';
measurement_alphas   = -11:20;

[short, long, elliptical] = getfiles(data_dir);
data = [short; long; elliptical];

moment_arm = (6.67 + 72.75)/1000; % m

wing_areas = [square_mm_to_m(14400); % Long
              square_mm_to_m(10800); % Short
              square_mm_to_m(2700*pi)]; % Elliptical
chord = 60/1000; % meters

for i = 1:length(data)
  fprintf('%d / %d\n', i, length(data));
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
  CMLE15 = -M15/(q15*S*chord);
  CL25   = L25/(q25*S);
  CD25   = D25/(q25*S);
  CMLE25 = -M25/(q25*S*chord);

  alpha25 = alpha(speed25);
  alpha15 = alpha(speed15);

  % ==== ATTENTION ====
  % This next section looks pretty gross - so here's a summarry
  % We have 500 measurements at each alpha, all loosely grouped around specific alpha value
  % Since we know where those measurements are supposed to actually go, we can pick out the relevant indexes
  % And then pick out the cooresponding data and average it into one data point
  idxs    = rm_empty_cells(group_into(alpha25, measurement_alphas, 0.5)); % Get the indexes
  tmp     = mean_by_idxs(idxs, [CL25, CD25, CMLE25, alpha25]);
  CL25    = tmp(:, 1);
  CD25    = tmp(:, 2);
  CMLE25  = tmp(:, 3);
  alpha25 = tmp(:, 4);

  idxs    = rm_empty_cells(group_into(alpha15, measurement_alphas, 0.5)); % Get the indexes
  tmp     = mean_by_idxs(idxs, [CL15, CD15, CMLE15, alpha15]);
  CL15    = tmp(:, 1);
  CD15    = tmp(:, 2);
  CMLE15  = tmp(:, 3);
  alpha15 = tmp(:, 4);
  % BLAM! Instead of 500 data points sort of milling about each location (and 16 different locations), we now have only one at each

  % Create and save an assload of plots. Don't bother displaying them to the user (it gets in the way)
  figure('visible', 'off'); hold on; grid on; % 15 m/s plot
  l = fit(alpha15, CL15,   'smoothingspline');
  d = fit(alpha15, CD15,   'smoothingspline');
  m = fit(alpha15, CMLE15, 'smoothingspline');
  alpha_range = linspace(min(alpha15), max(alpha15), 10000);

  idx = firstpeak(l, alpha_range);
  clmax = l(alpha_range(idx));
  plot(alpha_range(idx), clmax, 'ro', 'displayname', 'Stall')

  plot(alpha_range, l(alpha_range), 'displayname', 'C_L');
  plot(alpha_range, d(alpha_range), 'displayname', 'C_D');
  plot(alpha_range, m(alpha_range), 'displayname', 'C_{MLE}');
  title(sprintf('%s wing, %s angles, 15 m/s', ...
                file_info.type, file_info.angles));
  xlabel('Angle of Attack (degrees)');
  ylabel('C_l, C_d, C_m_{le} (unitless)');
  legend('show', 'location', 'southeast');
  print(['graphs/', gen_filename(file_info, 15)], '-dpng');

  figure('visible', 'off'); hold on; grid on; % 25 m/s plot
  l = fit(alpha25, CL25,   'smoothingspline');
  d = fit(alpha25, CD25,   'smoothingspline');
  m = fit(alpha25, CMLE25, 'smoothingspline');
  alpha_range = linspace(min(alpha25), max(alpha25), 10000);

  idx = firstpeak(l, alpha_range);
  clmax = l(alpha_range(idx));
  plot(alpha_range(idx), clmax, 'ro', 'displayname', 'Stall')

  plot(alpha_range, l(alpha_range), 'displayname', 'C_L');
  plot(alpha_range, d(alpha_range), 'displayname', 'C_D');
  plot(alpha_range, m(alpha_range), 'displayname', 'C_{MLE}');


  title(sprintf('%s wing, %s angles, 25 m/s', ...
  file_info.type, file_info.angles));
  xlabel('Angle of Attack (degrees)');
  ylabel('C_l, C_d, C_m_{le} (unitless)');
  legend('show', 'location', 'southeast');
  print(['graphs/', gen_filename(file_info, 25)], '-dpng');

  close all;

  zero_lift_alpha = fzero(l, 0); % degrees
  fprintf(outfile, 'Section %d Group %d\n',         file_info.section, file_info.group);
  fprintf(outfile, '    Type: %s, %s angles\n',       file_info.type, file_info.angles);
  fprintf(outfile, '    Zero lift at alpha %.2f\n',   zero_lift_alpha);
  fprintf(outfile, '    Lift slope %.2f pi\n',        differentiate(l, 0)*(180/pi^2)); % report as per-radian
  fprintf(outfile, '    CL max %.2f at alpha %.2f\n', clmax, alpha_range(idx));
end

disp('generated graphs in folder ./graphs/');
disp('Results printed to output file ./results.txt');

fclose(outfile);
