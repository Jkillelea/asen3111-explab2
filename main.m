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

elliptical_data15 = [];
elliptical_data25 = [];
long_data15       = [];
long_data25       = [];
short_data15      = [];
short_data25      = [];

wing_areas = [square_mm_to_m(14400); % Long
              square_mm_to_m(10800); % Short
              square_mm_to_m(2700*pi)]; % Elliptical

long_ar       = ((240/1000)^2)/wing_areas(1);
short_ar      = ((180/1000)^2)/wing_areas(2);
elliptical_ar = ((180/1000)^2)/wing_areas(3);
chord         = 60/1000; % meters

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

  cdi = d(5) - d(0);
  cl  = l(5);
  dl_dalpha = differentiate(l, 0)*(180/pi);
  zero_lift_alpha = fzero(l, 0); % degrees

  if strcmp(file_info.type, 'long')
    long_data15(i, :) = [zero_lift_alpha, dl_dalpha, clmax, alpha_range(idx), cdi, cl];
  elseif strcmp(file_info.type, 'short')
    short_data15(i, :) = [zero_lift_alpha, dl_dalpha, clmax, alpha_range(idx), cdi, cl];
  else
    elliptical_data15(i, :) = [zero_lift_alpha, dl_dalpha, clmax, alpha_range(idx), cdi, cl];
  end

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

  cdi = d(5) - d(0);
  cl  = l(5);
  dl_dalpha = differentiate(l, 0)*(180/pi);
  zero_lift_alpha = fzero(l, 0); % degrees

  if strcmp(file_info.type, 'long')
    long_data25(i, :) = [zero_lift_alpha, dl_dalpha, clmax, alpha_range(idx), cdi, cl];
  elseif strcmp(file_info.type, 'short')
    short_data25(i, :) = [zero_lift_alpha, dl_dalpha, clmax, alpha_range(idx), cdi, cl];
  else
    elliptical_data25(i, :) = [zero_lift_alpha, dl_dalpha, clmax, alpha_range(idx), cdi, cl];
  end

  close all;
end

disp(' ----------------------------------- ');

indexes = short_data15 ~= zeros(1, size(short_data15, 2));
zerolift   = mean(short_data15(indexes(:, 1), 1));
liftslope  = mean(short_data15(indexes(:, 2), 2));
stallangle = mean(short_data15(indexes(:, 4), 4));
clmax      = mean(short_data15(indexes(:, 3), 3));
cdi        = mean(short_data15(indexes(:, 3), 5));
cl         = mean(short_data15(indexes(:, 3), 6));
e = cl^2/(pi * cdi * short_ar);
disp('Short Wing 15');
fprintf('Zero lift %.2f degrees\n', zerolift);
fprintf('lift slope %.2f\n', liftslope);
fprintf('stall angle %.2f degrees\n', stallangle);
fprintf('clmax %.2f\n', clmax);
fprintf('cdi %.4f\n', cdi);
fprintf('e %.4f\n', e);

indexes = short_data25 ~= zeros(1, size(short_data25, 2));
zerolift   = mean(short_data25(indexes(:, 1), 1));
liftslope  = mean(short_data25(indexes(:, 2), 2));
stallangle = mean(short_data25(indexes(:, 4), 4));
clmax      = mean(short_data25(indexes(:, 3), 3));
cdi        = mean(short_data25(indexes(:, 3), 5));
cl         = mean(short_data25(indexes(:, 3), 6));
e = cl^2/(pi * cdi * short_ar);
disp('Short Wing 25');
fprintf('Zero lift %.2f degrees\n', zerolift);
fprintf('lift slope %.2f\n', liftslope);
fprintf('stall angle %.2f degrees\n', stallangle);
fprintf('clmax %.2f\n', clmax);
fprintf('cdi %.4f\n', cdi);
fprintf('e %.4f\n', e);

indexes = long_data15 ~= zeros(1, size(long_data15, 2));
zerolift   = mean(long_data15(indexes(:, 1), 1));
liftslope  = mean(long_data15(indexes(:, 2), 2));
stallangle = mean(long_data15(indexes(:, 4), 4));
clmax      = mean(long_data15(indexes(:, 3), 3));
cdi        = mean(long_data15(indexes(:, 3), 5));
cl         = mean(long_data15(indexes(:, 3), 6));
e = cl^2/(pi * cdi * long_ar);
disp('Long Wing 15');
fprintf('Zero lift %.2f degrees\n', zerolift);
fprintf('lift slope %.2f\n', liftslope);
fprintf('stall angle %.2f degrees\n', stallangle);
fprintf('clmax %.2f\n', clmax);
fprintf('cdi %.4f\n', cdi);
fprintf('e %.4f\n', e);

indexes = long_data25 ~= zeros(1, size(long_data25, 2));
zerolift   = mean(long_data25(indexes(:, 1), 1));
liftslope  = mean(long_data25(indexes(:, 2), 2));
stallangle = mean(long_data25(indexes(:, 4), 4));
clmax      = mean(long_data25(indexes(:, 3), 3));
cdi        = mean(long_data25(indexes(:, 3), 5));
cl         = mean(long_data25(indexes(:, 3), 6));
e = cl^2/(pi * cdi * long_ar);
disp('Long Wing 25');
fprintf('Zero lift %.2f degrees\n', zerolift);
fprintf('lift slope %.2f\n', liftslope);
fprintf('stall angle %.2f degrees\n', stallangle);
fprintf('clmax %.2f\n', clmax);
fprintf('cdi %.4f\n', cdi);
fprintf('e %.4f\n', e);

indexes = elliptical_data15 ~= zeros(1, size(elliptical_data15, 2));
zerolift   = mean(elliptical_data15(indexes(:, 1), 1));
liftslope  = mean(elliptical_data15(indexes(:, 2), 2));
stallangle = mean(elliptical_data15(indexes(:, 4), 4));
clmax      = mean(elliptical_data15(indexes(:, 3), 3));
cdi        = mean(elliptical_data15(indexes(:, 3), 5));
cl         = mean(elliptical_data15(indexes(:, 3), 6));
e = cl^2/(pi * cdi * elliptical_ar);
disp('Elliptical Wing 15');
fprintf('Zero lift %.2f degrees\n', zerolift);
fprintf('lift slope %.2f\n', liftslope);
fprintf('stall angle %.2f degrees\n', stallangle);
fprintf('clmax %.2f\n', clmax);
fprintf('cdi %.4f\n', cdi);
fprintf('e %.4f\n', e);

indexes = elliptical_data25 ~= zeros(1, size(elliptical_data25, 2));
zerolift   = mean(elliptical_data25(indexes(:, 1), 1));
liftslope  = mean(elliptical_data25(indexes(:, 2), 2));
stallangle = mean(elliptical_data25(indexes(:, 4), 4));
clmax      = mean(elliptical_data25(indexes(:, 3), 3));
cdi        = mean(elliptical_data25(indexes(:, 3), 5));
cl         = mean(elliptical_data25(indexes(:, 3), 6));
e = cl^2/(pi * cdi * elliptical_ar);
disp('Elliptical Wing 25');
fprintf('Zero lift %.2f degrees\n', zerolift);
fprintf('lift slope %.2f\n', liftslope);
fprintf('stall angle %.2f degrees\n', stallangle);
fprintf('clmax %.2f\n', clmax);
fprintf('cdi %.4f\n', cdi);
fprintf('e %.4f\n', e);


disp(' ----------------------------------- ');
disp('generated graphs in folder ./graphs/');

fclose(outfile);
