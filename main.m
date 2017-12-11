clear;
clc;
close all;

data_dir = 'data/';
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

  idxs = group_into(alpha25, measurement_alphas, 0.5);
  empty_place = zeros(1, length(idxs));
  for j = 1:length(idxs)
    if isempty(idxs{j})
      empty_place(j) = true;
    end
  end
  idxs = idxs(~logical(empty_place));

  tmp = zeros(length(idxs), 3);
  for j = 1:length(idxs)
    tmp(j, 1) = mean(CL25(idxs{j}));
    tmp(j, 2) = mean(CD25(idxs{j}));
    tmp(j, 3) = mean(CMLE25(idxs{j}));
    tmp(j, 4) = mean(alpha25(idxs{j}));
  end
  CL25    = tmp(:, 1);
  CD25    = tmp(:, 2);
  CMLE25  = tmp(:, 3);
  alpha25 = tmp(:, 4);


  idxs = group_into(alpha15, measurement_alphas, 0.5);
  empty_place = zeros(1, length(idxs));
  for j = 1:length(idxs)
    if isempty(idxs{j})
      empty_place(j) = true;
    end
  end
  idxs = idxs(~logical(empty_place));

  tmp = zeros(length(idxs), 3);
  for j = 1:length(idxs)
    tmp(j, 1) = mean(CL15(idxs{j}));
    tmp(j, 2) = mean(CD15(idxs{j}));
    tmp(j, 3) = mean(CMLE15(idxs{j}));
    tmp(j, 4) = mean(alpha15(idxs{j}));
  end
  CL15    = tmp(:, 1);
  CD15    = tmp(:, 2);
  CMLE15  = tmp(:, 3);
  alpha15 = tmp(:, 4);


  figure('visible', 'off'); hold on; grid on; % 25 m/s plot
  plot(alpha25, CL25, 'displayname', 'C_L');
  plot(alpha25, CD25, 'displayname', 'C_D');
  plot(alpha25, CMLE25, 'displayname', 'C_{MLE}');
  title(sprintf('%s wing, %s angles, 25 m/s', ...
                file_info.type, file_info.angles));
  xlabel('Angle of Attack (degrees)');
  ylabel('C_l, C_d, C_m_{le} (unitless)');
  legend('show', 'location', 'southeast');
  print(['graphs/', gen_filename(file_info, 25)], '-dpng');

  figure('visible', 'off'); hold on; grid on; % 15 m/s plot
  plot(alpha15, CL15, 'displayname', 'C_L');
  plot(alpha15, CD15, 'displayname', 'C_D');
  plot(alpha15, CMLE15, 'displayname', 'C_{MLE}');
  title(sprintf('%s wing, %s angles, 15 m/s', ...
                file_info.type, file_info.angles));
  xlabel('Angle of Attack (degrees)');
  ylabel('C_l, C_d, C_m_{le} (unitless)');
  legend('show', 'location', 'southeast');
  print(['graphs/', gen_filename(file_info, 15)], '-dpng');

  close all;
end
