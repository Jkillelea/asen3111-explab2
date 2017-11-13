function [airspeed, rho, alpha, N, A, M] = load_csv(filename)
  raw_data = csvread(filename, 1, 0); % filename, row, column

  % A -> Atmo Pressure
  % B -> Atmo Temp
  % C -> Atmo Density
  % D -> Airspeed
  % E -> Pitot q (unused?)
  % F -> Aux q (unused?)
  % G thru V -> Scanivalve Pressures 1 thru 16
  % W -> Alpha (degrees)
  % X -> Sting Normal Force (N)
  % Y -> Sting Axial Force (N)
  % Z -> Sting pitching moment (N*m)
  % AA -> ELD Probe X (mm)
  % AB -> ELD Probe Y (mm)

  % matlab is 1-indexed but the CSV is zero-indexed, hence the +1
  airspeed = raw_data(:, letter2col('D')+1);
  rho      = raw_data(:, letter2col('C')+1);
  alpha    = raw_data(:, letter2col('W')+1);
  N        = raw_data(:, letter2col('X')+1);
  A        = raw_data(:, letter2col('Y')+1);
  M        = raw_data(:, letter2col('Z')+1);
end
