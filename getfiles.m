function [short, long, elliptical] = getfiles(data_dir)
  short       = dir([data_dir '*Short.csv']);
  long        = dir([data_dir '*Long.csv']);
  elliptical  = dir([data_dir '*Elliptical.csv']);
end
