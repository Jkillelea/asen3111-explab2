function [short, long] = getfiles(data_dir)
  short = dir([data_dir '*Short.csv']);
  long  = dir([data_dir '*Long.csv']);
end
