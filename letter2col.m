function col = letter2col(str)
% Change CSV spreadsheet columns 'A', 'B', 'C' ... 'AA', 'ZZ' ...
% to a zero-based index           0,   1,   2  ...  26,   675 ...
  col = 0;
  for i = 1:length(str)-1
    col = col + 26*(str(i) - 64); % AA is position 26
  end
  col = col + str(end) - 65;
end
