function [out] = mean_by_idxs(idxs, data)
  % Params
  % idxs -> Cell array containing vectors of the indexes to average together
  % data -> matrix, where every row is a unique data set

  out = zeros(length(idxs), size(data, 2)); % Average all the data based on the index groupings
  for i = 1:length(idxs)
    for j = 1:size(data, 2);
      out(i, j) = mean(data(idxs{i}, j));
    end
  end

end
