function [out] = mean_by_idxs(idxs, data)
  % Params
  % idxs -> Cell array containing vectors of the indexes to average together
  % data -> matrix, where every row is a unique data set

  out = zeros(length(idxs), size(data, 2)); % Average all the data based on the index groupings
  for j = 1:length(idxs)
    out(j, 1) = mean(data(idxs{j}, 1));
    out(j, 2) = mean(data(idxs{j}, 2));
    out(j, 3) = mean(data(idxs{j}, 3));
    out(j, 4) = mean(data(idxs{j}, 4));
  end

end
