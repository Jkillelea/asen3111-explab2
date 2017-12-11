function [out] = rm_empty_cells(in)
  empty_place = zeros(1, length(in));
  for j = 1:length(in) % remove empty indexes (no data is supposed to go in these points)
    if isempty(in{j})
      empty_place(j) = true;
    end
  end
  out = in(~logical(empty_place));
end
