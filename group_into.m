function indicies = group_into(x, desired, tol)
  % groups points in 'x' by how close they are to the points in 'desired'
  %  params:
  %  x       -> vector of floats, raw data
  %  desired -> vector of floats, where we expect the points to show up
  %  tol     -> float, tolerance.
  indicies = cell(length(desired), 1);
  for i = 1:length(desired)
    d = desired(i);
    indicies{i} = find(abs(x - d) < tol);
  end
end
