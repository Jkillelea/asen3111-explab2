function out = groupdata(str)
  nums = sscanf(str, 'ExperimentalLab2_Section%d_Group%d_Odd_Long.csv');

  if contains(str, 'Odd')
    angles = 'odd';
  else
    angles = 'even';
  end
  if contains(str, 'Long')
    type = 'long';
  else
    type = 'short';
  end

  out = struct(       ...
  'section', nums(1), ...
  'group',   nums(2), ...
  'type',    type,    ...
  'angles',  angles   ...
  );
end
