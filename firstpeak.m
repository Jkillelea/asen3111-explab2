function idx = firstpeak(func, range)
  df = differentiate(func, range);
  for i = 1:length(df)
    if df(i) < 0
      idx = i;
      return;
    end
  end
  idx = -1;
end
