function inetsim(paramName, paramValues)


for k = 1:numel(paramValues)
  stats(k) = netsim(paramName, paramValues(k));  
end
% Close the waitbar

% Plot results
statsNames = fieldnames(stats);

for k = 1:numel(statsNames)
  statValues = zeros(1, numel(paramValues));
  for n = 1:numel(paramValues)
    statValues(n) = getfield(stats(n), statsNames{k});
  end
  figure
  plot(paramValues, statValues);
  ax(paramName, statsNames{k});  
end
 
