% Averages the xs and 'ds' for all runs. The ds are the lengths of the 
% projections in the yz plane.
function [xs, ds] = processEndToEndVects(allXs, allYs, allZs)

xs = mean(allXs')';
ds = mean(sqrt(allYs'.^2 + allZs'.^2))';
