function [out] = manifold_inflation(X, Y, embeddingDelay, embeddingDimension, nNeighbors, iReference, temporalNeighbors, modeChance, nEnsemble)
%% Cross Projection Method
% Estimating Causal Influence from time series data by crossprojecting neighbors
% Input:
%  X                      nx1 or 1xn dimensional time series
%  Y                      second time series 
%  nNeighbors             number of neighbors used for cross projection 
%  embeddingDelay         time-delay of phase space reconstruction
%  embeddingDimension     dimension of phase space reconstruction
%  iReference             time-index of reference point(s) - if empty, all points are chosen as reference
%  temporalNeighbors      should temporal neighbors be considered? 'remove' to remove neighbors, default is no removeable 
%  modeChance             chance-level estimation from temporal shifted neighbors ('temporal_shift') or random neighbors ('random')
%  nEnsemble              number of Ensembles generated for chance-level estimation
% Output:
%  d_i|j                  average log neighborhoodsizes
  
%% reviewing input
% check time-series
if size(X, 2) > size(X, 1)  
    X = X';
    if size(X, 2) > 1
        error('X is not a Nx1 time-series.')
    end
end

if size(Y, 2) > size(Y, 1)
    Y = Y';
    if size(Y, 2) > 1
        error('Y is not a Nx1 time-series.') 
    end
end

% check length of timeseries & determine length after embedding
if length(X) ~= length(Y)
    error('Careful the two time series have mismatching length!')
end
n = length(X) - (embeddingDimension - 1)*embeddingDelay;

% check if reference-points are provided
%  if not provide every point as reference-point
%  if enough data-points are available for the provided reference-points
if isempty(iReference) 
    iReference = 1:n;
elseif length(iReference) > n
    error('Provide a correct number of reference-points')    
end
nReference = length(iReference);

% set radius for removal of temporal neighbors
if strcmp(temporalNeighbors, 'remove')
    cutRadius = (embeddingDimension - 1)*embeddingDelay;
else
    cutRadius = 0;
end

%% calculation
% time-delay reconstruction of the time-series
xEmbedded = embed(X, embeddingDelay, embeddingDimension);
yEmbedded = embed(Y, embeddingDelay, embeddingDimension);

% generate kd-tree
xEmb_kd = createns(xEmbedded);
yEmb_kd = createns(yEmbedded);

% increase number of neighbors, if temporal neighbors are removed
k = nNeighbors + 1 + 2*cutRadius;

% search neighbors
[xNeighbors] = knnsearch(xEmb_kd, xEmbedded, 'k', k);
[yNeighbors] = knnsearch(yEmb_kd, yEmbedded, 'k', k);

% delete temporal neighbors
xNeighbors_h = zeros(nReference, nNeighbors);
yNeighbors_h = zeros(nReference, nNeighbors);
for ref = 1:1:nReference
    xNeighbors_ref = xNeighbors(iReference(ref), :);
    yNeighbors_ref = yNeighbors(iReference(ref), :);    
    xNeighbors_ref = setdiff(xNeighbors_ref, iReference(ref)-cutRadius:iReference(ref)+cutRadius, 'stable'); 
    yNeighbors_ref = setdiff(yNeighbors_ref, iReference(ref)-cutRadius:iReference(ref)+cutRadius, 'stable');       
    xNeighbors_h(ref, :) = xNeighbors_ref(1:nNeighbors);
    yNeighbors_h(ref, :) = yNeighbors_ref(1:nNeighbors);
end   

% calcuclate neighborhoodsizes
xSizexNeighbors = neighborhoodsize(xEmbedded, iReference, xNeighbors_h, nReference, nNeighbors);
ySizeyNeighbors = neighborhoodsize(yEmbedded, iReference, yNeighbors_h, nReference, nNeighbors);
ySizexNeighbors = neighborhoodsize(yEmbedded, iReference, xNeighbors_h, nReference, nNeighbors);
xSizeyNeighbors = neighborhoodsize(xEmbedded, iReference, yNeighbors_h, nReference, nNeighbors);

% calculate chance-level
if strcmp(modeChance, 'random') % random neighbors
    xSizeChance = chance_level_random(xEmbedded, iReference, nEnsemble, nReference, nNeighbors, cutRadius, n);
    ySizeChance = chance_level_random(yEmbedded, iReference, nEnsemble, nReference, nNeighbors, cutRadius, n);
elseif strcmp(modeChance, 'temporal_shift') % neighbors from random temporal shift
    xSizeChance = chance_level_temporal_shift(xEmbedded, yNeighbors, iReference, nEnsemble, nReference, nNeighbors, cutRadius, n);
    ySizeChance = chance_level_temporal_shift(yEmbedded, xNeighbors, iReference, nEnsemble, nReference, nNeighbors, cutRadius, n);
end

% derive mean & standard deviation of logarithmic neighborhoodsize
out.xSizexNeighbors = mean(log(xSizexNeighbors)); 
out.ySizeyNeighbors = mean(log(ySizeyNeighbors)); 
out.ySizexNeighbors = mean(log(ySizexNeighbors)); 
out.xSizeyNeighbors = mean(log(xSizeyNeighbors)); 
out.xSizeChance = mean(log(xSizeChance)); 
out.ySizeChance = mean(log(ySizeChance)); 

out.xSizexNeighbors_std = std(log(xSizexNeighbors)); 
out.ySizeyNeighbors_std = std(log(ySizeyNeighbors)); 
out.ySizexNeighbors_std = std(log(ySizexNeighbors)); 
out.xSizeyNeighbors_std = std(log(xSizeyNeighbors)); 
out.xSizeChance_std = std(log(xSizeChance)); 
out.ySizeChance_std = std(log(ySizeChance)); 

end

function [s_ab] = neighborhoodsize(aEmbedded, iReference, bNeighbors, nReference, nNeighbors)
% calculate size of the neighboorhood in space a using neighborrelations from space b
s_ab = zeros(nReference, nNeighbors);
for ref = 1:1:nReference
    s_ab(ref, :) = cummax(sqrt(sum((aEmbedded(iReference(ref), :) - aEmbedded(bNeighbors(ref,:), :)).^2 ,2)));  
end
end  

function s_xy_c = chance_level_temporal_shift(xEmbedded, yNeighbors, iReference, nEnsemble, nReference, nNeighbors, cutRadius, n)
% calculate chance-level from (randomly) temporal shifted neighborhoods
s_xy_c = zeros(nReference*nEnsemble, nNeighbors);
for E = 1:1:nEnsemble
    randomShift = randi(n, 1);
    yNeighborsShifted = (yNeighbors + randomShift).*(yNeighbors + randomShift <= n) + (yNeighbors + randomShift - n).*(yNeighbors + randomShift > n);
    yNeighborsShifted_h = zeros(nReference, nNeighbors);
    for ref = 1:1:nReference
     	yNeighborsShifted_ref = setdiff(yNeighborsShifted(ref,:), iReference(ref) - cutRadius:iReference(ref) + cutRadius, 'stable');
        yNeighborsShifted_h(ref, :) = yNeighborsShifted_ref(1:nNeighbors);
    end
    yNeighborsShifted = yNeighborsShifted_h;   
    s_xy_c(1+(E-1)*nReference:E*nReference,:) = neighborhoodsize(xEmbedded, iReference, yNeighborsShifted, nReference, nNeighbors);
end    
end
