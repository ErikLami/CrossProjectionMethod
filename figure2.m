%% figure 2 -  Logarithmic neighbourhood sizes unilateral coupled logistic maps

clear 
close all
rng('default')

%% parameters
% parameters for logistic map
R = [3.82; 3.82];  
n = 10^4;                                 % length of time-series generated from logistic map
W = [0,0;-0.3,0];                         % coupling weights 
C = [0, 0; 0, 0];                         % correlation of internal noise
X0 = rand(1,2);                           % initial conditions

nData = 10^3;                             % amount of data used for causality estimation
embeddingDimension = 4;                   % embedding dimension
embeddingDelay = 1;                       % embedding time-delay
nReference = 10^3;                        % number of reference points 
iReference = randperm(nData, nReference); % time index of reference points
nNeighbors = 20;                          % number of neighbors    
temporalNeighbors = 'remove';             % 'remove' or 'keep' temporal neighbors
nEnsemble = 100;                          % number of ensembles for chance-level estimation
modeChance = 'temporal_shift';            % 'random' or 'temporal_shift'
nSE = 5;                                  % nSE-fold standard error


%% calculations
% generated data from logistic maps (with reflecting borders )
x = logistic_reflect(n, W, R, C, X0);

% remove transient (some additonal data is selected, to have nData points after time-delay embedding)
X = x(end-nData-3*(embeddingDimension-1)*embeddingDelay-1:end, 1);
Y = x(end-nData-3*(embeddingDimension-1)*embeddingDelay-1:end, 2);

% derive <log(neighborhoodsizes)>
out = manifold_inflation(X, Y, embeddingDelay, embeddingDimension, nNeighbors, iReference, temporalNeighbors, modeChance, nEnsemble);


%% visualizing 
figure()
subplot(1,2,1)
hold on
fill([(psi(1:nNeighbors)-log(nData)), flip(psi(1:nNeighbors)-log(nData))], [out.xSizexNeighbors + nSE*out.xSizexNeighbors_std/sqrt(nReference), flip(out.xSizexNeighbors - nSE*out.xSizexNeighbors_std/sqrt(nReference))], [.6 .6 .6], ...
    'linestyle', 'none', 'facealpha', 1);
fill([(psi(1:nNeighbors)-log(nData)), flip(psi(1:nNeighbors)-log(nData))], [out.xSizeyNeighbors + nSE*out.xSizexNeighbors_std/sqrt(nReference), flip(out.xSizeyNeighbors - nSE*out.xSizeyNeighbors_std/sqrt(nReference))], [.6 .6 .6], ...
    'linestyle', 'none', 'facealpha', 1);
fill([(psi(1:nNeighbors)-log(nData)), flip(psi(1:nNeighbors)-log(nData))], [out.xSizeChance + nSE*out.xSizeChance_std/sqrt(nReference*nEnsemble), flip(out.xSizeChance - nSE*out.xSizeChance_std/sqrt(nReference*nEnsemble))], [.6 .6 .6], ...
    'linestyle', 'none', 'facealpha', 1);
plot(psi(1:nNeighbors)-log(nData), out.xSizexNeighbors, 'k+' , 'MarkerSize', 20)
plot(psi(1:nNeighbors)-log(nData), out.xSizeyNeighbors, 'k.' , 'MarkerSize', 20)
plot(psi(1:nNeighbors)-log(nData), out.xSizeChance, 'k--')
xlabel('\kappa(k)')
ylabel('d_i^j(k)')
axis([-7.5 -3.5 -8 0.5])
text(-0.3, 1.1,'(a)','Units', 'Normalized', 'FontWeight', 'bold')

subplot(1,2,2)
hold on
fill([(psi(1:nNeighbors)-log(nData)), flip(psi(1:nNeighbors)-log(nData))], [out.ySizeyNeighbors + nSE*out.ySizeyNeighbors_std/sqrt(nReference), flip(out.ySizeyNeighbors - nSE*out.ySizeyNeighbors_std/sqrt(nReference))], [.6 .6 .6], ...
    'linestyle', 'none', 'facealpha', 1);
fill([(psi(1:nNeighbors)-log(nData)), flip(psi(1:nNeighbors)-log(nData))], [out.ySizexNeighbors + nSE*out.ySizexNeighbors_std/sqrt(nReference), flip(out.ySizexNeighbors - nSE*out.ySizexNeighbors_std/sqrt(nReference))], [.6 .6 .6], ...
    'linestyle', 'none', 'facealpha', 1);
fill([(psi(1:nNeighbors)-log(nData)), flip(psi(1:nNeighbors)-log(nData))], [out.ySizeChance + nSE*out.ySizeChance_std/sqrt(nReference*nEnsemble), flip(out.ySizeChance - nSE*out.ySizeChance_std/sqrt(nReference*nEnsemble))], [.6 .6 .6], ...
    'linestyle', 'none', 'facealpha', 1);
plot(psi(1:nNeighbors) -log(nData), out.ySizeyNeighbors, 'k+' , 'MarkerSize', 20)
plot(psi(1:nNeighbors) -log(nData), out.ySizexNeighbors, 'k.' , 'MarkerSize', 20)
plot(psi(1:nNeighbors) -log(nData), out.ySizeChance, 'k--')
xlabel('\kappa(k)')
axis([-7.5 -3.5 -8 0.5])
text(-0.3, 1.1,'(b)','Units', 'Normalized', 'FontWeight', 'bold')
