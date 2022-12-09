%% figure 1 -  visualization of neighborhoodrelations for 2 (unilateral) coupled logistic maps 
% demonstrating the localisation of the neighbors in both origin spaces and the stronger 
% dispersion of neighbors in their respective projections  
% for additional visualisation the projection on the $y_t, y_{t+1}$-plane and 
% $x_t$-axis are shown 

clear 
close all
rng('default')

% parameters
R = [3.82; 3.82];
n = 10^6;                % number of generated data points
nData = 10^3;            % number of used data points
embeddingDimension = 2;  % embedding dimension
embeddingDelay = 1;      % embedding delay
nReference = 1;          % number of reference points
iReference = 7;          % index of refererence point
nNeighbors = 10;         % number of neighbors    
X0 = rand(1,2);          % initial conditions
C = [0,0;0,0];           % correlation of internal noise
W = [0,0;-0.3,0];        % coupling weights

% genrate data from logistic maps (with reflecting borders )
x = logistic_reflect(n, W, R, C, X0);

% remove transient
X = x(end-nData-3*(embeddingDimension-1)*embeddingDelay+1:end, 1);
Y = x(end-nData-3*(embeddingDimension-1)*embeddingDelay+1:end, 2);

% time-delay reconstruction
rx = embed(X, -embeddingDelay, embeddingDimension);
ry = embed(Y, -embeddingDelay, embeddingDimension);

% generate kd-tree
rx_c = createns(rx);
ry_c = createns(ry);

% find neighbors
[I_neix, distx_neix] = knnsearch(rx_c, rx(iReference, :), 'k', nData);
[I_neiy, disty_neiy] = knnsearch(ry_c, ry(iReference, :), 'k', nData);

%% visualize
figure()
hold on
scatter3(Y(1:end-1), Y(2:end), X(1:end-1), 2,'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerEdgeAlpha', .10, 'MarkerFaceAlpha', .10)
scatter3(Y(1:end-1), Y(2:end), -0.3+zeros(length(Y(2:end)),1), 10, 'k', 'filled', 'MarkerEdgeAlpha', .025, 'MarkerFaceAlpha', .025)
scatter3(ones(length(Y(2:end)),1), ones(length(Y(2:end)),1), X(1:end-1), 10, 'k', 'filled', 'MarkerEdgeAlpha', .025, 'MarkerFaceAlpha', .025)
scatter3(Y(I_neix(2:1+nNeighbors)), Y(I_neix(2:1+nNeighbors)+1), X(I_neix(2:1+nNeighbors)), 10, [0 1 0.7]*0.5, '*')
scatter3(Y(I_neiy(2:1+nNeighbors)), Y(I_neiy(2:1+nNeighbors)+1), X(I_neiy(2:1+nNeighbors)), 10, [1 0 0], 'diamond')
scatter3(Y(I_neix(2:1+nNeighbors)), Y(I_neix(2:1+nNeighbors)+1), -0.3+zeros(nNeighbors,1), 10, [0 1 0.7]*0.5, '*')
scatter3(Y(I_neiy(2:1+nNeighbors)), Y(I_neiy(2:1+nNeighbors)+1), -0.3+zeros(nNeighbors,1), 10, [1 0 0], 'diamond')
scatter3(ones(nNeighbors,1)-0.005, ones(nNeighbors,1)-0.005, X(I_neix(2:11)), 10, [0 1 0.7]*0.5, '*')
scatter3(ones(nNeighbors,1)-0.006, ones(nNeighbors,1)-0.006, X(I_neiy(2:11)), 10, [1 0 0], 'diamond')
zlabel('x_t', 'fontsize', 8)
xlabel('y_t', 'fontsize', 8)
ylabel('y_{t+1}', 'fontsize', 8)
zticks([0 0.5 1])
axis([0 1 0 1 -0.3 1])
set(gca, 'fontsize', 8)
grid('off')
view(70, 22)