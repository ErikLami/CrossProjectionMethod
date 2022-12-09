function Xm = embed(x, embeddingDelay, embeddingDimension) 
% time-delay reconstruction 
% x ... 1xn or nx1 time-series

if size(x,2) > size(x,1)
    x = x';
end

Xm = repmat(x,1,embeddingDimension);
for i=2:embeddingDimension
    Xm(:,i) = circshift(Xm(:,i),[embeddingDelay*(i-1),0]);
end

if embeddingDelay > 0
    Xm(1:(embeddingDimension-1)*embeddingDelay,:) = [];
elseif embeddingDelay < 0
    Xm(end+(embeddingDimension-1)*embeddingDelay:end,:) = [];
end

Xm = fliplr(Xm);
end
