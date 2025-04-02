function [means] = meanInd(X, label)

%% X:  input 2D data:   d* n 
% label: spuerpixel index  m*n

[M,N]=size(label);   % HSI 尺寸
Class = unique(label); % 超像素类标
m = length(Class);     % 超像素个数
means = zeros(size(X,1), m);  % d * m,   m表示超像素个数
labelVector = reshape(label,[1,M*N]); % (M*N) * 1 

for i=1:length(Class)
    sub_idx = find(labelVector==Class(i));
    means(:,i) = mean(X(:,sub_idx),2); % the means of the pixles in one cluster
end

end
