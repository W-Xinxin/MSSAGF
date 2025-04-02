function [means] = meanInd(X, label)

%% X:  input 2D data:   d* n 
% label: spuerpixel index  m*n

[M,N]=size(label);   % HSI �ߴ�
Class = unique(label); % ���������
m = length(Class);     % �����ظ���
means = zeros(size(X,1), m);  % d * m,   m��ʾ�����ظ���
labelVector = reshape(label,[1,M*N]); % (M*N) * 1 

for i=1:length(Class)
    sub_idx = find(labelVector==Class(i));
    means(:,i) = mean(X(:,sub_idx),2); % the means of the pixles in one cluster
end

end
