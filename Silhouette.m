%% compute silhouette coefficient of each moduality
%% 想用聚类结果评价每个模态中ERS 划分性能
%% coded by xinxin 07/09/2023

function [score] = Silhouette( spLabel, preLabel)
% spLabel : superpixel index 2D M*N
% preLabel ： predicted label of pixels  n*1, n=M*N

[M,N]=size(spLabel);
Gt=reshape(spLabel,[1,M*N]);
Class=unique(Gt);
Num=size(Class,2);
numC = zeros(Num,1);  % record the number of classes per ERS blocks

for i=1:Num
    index{1,i}=find(Gt==Class(i));
    numC(i) = length(unique(preLabel(index{1,i}))); 
end
score = sum(numC)/Num ;

