%% smooth the prediction of labels in each ERS via the voting strategy 
%% 每个ERS 块中的元素是同类的，因此用投票方法，将最频繁的类标作为ERS块中像素的类标
%% coded by xinxin 07/09/2023
function [smoothPre] = SmoothPredit(spLabel, preLabel)
% spLabel : superpixel index 2D M*N
% preLabel ： predicted label of pixels  n*1, n=M*N

[M,N]=size(spLabel);
Gt=reshape(spLabel,[1,M*N]);
Class=unique(Gt); % 超像素类别
Num=size(Class,2);
numC = zeros(Num,1);  % record the number of classes per ERS blocks

for i=1:Num
    index{1,i}=find(Gt==Class(i));  % 某个超像素块内元素对应的坐标
    labels = preLabel(index{1,i});  % 某个超像素块内元素对应的预测类标
    ind = unique(labels);   %这个超像素块内预测类标的类型
    num = length(ind);  % 这个超像素块内预测类标的数量
    count = zeros(1,num);
    for j = 1: num
          count(1,j) = length(find(labels==ind(j))); % 统计预测类标的频率
    end
   [~, order] = max(count); % 选出频率最大的类标
    MostLabel = ind(order);
    Gt(index{1,i}) = MostLabel; % 频率最高的类别作为整个超像素块内像素的预测标签
end
smoothPre = Gt;




