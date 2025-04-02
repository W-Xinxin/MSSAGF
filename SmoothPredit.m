%% smooth the prediction of labels in each ERS via the voting strategy 
%% ÿ��ERS ���е�Ԫ����ͬ��ģ������ͶƱ����������Ƶ���������ΪERS�������ص����
%% coded by xinxin 07/09/2023
function [smoothPre] = SmoothPredit(spLabel, preLabel)
% spLabel : superpixel index 2D M*N
% preLabel �� predicted label of pixels  n*1, n=M*N

[M,N]=size(spLabel);
Gt=reshape(spLabel,[1,M*N]);
Class=unique(Gt); % ���������
Num=size(Class,2);
numC = zeros(Num,1);  % record the number of classes per ERS blocks

for i=1:Num
    index{1,i}=find(Gt==Class(i));  % ĳ�������ؿ���Ԫ�ض�Ӧ������
    labels = preLabel(index{1,i});  % ĳ�������ؿ���Ԫ�ض�Ӧ��Ԥ�����
    ind = unique(labels);   %��������ؿ���Ԥ����������
    num = length(ind);  % ��������ؿ���Ԥ����������
    count = zeros(1,num);
    for j = 1: num
          count(1,j) = length(find(labels==ind(j))); % ͳ��Ԥ������Ƶ��
    end
   [~, order] = max(count); % ѡ��Ƶ���������
    MostLabel = ind(order);
    Gt(index{1,i}) = MostLabel; % Ƶ����ߵ������Ϊ���������ؿ������ص�Ԥ���ǩ
end
smoothPre = Gt;




