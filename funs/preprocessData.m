%% modified  from ICASSP2023， STRUCTURED-ANCHOR PROJECTED CLUSTERING FOR HYPERSPECTRAL IMAGES
%% Guozhu jiang  yongshan zhang
%% coded by xinxinwang   30/06/2023
function [newData,labels_data,labels, num_Pixel] = preprocessData(data3D,Ind)
%% HSI data preprocessing
% Input: 
%       data3D: 3D cube, HSI data.  multi-view datasets
%       Ind: valid index for instances
% Output:
%       newData:      new data of valid instances, 2D matrix. each column is a pixel
%       labels_data:  superpixel labels of valid instances
%       labels:  superpixel labels  of entire 2D data block
%       num_Pixel: the number of superpixel
%%

for v=1 : length(data3D)
  [nRow,nCol,dim] = size(data3D{v});
  data3D{v} = reshape(data3D{v},nRow*nCol,dim)
end

newData = cell2mat(data3D);
%[X,~] = mapminmax(newData);
X = newData;

p = 1;
coeff = pca(X);
Y_pca = X*coeff(:,1:p);

img = im2uint8(mat2gray(reshape(Y_pca, nRow, nCol, p)));

Tbase = 2000;
[num_Pixel] = pixelNum(img,Tbase);
fprintf('Superpixels number : %d\n',num_Pixel);

% ERS super-pixel segmentation.
labels = mex_ers(double(img),num_Pixel);
labels = labels + 1;

tic;
%newData = S3_PCA(data3D,k,labels);
% construct the affinity matrix for samples
labels_colA = reshape(labels,nRow*nCol,1);  % 超像素编号
Data = X(Ind,:);  % 有效样本
labels_data = labels_colA(Ind,:); % 样本的超像素编号
Sa = Cubseg_Gen_adj_2D( Data ,labels_data ); % 构造spatial-spectral graph
num = length(Ind); % 有效样本数量
Kernel = spdiags(ones(num,1),0,num,num) - 0.5* (diag(sum(Sa,2)) - Sa);
k =1 ;
newData = graphfilter(Data,Kernel,k); 

time1 = toc;
fprintf('denoising time = %f\n',time1);

end

function [num]=pixelNum(img,Tbase)
% Calculate the number of superpixels by RLPA
% https://github.com/junjun-jiang/RLPA
[m,n] = size(img);
% img =  rgb2gray(img);
BW = edge(img,'log');
% figure,imshow(BW);
ind = find(BW~=0);
Len = length(ind);
Ratio = Len/(m*n);
num = fix(Ratio * Tbase);
end