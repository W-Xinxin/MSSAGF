%% modified  from ICASSP2023£¬ STRUCTURED-ANCHOR PROJECTED CLUSTERING FOR HYPERSPECTRAL IMAGES
%% Guozhu jiang  yongshan zhang
%% coded by xinxinwang   30/06/2023
function [data2D_ds,labels, num_Pixel] = preprocessData2(data3D,dk)
%% HSI data preprocessing
% Input: 
%       data3D: 3D cube, HSI data.  multi-view datasets
%       dk: neighbors number
% Output:
%       newData:      new data , 2D matrix. each column is a pixel
%       labels_data:  superpixel labels of valid instances
%       labels:  superpixel labels  of entire 2D data block
%       num_Pixel: the number of superpixel
%%

numV = length(data3D);
for v=1 : numV
  [nRow,nCol,dim] = size(data3D{v});
  data2D{v} = reshape(data3D{v},nRow*nCol,dim);
end

newData2D = cell2mat(data2D);
%[X,~] = mapminmax(newData);
X = newData2D;

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

data2D_ds = cell(numV,1);  
for iv = 1 : numV
   newdata3D = S3_PCA(data3D{iv},dk,labels);
   [nRow,nCol,dim] = size(data3D{iv});
   temp = reshape(newdata3D,nRow*nCol,dim);
   data2D_ds{iv} = temp;
end

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