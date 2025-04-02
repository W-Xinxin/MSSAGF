function [X,labels,num_Pixel] = preData(data3D,lambda,d_nk)
%% HSI data preprocessing
% Input: 
%       data3D: 3D cube, HSI data.
%       lambda :   the number of neighbors for denoising (parameter for manifold term )
% Output:
%       X:      new data, 2D matrix. each column is a pixel
%       labels: superpixel labels
%       num_ERS: the number of superpixel
%%

[nRow,nCol,dim] = size(data3D);
X = reshape(data3D,nRow*nCol,dim);
[X,~] = mapminmax(X);

p = 1;
coeff = pca(X);
Y_pca = X*coeff(:,1:p);

img = im2uint8(mat2gray(reshape(Y_pca, nRow, nCol, p)));

Tbase = 2000;
[num_Pixel] = pixelNum(img,Tbase);
fprintf('Superpixels number : %d\n',num_Pixel);

%num_Pixel = num_ERS;  %%%%% added 

% ERS super-pixel segmentation.
labels = mex_ers(double(img),num_Pixel);
labels = labels + 1;

tic;
%% 不同的去噪方法
% newData = S3_PCA(data3D,k,labels); % zhangxin's method fo
  newData = LED(data3D,lambda,labels,d_nk);

time1 = toc;
fprintf('denoising time = %f\n',time1);

X = reshape(newData,nRow*nCol,dim);
X = X';
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