%% local embedding denoising in each ERS block
%% ||X-Y||_F^2 + beta * tr(Y'LY)
%% Y = (I+ lambda * L)^-1 X

function [dataNew] = LED(data,lambda,labels,d_nk)
%%
%  lambda :   the parameter
%  labels : superpixel index
%  data :  3D cube
%%

[nRow,nCol,dim] = size(data);
Results_segment= seg_im_class(data,labels);

Num = size(Results_segment.Y,2);
A = zeros(nRow*nCol,dim);

for i=1:Num
    sampleInERS = Results_segment.Y{1,i};
    W = Gen_Achor_Adj(sampleInERS',sampleInERS',d_nk,1);  % 10 
    L = diag(sum(W,2))- W;
    tmpY = (lambda*L + eye(size(sampleInERS,1)))^(-1)* sampleInERS;  %(A\B)  or A^-1B
    A(Results_segment.index{1,i},:) = tmpY;
end

dataNew = reshape(A,[nRow,nCol,dim]);




