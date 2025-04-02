%% coded by xinxinwang   07/07/2023
function [newData] = graphfilter(Data,kernel,k)
% Input: 
%       Data: 2D cube, HSI data.  num * dim
%       kernel:  graph filter with low-pass function
%       k :   k-order graph filter
% Output:
%       newData:      filtered data

%%
newData = kernel^k * Data;
