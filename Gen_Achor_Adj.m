function [S, D] = Gen_Achor_Adj(X,anchor,k, issymmetric)
%% code by wang xinxin 2022/08/20
% X: each column is a data point £ºd*n
% anchor: each column is a data point £º d*m
% k: number of neighbors
% issymmetric: set S = (S+S')/2 if issymmetric=1
% S: similarity matrix, each row is a data point
% Ref: F. Nie, X. Wang, M. I. Jordan, and H. Huang, The constrained
% Laplacian rank algorithm for graph-based clustering, in AAAI, 2016.

if nargin < 3
    issymmetric = 1;
end;
if nargin < 2
    k = 5;
end;

[~, n] = size(X);
[~, m] = size(anchor);
D = EuDist2(X', anchor',0);
[~, idx] = sort(D, 2); % sort each row

S = zeros(n,m);
for i = 1:n
    id = idx(i,2:k+2);
    di = D(i, id);
    S(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end;

if issymmetric == 1
    S = (S+S')/2;
end;