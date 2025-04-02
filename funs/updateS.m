function S = updateS(X, F, gamma, lambda, k)

num = size(X,2);

distf = L2_distance_1(F',F');
distX = L2_distance_1(X,X);
[~, idx] = sort(distX,2); %sort each row, increase

S = zeros(num);
for i=1:num
    idxa0 = idx(i,2:k+1);
    dfi = distf(i,idxa0);
    dxi = distX(i,idxa0);
    ad = -(dxi + lambda*dfi)/(2*gamma);
    S(i,idxa0) = EProjSimplex_new(ad);
end
S = (S+S')/2;

end