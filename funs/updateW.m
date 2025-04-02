function newW = updateW(X,A,Z,St,r)

H2 = X*Z*A';
H3 = H2+H2';
Q = St - H3 + A*diag(sum(Z))*A';
Q = (Q+Q')/2;
H = St\Q;
W = eig1(H,r,0,0);
newW = W*diag(1./sqrt(diag(W'*W)));

end