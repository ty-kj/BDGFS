function U  = UpdateU(Z,U,k,gammaoverlambda)  
n = size(Z,2);
one = ones(n,1);

L = diag(U*one)-U;    
% update W
[V, D] = eig(L);
D = diag(D);
[~, ind] = sort(D);    
W = V(:,ind(1:k))*V(:,ind(1:k))';


%Uk = U;
U = Z-gammaoverlambda*(repmat(diag(W),1,n)-W);
U = max(0,(U+U')/2);
U = U-diag(diag(U));
    
end