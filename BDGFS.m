% min_W,Z||Y-YZ||_F,2+alpha||U||_k+\lambda||W||_2,1+\betaTr(W^TX^TLXW),
% \gamma*||U-Z||_F^2
% s.t. Z^t1=1,diag(Z)=0,W^tW=I

function [W,obj]=BDGFS(X,alpha,gamma,lambda,beta,nc)
iter_num = 100;%50;
[n, m] = size(X);  %[n d]

%Initialize W and Z
W=eye(m);
Z = ones(n,n)/n;
U =  ones(n,n)/n;
for iter = 1:iter_num
    F = X*W;
    Y = X';
    Z = UpdateZ(Z,Y,F,U,beta);%add U
    
    over = alpha/gamma;
    U = UpdateU(Z,U,nc,over);%block-diagonal

    % Fix Z,update W
    L = Laplacian(U);    
    [W,wc]=UpdateW(X,L,lambda/beta,floor(m/3));
    
    obj(iter) = norm(Y-Y*Z,'fro').^2+lambda*sum(wc)+beta*trace(W'*X'*L*X*W);
%     disp(['iter: ',num2str(iter),' obj:']);
%     disp(obj(iter));
end
% disp(obj);
