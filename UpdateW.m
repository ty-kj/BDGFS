function [A,v]=UpdateW(X,K,alpha,d)
%alpha: coefficient of l21
%X: num*dim
%d=nc-1: projection dimension
%A: projection matrix
%v: l2 norm of row


NIter=100;
% n=size(S,2);
% K=eye(n)-S-S'+S'*S;
dd=ones(size(X,2),1);
for iter=1:NIter
    %update A
    U=diag(dd);
    H=(X'*K*X+alpha*U);
%     [evec,eval]=eig(H);
%     [~,idx]=sort(diag(eval));
%     A=evec(:,idx(1:d));

    A=eig1(H,d,0);
    % update v: row-wise l_2 
    v=sqrt(sum(A.*A,2)+eps);
    dd=0.5./(v);
    
    obj(iter)=trace(A'*X'*K*X*A)+alpha*sum(v);
    if iter>1
        if abs(obj(iter)-obj(iter-1))<0.01
            break
        end
    end
end
%plot(obj);
end


