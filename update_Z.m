function [ Z ] = update_Z( X, A, Y, H, C, k, beta, mu)

V = length(X);
Z = cell(V,1);
XtA=cell(V,1);
%beta=ones(V,1);
for v=1:V
    XtA{v} = X{v}'*A{v};
    YtXtA=zeros(k,size(A{1},2));
    for i = 1: k
        YtXtA(i,:) = sum(XtA{v}(Y==i,:));
    end
    YtXtA = beta(v)^2*YtXtA;
    
    YtYvec=zeros(k,1);
    for i=1:k
       YtYvec(i)=length(find(Y==i)); 
    end
    YtYvec = 2*beta(v)^2*YtYvec + mu;
    %YtYvec = 2*YtYvec + mu;
    invYtYvec = 1./YtYvec;
    Z{v}=(2*YtXtA'+H{v}+mu*C{v}).*repmat(invYtYvec',size(H{1},1),1);
end
end


