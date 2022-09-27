function [Y] = OSMvSC(X, label,dim,lambda)
%dataset(v): d_v*N(label_m+unlabel_n)
    if isempty(lambda)
    lambda=0;
    end
    iter = 1;
    thrsh = 0.00001;
    mu = 1e-4;
    pho = 1.1;
    max_mu = 1e6;
    view_num = length(X);
    IsConverge = 0;
    max_iter = 100;
    N = size(X{1},2);
    k=length(unique(label));
   
    beta = ones(view_num,1)./view_num;
    d=zeros(view_num,1);
    A=cell(view_num,1);
    C=cell(view_num,1);
    Z=cell(view_num,1);
    H=cell(view_num,1);
    Y = randi([1,k],N,1);
    
    for i = 1 : view_num
       d(i)=size(X{i},1);
       A{i}=rand(d(i),dim);
       Z{i}=zeros(dim,k);
       H{i}=zeros(dim,k);
    end

    while (IsConverge == 0&&iter<max_iter+1)

        [obj(iter)] = cal_obj(X, Y, Z, A,beta,lambda);              
        [A] = update_A( X, Y, Z, dim, beta);
        for i = 1: view_num
            [C{i}] = softth(Z{i}-H{i}/mu,lambda/mu);
        end
        [Z] = update_Z( X, A, Y, H, C, k, beta,mu);
        
        [beta] = update_beta( X, Y, Z, A);
        [Y] = update_Y( X, A, Z, k, beta);
       
        count = 0;
        for i = 1: view_num
            H{i} = H{i} + mu * (C{i}-Z{i});
            if norm(C{i}-Z{i},inf)<thrsh
                count = count + 1;
            end
        end
        
        mu = min(pho*mu, max_mu);
        
        if iter>2&&abs(obj(iter)-obj(iter-1))/obj(iter)<thrsh
            IsConverge = 1;
        end
        iter = iter + 1;
    end
end

function loss_sum = cal_obj(X, Y, C, A, beta,lambda)
    V = length(X);
    loss = zeros(V, 1);
    for v=1:V
        temp = C{v}'*A{v}';
        loss(v) = sum(sum((X{v}'-temp(Y,:)).^2))+lambda*trace(sqrt(C{v}'*C{v}));
    end
    loss_sum = beta'*loss;
end


