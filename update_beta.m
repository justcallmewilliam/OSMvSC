function [ beta ] = update_beta( X, Y, Z, A )
r=2;
V = length(X);

loss = zeros(V, 1);

for v=1:V
    temp = Z{v}'*A{v}';
    loss(v) = sum(sum((X{v}'-temp(Y,:)).^2));
end

beta = loss.^(1/(1-r));
beta = loss/sum(loss);

end
