function [ A ] = update_A( X, Y, C, anchor, beta)

V = length(X);
Z = cell(V, 1);
A = cell(V, 1);
G = cell(V, 1);
%beta=ones(V,1);
for v=1:V
    temp = C{v}';
    [G{v},~,Z{v}] = svds(beta(v)^2*X{v}*temp(Y,:),anchor);
    %[G{v},~,Z{v}] = svds(X{v}*temp(Y,:),anchor);
    A{v} = G{v}*Z{v}';
end
end




