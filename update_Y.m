function [ Y ] = update_Y( X, A, C, k, beta)

V = length(X);
n = size(X{1}, 2);

AtCt = cell(V, 1);
for v=1:V
    AtCt{v} = C{v}'*A{v}';
end

loss = zeros(n, k);
for v=1:V
    loss = loss + beta(v)^2 * EuDist2(X{v}', AtCt{v}, 0);
end
[~, Y] = min(loss, [], 2);

end

