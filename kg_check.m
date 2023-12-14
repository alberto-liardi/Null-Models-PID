X = [1];
Y = [4];

A = [0.8147, 0.0126; 1.2699, 0.9058];

Z = vertcat(X, Y);

k = 0.99;
g = 100;
for t = 1:1000
    V = [g k; k g];
    res = mvnrnd([0,0], V,1);
    Z(:,end+1) = A*Z(:,end) + res';
end

fig = figure();
plot(1:1001, Z(1,:), 1:1001, Z(2,:));