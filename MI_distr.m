rng('default') % for reproducibility of random variables

N = 20;

MI = N; 
MI_X = 0:0.5:N;
MI_Y = 0:0.5:N;

Red = zeros(1,length(MI_X));
Syn = zeros(1,length(MI_X));
Un_X = zeros(1,length(MI_X));
Un_Y = zeros(1,length(MI_X));

for n = 1:length(MI_X)
    Red(n) = min(MI_X(n), MI_Y(n));
    Un_X(n) = MI_X(n) - Red(n);
    Un_Y(n) = MI_Y(n) - Red(n);
    Syn(n) = MI - Red(n) - Un_X(n) - Un_Y(n);
end

fig = figure();
edges1=linspace(0, N, 20);
histogram(Syn,"BinEdges", edges1);
fig = figure();
histogram(Syn,"BinEdges", edges1);
fig = figure();
histogram(Syn,"BinEdges", edges1);
fig = figure();
histogram(Syn,"BinEdges", edges1);

% it obviously gives stupid results... think about what you do please...