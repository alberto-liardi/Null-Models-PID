%% 3-variable case

% Covariance Matrix is:
%         1    a    c 
% Sigma = a    1    b
%         c    b    1

rng('default') % for reproducibility of random variables

% Fix total Mutual Information

I = rand(1);
% I = 2.5397; 

% CASE (a): fixing (b,c), calculating a

% create a grid of (b,c)
B = -0.999:0.001:0.999;
C = -0.999:0.001:0.999;

% b = 0.4; 
% c = 0.1;
I_list = 1:0.5:3;
Synergy = NaN(length(I_list),length(B),length(C));
num_bin = 0;

count = 0;
count2 = 0;
for i = 1:length(I_list)
    
    I = I_list(i);

    C = -sqrt(1-exp(-2*I)):0.001:sqrt(1-exp(-2*I));

    % preparing variables
    Unique_X = zeros(length(B),length(C));
    Unique_Y = zeros(length(B),length(C));
    Redundancy = zeros(length(B),length(C));
    I_X = zeros(length(B),length(C));
    I_Y = zeros(length(B),length(C));
    % CSER = zeros(length(B),length(C));
    % corr = zeros(length(B),length(C));
    
    A = zeros(length(B),length(C));
    
    for n = 1:length(B)
            for l = 1:length(C)
%                 disp(n); disp(l);
                
                b = B(n); c = C(l);

                % calculate a
                delta = b^2*c^2-(b^2+c^2-1+(1-b^2)*exp(-2*I));
                if delta < -1e-10
                    count2 = count2 + 1;
                    continue;
                end
                a = b*c + sqrt(delta);
%                 a = b*c - sqrt(delta);
                A(n,l) = a;
    
                if(1-(a^2+b^2+c^2)+2*a*b*c<1e-12)
                    count = count + 1; 
                    continue
                end
    
                I_X(n,l) = 0.5*log(1/(1-a^2));
                I_Y(n,l) = 0.5*log(1/(1-c^2));
                Redundancy(n,l) = min(I_X(n,l), I_Y(n,l));
                Unique_X(n,l) = I_X(n,l) - Redundancy(n,l);
                Unique_Y(n,l) = I_Y(n,l) - Redundancy(n,l);
                Synergy(i,n,l) = I - Redundancy(n,l) - Unique_Y(n,l) - Unique_X(n,l);   
%                 Synergy(i,n,l) = 0.5*log((1-b^2)*(1-max(abs(a),abs(c))^2)/(1-(a^2+b^2+c^2)+2*a*b*c));
    
            end
    end
    
    % heatmaps(C,B,Unique_X,Unique_Y,Synergy,Redundancy,'c','b');
    
    % histogram(Synergy(i,:));
    [~,edges] = histcounts(Synergy(i,:));
    if length(edges)>num_bin, num_bin=length(edges); end
    
end

disp("done");
disp(count2);

X = zeros(length(I_list),num_bin);
Y = zeros(length(I_list),num_bin);
for j = 1:length(I_list)
    Y(j,:) = I_list(j);
end
Z = zeros(length(I_list),num_bin);

for i = 1:length(I_list)
    
    [N,edges] = histcounts(Synergy(i,:), num_bin);
    N = N./sum(N);              % normalise w.r.t. N_tot
    centers = movmean(edges,2);
    X(i,:) = centers(2:end);
    Z(i,:) = N;

end

%surf(real(C),real(B),real(A));
% xlabel("c");
% ylabel("b");
% zlabel("a");

% waterfall
fig = figure();
waterfall(X,Y,Z);
xlabel('Synergy');
ylabel("Mutual Information");
zlabel("Counts");
title("Waterfall with (b,c) mesh",'FontSize',15,'interpreter','latex');
save = "Waterfall_bc.png";
exportgraphics(fig,save,'Resolution',300);
%%
% calculate and plot means
fig = figure();
averages = zeros(1,length(I_list));
for i = 1:length(I_list)
    averages(i) = mean(Synergy(i,:),'omitnan');
end

plot(I_list, averages, 'LineWidth',1);
title("Averages with (b,c) mesh",'FontSize',15,'interpreter','latex');
ylabel("Synergy averages");
xlabel("Mutual Information");
save = "Averages_bc.png";
exportgraphics(fig,save,'Resolution',300);
% hold on
% fplot(@(x) x,[0,10]);
% fimplicit(@(x,y)  - 0.8*x.^2 + 0.8*(y+1).^2 - 1);

%%
x = I_list;
y = averages;

hyprb = @(b,x) b(4)+b(1)*sqrt((x-b(2)).^2/b(3) +1);                 % Generalised Hyperbola
NRCF = @(b) norm(y - hyprb(b,x));                       % Residual Norm Cost Function
B0 = [1; 1; 1; 1];
options = optimset('MaxFunEvals',2e3);
B = fminsearch(NRCF, B0, options);                               % Estimate Parameters
figure(1)
plot(x, y)
hold on
plot(x, hyprb(B,x), '-r')
hold off
grid
% h = chi2gof(averages,hyprb(B,x));
%%


% CASE (b): fixing (a,c), calculating b

% create a grid of (a,c)
A = -0.999:0.001:0.999;
C = -0.999:0.001:0.999;

% a = 0.4; 
% c = 0.1;

I_list = 1:0.1:3;
Synergy = NaN(length(I_list),length(A),length(C));
num_bin = 0;

count = 0;
count2 = 0;
for i = 1:length(I_list)
    
    I = I_list(i);

    % create a grid of (a,c)
    A = -sqrt(1-exp(-2*I)):0.001:sqrt(1-exp(-2*I));
    C = -sqrt(1-exp(-2*I)):0.001:sqrt(1-exp(-2*I));

    % preparing variables
    Unique_X = zeros(length(A),length(C));
    Unique_Y = zeros(length(A),length(C));
    Redundancy = zeros(length(A),length(C));
    I_X = zeros(length(A),length(C));
    I_Y = zeros(length(A),length(C));
    % CSER = zeros(length(A),length(C));
    % corr = zeros(length(A),length(C));

    B = zeros(length(A), length(C));
    
    for n = 1:length(A)
            for l = 1:length(C)
    %             disp(n); disp(l);
                
                a = A(n); c = C(l);
    
                % calculate a
                delta = a^2*c^2-(a^2+c^2-1+exp(-2*I))*(1-exp(-2*I));
                if delta < -1e-10
%                     disp("delta not working");
%                     disp(count2);
                    count2 = count2 + 1;
                    continue;

                end
                b = (a*c + sqrt(delta))/(1-exp(-2*I));
%                 b = (a*c - sqrt(delta))/(1-exp(-2*I));
                B(n,l) = b;
    
                if(1-(a^2+b^2+c^2)+2*a*b*c<1e-12)
                    count = count + 1; 
    %                 disp("break for cycle");
    %                 disp(n); disp(l);
                    continue
                end
    
    %             disp("here");
                I_X(n,l) = 0.5*log(1/(1-a^2));
                I_Y(n,l) = 0.5*log(1/(1-c^2));
                Redundancy(n,l) = min(I_X(n,l), I_Y(n,l));
                Unique_X(n,l) = I_X(n,l) - Redundancy(n,l);
                Unique_Y(n,l) = I_Y(n,l) - Redundancy(n,l);
                Synergy(i,n,l) = I - Redundancy(n,l) - Unique_Y(n,l) - Unique_X(n,l);   
    %             S_theor = 0.5*log((1-b^2)*(1-max(abs(a),abs(c))^2)/(1-(a^2+b^2+c^2)+2*a*b*c))
    
            end

    end

    [~,edges] = histcounts(Synergy(i,:));
    if length(edges)>num_bin, num_bin=length(edges); end

end

disp("done");
disp(count2); % on the diagonal

% heatmaps(C,A,Unique_X,Unique_Y,Synergy,Redundancy,'c','a');

% histogram(Synergy(:));

X = zeros(length(I_list),num_bin);
Y = zeros(length(I_list),num_bin);
for j = 1:length(I_list)
    Y(j,:) = I_list(j);
end
Z = zeros(length(I_list),num_bin);

for i = 1:length(I_list)
    
    [N,edges] = histcounts(Synergy(i,:), num_bin);
    N = N./sum(N);              % normalise w.r.t. N_tot
    centers = movmean(edges,2);
    X(i,:) = centers(2:end);
    Z(i,:) = N;

end

% fig = figure();
%surf(real(A),real(C),real(B));
% xlabel("a");
% ylabel("c");
% zlabel("b");

% waterfall
fig = figure();
waterfall(X,Y,Z);
xlabel('Synergy');
ylabel("Mutual Information");
zlabel("Counts");
title("Waterfall with (a,c) mesh",'FontSize',15,'interpreter','latex');
save = "Waterfall_ac.png";
exportgraphics(fig,save,'Resolution',300);

%%
% calculate and plot means
fig = figure();
averages = zeros(1,length(I_list));
for i = 1:length(I_list)
    averages(i) = mean(Synergy(i,:),'omitnan');
end

plot(I_list, averages, 'LineWidth', 1);
% hold on
% fplot(@(x) 0.95*x-0.35,[1,3]);
ylabel("Synergy averages");
xlabel("Mutual Information");
save = "Averages_ac.png";
exportgraphics(fig,save,'Resolution',300);


function [] = heatmaps(k,g,Unique_X,Unique_Y,Synergy,Redundancy,x_lab,y_lab)
%     path = "Results/Minimal_Example_Theor/";

    kkk = [k(1), k(5:5:end-4),k(end)];
    xlabel = string(kkk);
    kind = 1:4:length(kkk);
    xlabelfin = strings(1,length(kkk));
    xlabelfin(kind) = xlabel(kind);
    dp = compose("%.2f", xlabelfin);
    xlabelfin(kind) = dp(kind);

    ggg = flip([g(1), g(5:5:end-4),g(end)]);
    ylabel = string(ggg);
    gind = 1:4:length(ggg);
    ylabelfin = strings(1,length(ggg));
    ylabelfin(gind) = ylabel(gind);
    dp = compose("%.2f", ylabelfin);
    ylabelfin(gind) = dp(gind);

    Syn_lab = "Synergy ";
    Red_lab = "Redundancy ";
    Un_x_lab = "Unique X ";
    Un_y_lab = "Unique Y ";
    
    % MMI
    fig = figure('Position', [100 100 550 500]);
    ff = vertcat(Unique_X(1,1), Unique_X(5:5:end-4,1), Unique_X(end,1));
    ss = vertcat(Unique_X(1,5:5:end-4), Unique_X(5:5:end-4,5:5:end-4), Unique_X(end,5:5:end-4));
    ll = vertcat(Unique_X(1,end), Unique_X(5:5:end-4,end), Unique_X(end,end));
    b = heatmap(kkk,ggg, flip([ff, ss, ll]), 'Colormap',jet);
    b.Title = Un_x_lab+'MMI';
    grid off;
    b.XLabel = x_lab;
    b.YLabel = y_lab;
%     save = path+Un_x_lab+'MMI.pdf';
    b.XDisplayLabels = xlabelfin;        
    b.YDisplayLabels = ylabelfin;
    set(gca,'FontName','CMU serif','FontSize',15);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
%     exportgraphics(fig,save,'Resolution',300);

    fig = figure('Position', [100 100 550 500]);
    ff = vertcat(Unique_Y(1,1), Unique_Y(5:5:end-4,1), Unique_Y(end,1));
    ss = vertcat(Unique_Y(1,5:5:end-4), Unique_Y(5:5:end-4,5:5:end-4), Unique_Y(end,5:5:end-4));
    ll = vertcat(Unique_Y(1,end), Unique_Y(5:5:end-4,end), Unique_Y(end,end));
    b = heatmap(kkk,ggg, flip([ff, ss, ll]), 'Colormap',jet);
    b.Title = Un_y_lab+'MMI';
    grid off;
    b.XLabel = x_lab;
    b.YLabel = y_lab;
%     save = path+Un_x_lab+'MMI.pdf';
    b.XDisplayLabels = xlabelfin;        
    b.YDisplayLabels = ylabelfin;
    set(gca,'FontName','CMU serif','FontSize',15);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
%     exportgraphics(fig,save,'Resolution',300);

    fig = figure('Position', [100 100 550 500]);
    ff = vertcat(Synergy(1,1), Synergy(5:5:end-4,1), Synergy(end,1));
    ss = vertcat(Synergy(1,5:5:end-4), Synergy(5:5:end-4,5:5:end-4), Synergy(end,5:5:end-4));
    ll = vertcat(Synergy(1,end), Synergy(5:5:end-4,end), Synergy(end,end));
    b = heatmap(kkk,ggg, flip([ff, ss, ll]), 'Colormap',jet);
    b.Title = Syn_lab+'MMI';
    grid off;
    b.XLabel = x_lab;
    b.YLabel = y_lab;
%     save = path+Syn_lab+'MMI.pdf';
    b.XDisplayLabels = xlabelfin;        
    b.YDisplayLabels = ylabelfin;
    set(gca,'FontName','CMU serif','FontSize',15);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
%     exportgraphics(fig,save,'Resolution',300);
    
    fig = figure('Position', [100 100 550 500]);
    ff = vertcat(Redundancy(1,1), Redundancy(5:5:end-4,1), Redundancy(end,1));
    ss = vertcat(Redundancy(1,5:5:end-4), Redundancy(5:5:end-4,5:5:end-4), Redundancy(end,5:5:end-4));
    ll = vertcat(Redundancy(1,end), Redundancy(5:5:end-4,end), Redundancy(end,end));
    c = heatmap(kkk,ggg, flip([ff, ss, ll]), 'Colormap',jet);
    c.Title = Red_lab+ 'MMI';
    grid off;
    c.XLabel = x_lab;
    c.YLabel = y_lab;
%     save = path+Red_lab+'MMI.pdf';
    c.XDisplayLabels = xlabelfin;        
    c.YDisplayLabels = ylabelfin;
    set(gca,'FontName','CMU serif','FontSize',15);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
%     exportgraphics(fig,save,'Resolution',300);

end