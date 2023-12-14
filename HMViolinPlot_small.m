function [] = HMViolinPlot_small(X1_M,X2_M,X3_M,X1_H,X2_H,X3_H,save,tit,label)

L_M = length(X1_M);
labels_M = strings(1,3*L_M);

labels_M(1,1:L_M) = label(1);
labels_M(2,L_M+1:2*L_M) = label(2);
labels_M(3,2*L_M+1:3*L_M) = label(3);

L_H = length(X1_H);
labels_H = strings(1,3*L_H);

labels_H(1,1:L_H) = label(1);
labels_H(2,L_H+1:2*L_H) = label(2);
labels_H(3,2*L_H+1:3*L_H) = label(3);

x_M = categorical(labels_M, label);
y_M = [ X1_M X2_M X3_M ];
x_H = categorical(labels_H, label);
y_H = [ X1_H X2_H X3_H ];

map = [0 0.4470 0.7410;
    0.6350 0.0780 0.1840;
    0.9290, 0.6940, 0.1250];

map2_H = cell(3*L_H,1);
for l = 1:L_H
    map2_H{l} = map(1,:);
    map2_H{L_H+l} = map(2,:);
    map2_H{2*L_H+l} = map(3,:);
end
 
map_H = cell2mat(map2_H);

map2_M = cell(3*L_M,1);
for l = 1:L_M
    map2_M{l} = map(1,:);
    map2_M{L_M+l} = map(2,:);
    map2_M{2*L_M+l} = map(3,:);
end
 
map_M = cell2mat(map2_M);

fig = figure(); 
swarmchart(x_M,y_M,100,'^','CData', map_M, 'LineWidth',0.75); %blue
hold on
swarmchart(x_H,y_H,100,'o','CData', map_H, 'LineWidth',0.75); %red

title(tit,'FontSize',15,'interpreter','latex'); 
yline(0, 'k--'); ylabel("Bias",'FontSize',15,'interpreter','latex');
LH(1) = plot(nan, nan,'black ^','LineWidth',0.75); %blue
L{1} = 'Monkeys';
LH(2) = plot(nan, nan, 'black o','LineWidth',0.75); %red
L{2} = 'Humans';
legend(LH, L, 'Location','southeast');
hold off

set(gca,'FontName','CMU serif','FontSize',15);
exportgraphics(fig,save,'Resolution',300);
% saveas(fig, save);

end