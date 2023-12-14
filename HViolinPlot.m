function [] = HViolinPlot(UX_H,UY_H,R_H,S_H,MI_H, tit,save)

L_H = length(UX_H);
labels_H = strings(1,5*L_H);
labels_H(1,1:L_H) = "Ux";
labels_H(2,L_H+1:2*L_H) = "Uy";
labels_H(3,2*L_H+1:3*L_H) = "R";
labels_H(4,3*L_H+1:4*L_H) = "S";
labels_H(5,4*L_H+1:end) = "MI";

x_H = categorical(labels_H, ["MI" "Ux" "Uy" "R" "S"]);
y_H = [UX_H UY_H R_H S_H MI_H];

map = [0 0.4470 0.7410;
    0.6350 0.0780 0.1840;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880];

map2_H = cell(5*L_H,1);
for l = 1:L_H
    map2_H{l} = map(1,:);
    map2_H{L_H+l} = map(2,:);
    map2_H{2*L_H+l} = map(3,:);
    map2_H{3*L_H+l} = map(4,:);
    map2_H{4*L_H+l} = map(5,:);
end
 
map_H = cell2mat(map2_H);

fig = figure(); 
swarmchart(x_H,y_H,100,'o','CData', map_H, 'LineWidth',0.75); %red

%'Color','#0072BD',
% ,'Color','#A2142F'

title(tit);
yline(0, 'k--'); ylabel("Drug - Placebo",'FontSize',15,'interpreter','latex');
hold off

%saveas(fig, save);
set(gca,'FontName','CMU serif','FontSize',15);
exportgraphics(fig,save,'Resolution',300);

end