function [] = ViolinPlot(UX,UY,R,S,MI,tit,save)

L = length(UX);
labels = strings(1,5*L);
labels(1,1:L) = "Ux";
labels(2,L+1:2*L) = "Uy";
labels(3,2*L+1:3*L) = "R";
labels(4,3*L+1:4*L) = "S";
labels(5,4*L+1:end) = "MI";

map = [0 0.4470 0.7410;
    0.6350 0.0780 0.1840;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880];
% colormap(map);

map2 = cell(5*L,1);
for l = 1:L
    map2{l} = map(1,:);
    map2{L+l} = map(2,:);
    map2{2*L+l} = map(3,:);
    map2{3*L+l} = map(4,:);
    map2{4*L+l} = map(5,:);
end
 
map = cell2mat(map2);

x = categorical(labels, ["MI", "Ux","Uy", "R","S"]);
y = [ UX UY R S MI];

fig = figure(); swarmchart(x,y,'.', 'SizeData',200, 'CData', map);
title(tit); 
ax = gca;
% ax.FontName = 'CMU serif'; 
set(gca,'FontName','CMU serif','FontSize',15);
yline(0, 'k--'); ylabel("Conscious-Unconscious",'FontSize',15,'interpreter','latex');
% xlabel("", 'interpreter','latex');
%saveas(fig, save);
exportgraphics(fig,save,'Resolution',300);

end