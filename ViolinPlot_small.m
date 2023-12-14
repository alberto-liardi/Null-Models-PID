function [] = ViolinPlot_small(X1,X2,X3,save,tit,label)

L = length(X1);
labels = strings(1,3*L);

labels(1,1:L) = label(1);
labels(2,L+1:2*L) = label(2);
labels(3,2*L+1:3*L) = label(3);

x = categorical(labels, label);
y = [ X1 X2 X3 ];

fig = figure(); swarmchart(x,y,'o','filled');
title(tit); yline(0, 'k--'); ylabel("Bias");

% set(fig, 'MenuBar', 'none');
% set(fig, 'ToolBar', 'none');

if ~isempty(save)
    saveas(fig, save);
end

end