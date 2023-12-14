


save = [MI_list; PIDs(:,:,2)];
save = save';
T = array2table(save);
T.Properties.VariableNames(1:5) = {'MI','Unique X', 'Unique Y', 'Redundancy', 'Synergy'};
save22csv(T,"../Null_model_figures/Different_PID_study/", "PIDs_Iccs");


% writematrix(null_Reds,"../Null_model_figures/max_discrete/max_redundancy.csv")
% writematrix(null_Syns,"../Null_model_figures/max_discrete/max_synergy.csv")
% writematrix(null_Uns,"../Null_model_figures/max_discrete/max_unique.csv")

%% random quantile calculation for presentation

rng('default')
MI = 3.0;
null_PIDs = MI_null_model_Gauss(MI);
null_Reds = null_PIDs(3,:);
null_Syns = null_PIDs(4,:);
null_Uns = null_PIDs(1,:) + null_PIDs(2,:);

fig = figure();
histogram(null_Reds);
hold on
xline(mean(null_Reds),'--', 'LineWidth', 2, 'Color', 'black');
xlabel('Redundancies');
ylabel('Counts');
title("Redundancy distribution for MI = "+num2str(MI));
legend('Red', "mean = "+sprintf('%.2f',mean(null_Reds)));

fig = figure();
histogram(null_Syns);
hold on
xline(mean(null_Syns),'--', 'LineWidth', 2, 'Color', 'black');
xlabel('Synergies');
ylabel('Counts');
title("Synergy distribution for MI = "+num2str(MI));
legend('Red', "mean = "+sprintf('%.2f',mean(null_Syns)));

fig = figure();
histogram(null_Uns);
hold on
xline(mean(null_Uns),'--', 'LineWidth', 2, 'Color', 'black');
xlabel('Unique information');
ylabel('Counts');
title("Unique information distribution for MI = "+num2str(MI));
legend('Un', "mean = "+sprintf('%.2f',mean(null_Uns)));

mean(null_Reds)
mean(null_Syns)
mean(null_Uns)

function [] = save22csv(data, dir, name)

    writetable(data,dir+name+".csv");

end
