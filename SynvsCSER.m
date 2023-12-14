%% CSER - Synergy plots

%% MONKEYS 

monkeys = ["Kin2", "Chibi", "Su", "George"];

fin_CSER_W = cell(1,4); fin_CSER_S = cell(1,4);
fin_Syn_W = cell(1,4); fin_Syn_S = cell(1,4);

for i = 1:4

    as = "Results/"+monkeys(i)+"_c20/Sleep/";
    aw = "Results/"+monkeys(i)+"_c20/Awake/";
    subject = monkeys(i);
    
    % read from file all the necessary info
    % PID atoms
    [Un_X_S, Un_Y_S, Red_S, Syn_S, MI_S] = MonkeyPIDLoader(as + "MMI/data.txt");
    [Un_X_W, Un_Y_W, Red_W, Syn_W, MI_W] = MonkeyPIDLoader(aw + "MMI/data.txt");
    
    % CSER
    InfoDyn_S = readmatrix(as+"InfoDyn.txt");
    InfoDyn_W = readmatrix(aw+"InfoDyn.txt");
    CSER_S = InfoDyn_S(:,1:2); 
    CSER_W = InfoDyn_W(:,1:2);
    
    % using only Var(p)
    Syn_S = Syn_S{2};
    Syn_W = Syn_W{2};
    CSER_S = CSER_S(:,2);
    CSER_W = CSER_W(:,2);
    
    % sorting by noise (CSER)
    [CSER_S, id_S] = sort(CSER_S);
    [CSER_W, id_W] = sort(CSER_W);
    Syn_S = Syn_S(id_S);
    Syn_W = Syn_W(id_W);

    % plotting
    % fig = figure();
    % plot(CSER_S, Syn_S);
    % fig = figure();
    % plot(CSER_W, Syn_W);
    fig = figure();
    plot(CSER_S, Syn_S, CSER_W, Syn_W);
    title(subject);
    xlabel('CSER');
    ylabel('Synergy');
    legend("Sleep/Sedated", "Awake");

    fin_CSER_W{i} = CSER_W; fin_CSER_S{i} = CSER_S;
    fin_Syn_W{i} = Syn_W; fin_Syn_S{i} = Syn_S;

end

% plotting
% Awake
fig = figure();
hold on;
for i = 1:4
    plot(fin_CSER_W{i}, fin_Syn_W{i}, 'LineWidth', 1);
end
title('Awake');
xlabel('CSER');
ylabel('Synergy');
legend('Kin2', 'Chibi', 'Su', 'George');

% Sedated
fig = figure();
hold on;
for i = 1:4
    plot(fin_CSER_S{i}, fin_Syn_S{i}, 'LineWidth', 1);
end
title('Sedated');
xlabel('CSER');
ylabel('Synergy');
legend('Kin2', 'Chibi', 'Su', 'George');

%% smoothing everything with moving average
 
% mm_CSER_W = cell(1,4); mm_CSER_S = cell(1,4);
mm_Syn_W = cell(1,4); mm_Syn_S = cell(1,4);

for i = 1:4
    mm_Syn_W{i} = movmean(fin_Syn_W{i},7);
    mm_Syn_S{i} = movmean(fin_Syn_S{i},7);
end

% plotting
% Awake
fig = figure();
hold on;
for i = 1:4
    plot(fin_CSER_W{i}, mm_Syn_W{i}, 'LineWidth', 1);
end
title('Monkeys awake','FontSize',15,'interpreter','latex');
xlabel('CSER');
ylabel('Synergy');
legend('Kin2', 'Chibi', 'Su', 'George');

% Sedated
fig = figure();
hold on;
for i = 1:4
    plot(fin_CSER_S{i}, mm_Syn_S{i}, 'LineWidth', 1);
end
title('Monkeys sedated','FontSize',15,'interpreter','latex');
xlabel('CSER');
ylabel('Synergy');
legend('Kin2', 'Chibi', 'Su', 'George');
save = "SynvsCSER_M.png";
exportgraphics(fig,save,'Resolution',300);

%% HUMANS 

humans = ["S2", "S4", "S5", "S8", "S11", "S12", "S13", "S15", "S16"];

fin_CSER_W = cell(1,9); fin_CSER_S = cell(1,9);
fin_Syn_W = cell(1,9); fin_Syn_S = cell(1,9);

for i = 1:9
    
    as = "Results/"+humans(i)+"_c20/Sleep/";
    if(~exist(as, 'dir'))
        as = "Results/"+humans(i)+"_c16/Sleep/";
        if(~exist(as, 'dir'))
            as = "Results/"+humans(i)+"_c12/Sleep/";
        end
    end
    
    aw = "Results/"+humans(i)+"_c20/Awake/";
    if(~exist(aw, 'dir'))
        aw = "Results/"+humans(i)+"_c16/Sleep/";
        if(~exist(aw, 'dir'))
            aw = "Results/"+humans(i)+"_c12/Sleep/";
        end
    end

    subject = humans(i);


    % read from file all the necessary info
    % PID atoms
    [Un_X_S, Un_Y_S, Red_S, Syn_S, MI_S] = MonkeyPIDLoader(as + "MMI/data.txt");
    [Un_X_W, Un_Y_W, Red_W, Syn_W, MI_W] = MonkeyPIDLoader(aw + "MMI/data.txt");
    
    % CSER
    InfoDyn_S = readmatrix(as+"InfoDyn.txt");
    InfoDyn_W = readmatrix(aw+"InfoDyn.txt");
    CSER_S = InfoDyn_S(:,1:2); 
    CSER_W = InfoDyn_W(:,1:2);
    
    % using only Var(p)
    Syn_S = Syn_S{2};
    Syn_W = Syn_W{2};
    CSER_S = CSER_S(:,2);
    CSER_W = CSER_W(:,2);
    
    % sorting by noise (CSER)
    [CSER_S, id_S] = sort(CSER_S);
    [CSER_W, id_W] = sort(CSER_W);
    Syn_S = Syn_S(id_S);
    Syn_W = Syn_W(id_W);

    % plotting
    % fig = figure();
    % plot(CSER_S, Syn_S);
    % fig = figure();
    % plot(CSER_W, Syn_W);
    fig = figure();
    plot(CSER_S, Syn_S, CSER_W, Syn_W);
    title(subject);
    xlabel('CSER');
    ylabel('Synergy');
    legend("Sleep/Sedated", "Awake");

    fin_CSER_W{i} = CSER_W; fin_CSER_S{i} = CSER_S;
    fin_Syn_W{i} = Syn_W; fin_Syn_S{i} = Syn_S;

end

% plotting
% Awake
fig = figure();
hold on;
for i = 1:9
    plot(fin_CSER_W{i}, fin_Syn_W{i}, 'LineWidth', 1);
end
title('Awake');
xlabel('CSER');
ylabel('Synergy');
legend("S2", "S4", "S5", "S8", "S11", "S12", "S13", "S15", "S16");

% Sedated
fig = figure();
hold on;
for i = 1:9
    plot(fin_CSER_S{i}, fin_Syn_S{i}, 'LineWidth', 1);
end
title('Sleep');
xlabel('CSER');
ylabel('Synergy');
legend("S2", "S4", "S5", "S8", "S11", "S12", "S13", "S15", "S16");

%% smoothing everything with moving average
 
% mm_CSER_W = cell(1,9); mm_CSER_S = cell(1,9);
mm_Syn_W = cell(1,9); mm_Syn_S = cell(1,9);

for i = 1:9
    mm_Syn_W{i} = movmean(fin_Syn_W{i},7);
    mm_Syn_S{i} = movmean(fin_Syn_S{i},7);
end

% plotting
% Awake
fig = figure();
hold on;
for i = 1:9
    plot(fin_CSER_W{i}, mm_Syn_W{i}, 'LineWidth', 1);
end
title('Humans awake','FontSize',15,'interpreter','latex');
xlabel('CSER');
ylabel('Synergy');
legend("S2", "S4", "S5", "S8", "S11", "S12", "S13", "S15", "S16");

% Sedated
fig = figure();
hold on;
for i = 1:9
    plot(fin_CSER_S{i}, mm_Syn_S{i}, 'LineWidth', 1);
end
title('Humans asleep','FontSize',15,'interpreter','latex');
xlabel('CSER');
ylabel('Synergy');
legend("S2", "S4", "S5", "S8", "S11", "S12", "S13", "S15", "S16");
save = "SynvsCSER_H.png";
exportgraphics(fig,save,'Resolution',300);
