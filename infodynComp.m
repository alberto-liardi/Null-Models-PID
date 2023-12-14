%% CSER - TOTAL CORRELATION - O-INFORMATION

% as = "Results/Kin2_c20/Sleep/";
% aw = "Results/Kin2_c20/Awake/";
% aw = "Results/S2_c20/Awake/";
% as = "Results/S2_c20/Sleep/";
% path = as+"/../";
% subject = "Kin2";

chan = '10';
drug = "LSD";
subject = "5";
aw = "Results_psychedelics/"+drug+"/"+subject+"_c"+chan+"/"+drug+"/";
as = "Results_psychedelics/"+drug+"/"+subject+"_c"+chan+"/Placebo/";
path = as+"../";

summary = readmatrix(as+"MMI/summary.txt");
N_tot = summary(1,4);
hun = round(log10(N_tot)-0.5,0);
N_tot = 10^hun;

InfoDyn_S = readmatrix(as+"InfoDyn.txt");
InfoDyn_W = readmatrix(aw+"InfoDyn.txt");

% CSER
CSER_S = InfoDyn_S(:,1:2); 
CSER_W = InfoDyn_W(:,1:2); 

m = min(min(CSER_S(:,1)), min(CSER_W(:,1)));
M = max(max(CSER_S(:,1)), max(CSER_W(:,1)));
edges1=linspace(m, M, 50);

[N_S,~] = histcounts(CSER_S(:,1), edges1);
[N_W,~] = histcounts(CSER_W(:,1), edges1);
N = round(max(max(N_S), max(N_W))+5*10^(hun-2), -(hun-1));

fig = figure();
mu_S = mean(CSER_S(:,1)); mu_W = mean(CSER_W(:,1));
subplot(2,1,1)
histogram(CSER_S(:,1), 'BinEdges', edges1, "FaceColor", "r");
ylim([0 N]);
labels = num2cell([0:N_tot/20:N]); labels{1} = [];
set(gca,'xtick',[],'YTick',[0:N_tot/20:N],'YTickLabel',labels);
xline(mu_S, 'Color', '#D95319', 'LineWidth', 2);
title("CSER "+subject+" Var(1)");
legend("Sedated","\mu= "+sprintf('%.2f',mu_S),'Location','northwest');
subplot(2,1,2)
histogram(CSER_W(:,1), 'BinEdges',edges1, "FaceColor", "g");
ylim([0 N]);
set(gca,'YTick',[0:N_tot/20:N]);
ha=get(gcf,'children');
set(ha(1),'position',[.1 .1 .8 .4]);
set(ha(3),'position',[.1 .5 .8 .4]);
xline(mu_W, 'Color', '#77AC30', 'LineWidth', 2);
legend("Awake", "\mu= "+sprintf('%.2f',mu_W),'Location','northwest');
exportgraphics(fig,path+"/CSER_Var(1).pdf",'Resolution',300);
% saveas(fig, path+"/CSER_Var(1).png");

m = min(min(CSER_S(:,2)), min(CSER_W(:,2)));
M = max(max(CSER_S(:,2)), max(CSER_W(:,2)));
edges2=linspace(m, M, 50);

[N_S,~] = histcounts(CSER_S(:,2), edges2);
[N_W,~] = histcounts(CSER_W(:,2), edges2);
N = round(max(max(N_S), max(N_W))+5*10^(hun-2), -(hun-1));

fig = figure();
mu_S = mean(CSER_S(:,2)); mu_W = mean(CSER_W(:,2));
subplot(2,1,1)
histogram(CSER_S(:,2), 'BinEdges', edges2, "FaceColor", '#EDB120', "EdgeColor", 'none');
ylim([0 N]);
labels = num2cell([0:N_tot/20:N]); labels{1} = [];
set(gca,'xtick',[],'YTick',[0:N_tot/20:N],'YTickLabel',labels);
xline(mu_S, '--', 'Color', 'black', 'LineWidth', 1); %#D95319
title("CSER "+subject+" Var(p)");
legend("Sedated","\mu= "+sprintf('%.2f',mu_S),'Location','northwest');
subplot(2,1,2)
histogram(CSER_W(:,2), 'BinEdges',edges2, "FaceColor", '#77AC30', "EdgeColor", 'none');
ylim([0 N]);
set(gca,'YTick',[0:N_tot/20:N]);
ha=get(gcf,'children');
set(ha(1),'position',[.1 .1 .8 .4]);
set(ha(3),'position',[.1 .5 .8 .4]);
xline(mu_W, '--', 'Color', 'black', 'LineWidth', 1); %'#77AC30'
legend("Awake", "\mu= "+sprintf('%.2f',mu_W),'Location','northwest');
exportgraphics(fig,path+"/CSER_Var(p).pdf",'Resolution',300);
% saveas(fig, path+"/CSER_Var(p).png");

%% CSER both in 1 plot

fig = figure('Position', [100 100 750 500]);
mu_S = mean(CSER_S(:,1)); mu_W = mean(CSER_W(:,1));
histogram(CSER_S(:,1), 'BinEdges', edges1, "FaceColor", '#EDB120', "EdgeColor", 'none');
xline(mu_S, '--', 'Color', 'black', 'LineWidth', 1.5); %#D95319
hold on;
histogram(CSER_W(:,1), 'BinEdges',edges1, "FaceColor", '#77AC30', "EdgeColor", 'none');
xline(mu_W, '--', 'Color', 'black', 'LineWidth', 1.5); %'#77AC30'
ylim([0 17]);
xlim([-30 5]);
xlabel("CSER", 'FontSize',15,'interpreter','latex');
ylabel("Counts");
set(gca,'FontName','CMU serif','FontSize',15);
legend("Sedated","\mu_S= "+sprintf('%.2f',mu_S),"Awake", "\mu_W= "+sprintf('%.2f',mu_W),'Location','northwest');
exportgraphics(fig,path+"/CSER_Var(1).pdf",'Resolution',300);

fig = figure('Position', [100 100 750 500]);
mu_S = mean(CSER_S(:,2)); mu_W = mean(CSER_W(:,2));
histogram(CSER_S(:,2), 'BinEdges', edges2, "FaceColor", '#EDB120', "EdgeColor", 'none');
xline(mu_S, '--', 'Color', 'black', 'LineWidth', 1.5); %#D95319
hold on;
histogram(CSER_W(:,2), 'BinEdges',edges2, "FaceColor", '#77AC30', "EdgeColor", 'none');
xline(mu_W, '--', 'Color', 'black', 'LineWidth', 1.5); %'#77AC30'
ylim([0 17]);
xlim([-30 0]);
xlabel("CSER", 'FontSize',15,'interpreter','latex');
ylabel("Counts");
set(gca,'FontName','CMU serif','FontSize',15);
legend("Sedated","\mu_S= "+sprintf('%.2f',mu_S),"Awake", "\mu_W= "+sprintf('%.2f',mu_W),'Location','northwest');
exportgraphics(fig,path+"/CSER_Var(p).pdf",'Resolution',300);

%%
% TOTAL CORRELATION
TC_S = InfoDyn_S(:,3:4);
TC_W = InfoDyn_W(:,3:4);

m = min(min(TC_S(:,1)), min(TC_W(:,1)));
M = max(max(TC_S(:,1)), max(TC_W(:,1)));
edges1=linspace(m, M, 50);

[N_S,~] = histcounts(TC_S(:,1), edges1);
[N_W,~] = histcounts(TC_W(:,1), edges1);
N = round(max(max(N_S), max(N_W))+5*10^(hun-2), -(hun-1));

fig = figure();
mu_S = mean(TC_S(:,1)); mu_W = mean(TC_W(:,1));
subplot(2,1,1)
histogram(TC_S(:,1), 'BinEdges', edges1, "FaceColor", "r");
ylim([0 N]);
labels = num2cell([0:N_tot/20:N]); labels{1} = [];
set(gca,'xtick',[],'YTick',[0:N_tot/20:N],'YTickLabel',labels);
xline(mu_S, 'Color', '#D95319', 'LineWidth', 2);
title("Total Correlation "+subject+" Var(1)");
legend("Sedated","\mu= "+sprintf('%.2f',mu_S),'Location','northwest');
subplot(2,1,2)
histogram(TC_W(:,1), 'BinEdges',edges1, "FaceColor", "g");
ylim([0 N]);
set(gca,'YTick',[0:N_tot/20:N]);
ha=get(gcf,'children');
set(ha(1),'position',[.1 .1 .8 .4]);
set(ha(3),'position',[.1 .5 .8 .4]);
xline(mu_W, 'Color', '#77AC30', 'LineWidth', 2);
legend("Awake", "\mu= "+sprintf('%.2f',mu_W),'Location','northwest');
saveas(fig, path+"/TC_Var(1).png");

m = min(min(TC_S(:,2)), min(TC_W(:,2)));
M = max(max(TC_S(:,2)), max(TC_W(:,2)));
edges2=linspace(m, M, 50);

[N_S,~] = histcounts(TC_S(:,2), edges2);
[N_W,~] = histcounts(TC_W(:,2), edges2);
N = round(max(max(N_S), max(N_W))+5*10^(hun-2), -(hun-1));

fig = figure();
mu_S = mean(TC_S(:,2)); mu_W = mean(TC_W(:,2));
subplot(2,1,1)
histogram(TC_S(:,2), 'BinEdges', edges2, "FaceColor", "r");
ylim([0 N]);
labels = num2cell([0:N_tot/20:N]); labels{1} = [];
set(gca,'xtick',[],'YTick',[0:N_tot/20:N],'YTickLabel',labels);
xline(mu_S, 'Color', '#D95319', 'LineWidth', 2);
title("Total Correlation "+subject+" Var(p)");
legend("Sedated","\mu= "+sprintf('%.2f',mu_S),'Location','northwest');
subplot(2,1,2)
histogram(TC_W(:,2), 'BinEdges',edges2, "FaceColor", "g");
ylim([0 N]);
set(gca,'YTick',[0:N_tot/20:N]);
ha=get(gcf,'children');
set(ha(1),'position',[.1 .1 .8 .4]);
set(ha(3),'position',[.1 .5 .8 .4]);
xline(mu_W, 'Color', '#77AC30', 'LineWidth', 2);
legend("Awake", "\mu= "+sprintf('%.2f',mu_W),'Location','northwest');
saveas(fig, path+"/TC_Var(p).png");


%O-INFORMATION
Oin_S = InfoDyn_S(:,5:6);
Oin_W = InfoDyn_W(:,5:6);

m = min(min(Oin_S(:,1)), min(Oin_W(:,1)));
M = max(max(Oin_S(:,1)), max(Oin_W(:,1)));
edges1=linspace(m, M, 50);

[N_S,~] = histcounts(Oin_S(:,1), edges1);
[N_W,~] = histcounts(Oin_W(:,1), edges1);
N = round(max(max(N_S), max(N_W))+5*10^(hun-2), -(hun-1));

fig = figure();
mu_S = mean(Oin_S(:,1)); mu_W = mean(Oin_W(:,1));
subplot(2,1,1)
histogram(Oin_S(:,1), 'BinEdges', edges1, "FaceColor", "r");
ylim([0 N]);
labels = num2cell([0:N_tot/20:N]); labels{1} = [];
set(gca,'xtick',[],'YTick',[0:N_tot/20:N],'YTickLabel',labels);
xline(mu_S, 'Color', '#D95319', 'LineWidth', 2);
title("O Information "+subject+" Var(1)");
legend("Sedated","\mu= "+sprintf('%.2f',mu_S),'Location','northwest');
subplot(2,1,2)
histogram(Oin_W(:,1), 'BinEdges',edges1, "FaceColor", "g");
ylim([0 N]);
set(gca,'YTick',[0:N_tot/20:N]);
ha=get(gcf,'children');
set(ha(1),'position',[.1 .1 .8 .4]);
set(ha(3),'position',[.1 .5 .8 .4]);
xline(mu_W, 'Color', '#77AC30', 'LineWidth', 2);
legend("Awake", "\mu= "+sprintf('%.2f',mu_W),'Location','northwest');
exportgraphics(fig,path+"/Oinfo_Var(1).pdf",'Resolution',300);
% saveas(fig, path+"/Oinfo_Var(1).png");

m = min(min(Oin_S(:,2)), min(Oin_W(:,2)));
M = max(max(Oin_S(:,2)), max(Oin_W(:,2)));
edges2=linspace(m, M, 50);

[N_S,~] = histcounts(Oin_S(:,2), edges2);
[N_W,~] = histcounts(Oin_W(:,2), edges2);
N = round(max(max(N_S), max(N_W))+5*10^(hun-2), -(hun-1));

fig = figure();
mu_S = mean(Oin_S(:,2)); mu_W = mean(Oin_W(:,2));
subplot(2,1,1)
histogram(Oin_S(:,2), 'BinEdges', edges2, "FaceColor", "r");
ylim([0 N]);
labels = num2cell([0:N_tot/20:N]); labels{1} = [];
set(gca,'xtick',[],'YTick',[0:N_tot/20:N],'YTickLabel',labels);
xline(mu_S, 'Color', '#D95319', 'LineWidth', 2);
title("O Information "+subject+" Var(p)");
legend("Sedated","\mu= "+sprintf('%.2f',mu_S),'Location','northwest');
subplot(2,1,2)
histogram(Oin_W(:,2), 'BinEdges',edges2, "FaceColor", "g");
ylim([0 N]);
set(gca,'YTick',[0:N_tot/20:N]);
ha=get(gcf,'children');
set(ha(1),'position',[.1 .1 .8 .4]);
set(ha(3),'position',[.1 .5 .8 .4]);
xline(mu_W, 'Color', '#77AC30', 'LineWidth', 2);
legend("Awake", "\mu= "+sprintf('%.2f',mu_W),'Location','northwest');
exportgraphics(fig,path+"/Oinfo_Var(p).pdf",'Resolution',300);
saveas(fig, path+"/Oinfo_Var(p).png");

%% Model Order distribution

% as = "Results/Chibi_c2/Sedated/";
% aw = "Results/Chibi_c2/Awake/";
% aw = "Results/S2_c20/Awake/";
% as = "Results/S2_c20/Sleep/";
chan = '2';
drug = "PSIL";
% subjects = "1";

for subject = ["1", "2", "3", "4", "5", "6"]
aw = "Results_psychedelics/"+drug+"/"+subject+"_c"+chan+"/"+drug+"/";
as = "Results_psychedelics/"+drug+"/"+subject+"_c"+chan+"/Placebo/";
path = as+"../";


summary = readmatrix(as+"MMI/summary.txt");
N_tot = summary(1,4);
hun = round(log10(N_tot)-0.5,0);
N_tot = 10^hun;

InfoDyn_S = readmatrix(as+"InfoDyn.txt");
InfoDyn_W = readmatrix(aw+"InfoDyn.txt");

% model order history
pHisto_S = InfoDyn_S(:,7:8);
pHisto_W = InfoDyn_W(:,7:8);

m = min(min(pHisto_S(:,2)), min(pHisto_W(:,2)));
M = max(max(pHisto_S(:,2)), max(pHisto_W(:,2)));
edges1=linspace(m, M, 50);

[N_S,~] = histcounts(pHisto_S(:,2), edges1);
[N_W,~] = histcounts(pHisto_W(:,2), edges1);
N = round(max(max(N_S), max(N_W))+5*10^(hun-2), -(hun-1));

fig = figure();
mu_S = mean(pHisto_S(:,2)); mu_W = mean(pHisto_W(:,2));
subplot(2,1,1)
histogram(pHisto_S(:,2), 'BinEdges', edges1, "FaceColor", "r");
ylim([0 N]);
labels = num2cell([0:N_tot/20:N]); labels{1} = [];
set(gca,'xtick',[],'YTick',[0:N_tot/20:N],'YTickLabel',labels);
xline(mu_S, 'Color', '#D95319', 'LineWidth', 2);
title("Model order history subject "+subject+" Var(p)");
legend("Sedated","\mu= "+sprintf('%.2f',mu_S),'Location','northwest');
subplot(2,1,2)
histogram(pHisto_W(:,2), 'BinEdges',edges1, "FaceColor", "g");
ylim([0 N]);
set(gca,'YTick',[0:N_tot/20:N]);
ha=get(gcf,'children');
set(ha(1),'position',[.1 .1 .8 .4]);
set(ha(3),'position',[.1 .5 .8 .4]);
xline(mu_W, 'Color', '#77AC30', 'LineWidth', 2);
legend("Awake", "\mu= "+sprintf('%.2f',mu_W),'Location','northwest');
saveas(fig, path+"/pHisto_Var(p)_filt.png");

end