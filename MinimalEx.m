rng('default') % for reproducibility of random variables
javaaddpath("infodynamics.jar");

A = RandEvolMatrix(1,1);

% A = [0.8147, 0.0126; 1.2699, 0.9058];


% (k,g, PID) 3D plot

path = "Results/Minimal_Example_Theor/";

% try to load pre-saved data, otherwise do the calculation

try
    k = load(path+"k", 'k').k;
    g = load(path+"g", 'g').g;
    Unique_X = load(path+"Unique_X", "Unique_X").Unique_X;
    Unique_Y = load(path+"Unique_Y", "Unique_Y").Unique_Y;
    Synergy = load(path+"Synergy", "Synergy").Synergy;
    Redundancy = load(path+"Redundancy", "Redundancy").Redundancy;
    MI = load(path+"MI", "MI").MI;
    CSER = load(path+"CSER", "CSER").CSER;
    corr = load(path+"Correlation", "corr").corr;
    Oinfo = load(path+'Oinfo','Oinfo').Oinfo;
    Sinfo = load(path+'Sinfo','Sinfo').Sinfo;

catch
    k = -0.99:0.01:0.99;
    g = 1:0.01:3;
    
    Unique_X = zeros(3,length(g),length(k));
    Unique_Y = zeros(3,length(g),length(k));
    Synergy = zeros(3,length(g),length(k));
    Redundancy = zeros(3, length(g),length(k));
    MI = zeros(3, length(g),length(k));
    CSER = zeros(length(g),length(k));
    corr = zeros(length(g),length(k));

    Oinfo = zeros(length(g),length(k));
    Sinfo = zeros(length(g),length(k));
    oCalc = infodynamics.measures.continuous.gaussian.OInfoCalculatorGaussian();
    sCalc = infodynamics.measures.continuous.gaussian.SInfoCalculatorGaussian();
    
    for n = 1:length(g)
        for l = 1:length(k)
            disp(n); disp(l);

            V = [g(n) k(l); k(l) g(n)];
            
            PID_atoms = PIDcalc_Mult(1,A,V,1,1,"joint",['y', 'y', 'y']);
            
            Unique_X(:,n,l) = PID_atoms(:,1);
            Unique_Y(:,n,l) = PID_atoms(:,2);
            Synergy(:,n,l) = PID_atoms(:,4);
            Redundancy(:,n,l) = PID_atoms(:,3);
            MI(:,n,l) = sum(PID_atoms,2);
            CSER(n,l) = .5*log(det(2*pi*exp(1)*V));
            temp = diag(diag(V).^(-1/2))*V*diag(diag(V).^(-1/2));
            corr(n,l) = temp(1,2);
            
            oCalc.initialise(2);     oCalc.setCovariance(V);
            Oinfo(n,l) = oCalc.computeAverageLocalOfObservations();
            sCalc.initialise(2);     sCalc.setCovariance(V);
            Sinfo(n,l) = sCalc.computeAverageLocalOfObservations();
        end
    end

    % save results
    save(path+"k", 'k');
    save(path+"g", 'g');
    save(path+'Unique_X.mat','Unique_X');
    save(path+'Unique_Y.mat','Unique_Y');
    save(path+'Synergy.mat','Synergy');
    save(path+'Redundancy.mat','Redundancy');
    save(path+'MI.mat','MI');
    save(path+'CSER.mat','CSER');
    save(path+'Correlation.mat','corr');
    save(path+'Oinfo.mat','Oinfo');
    save(path+'Sinfo.mat','Sinfo');

end

disp("done");

%% HEATMAPS

% heatmaps for PID
% heatmaps(k,g,Unique_X,Unique_Y,Synergy,Redundancy,"abs",CSER,corr,MI);

% normalise Synergy and Redundancy
Synergy_norm = Synergy./MI;
Redundancy_norm = Redundancy./MI;
Unique_X_norm = Unique_X./MI;
Unique_Y_norm = Unique_Y./MI;

heatmaps(k,g,Unique_X_norm,Unique_Y_norm,Synergy_norm,Redundancy_norm,"norm");

% heatmaps for Information Dynamics (O information)
% heatmaps_dyn(k,g,Oinfo,Sinfo);


%% 3D plots
% Synergy = Synergy./MI;
% Redundancy = Redundancy./MI;

% O-info, S-info, and O-info normalised
[X,Y] = meshgrid(k,g);

fig = figure();
surf(X,Y,Oinfo);
xlabel("k"); ylabel("g");
title("O Information");
fig = figure();
surf(X,Y,Sinfo);
title("S Information");
xlabel("k"); ylabel("g");

Oinfo_norm = Oinfo ./ Sinfo;
fig = figure();
surf(X,Y,Oinfo_norm);
title("O Information Normalised");
xlabel("k"); ylabel("g");

% TC and DRM
TC = (Oinfo + Sinfo) ./ 2;
DTC = (- Oinfo + Sinfo) ./ 2;

[X,Y] = meshgrid(k,g);

fig = figure();
surf(X,Y,TC);
xlabel("k"); ylabel("g");
title("Total correlation");
fig = figure();
surf(X,Y,DTC);
title("Shared randomness");
xlabel("k"); ylabel("g");

% MI
MI_pl = squeeze(MI(1,:,:));

fig = figure();
surf(X,Y,MI_pl);
title("Mutual Information");
xlabel("k"); ylabel("g");


% CSER and Correlation
[X,Y] = meshgrid(k,g);

fig = figure();
surf(X,Y,CSER);
xlabel("k"); ylabel("g");
title("CSER");
fig = figure();
surf(X,Y,corr);
title("Residual Correlation");
xlabel("k"); ylabel("g");


% MMI
S = squeeze(Synergy(1,:,:));
R = squeeze(Redundancy(1,:,:));

fig = figure();
surf(X,Y,S);
title("MMI - Synergy");
xlabel("k"); ylabel("g");
fig = figure();
surf(X,Y,R);
xlabel("k"); ylabel("g");
title("MMI - Redundancy");


% IDEP
[X,Y] = meshgrid(k,g);
S = squeeze(Synergy(2,:,:));
R = squeeze(Redundancy(2,:,:));

fig = figure();
surf(X,Y,S);
title("IDEP- Synergy");
xlabel("k"); ylabel("g");
fig = figure();
surf(X,Y,R);
xlabel("k"); ylabel("g");
title("IDEP - Redundancy");

% ICCS
[X,Y] = meshgrid(k,g);
S = squeeze(Synergy(3,:,:));
R = squeeze(Redundancy(3,:,:));

fig = figure();
surf(X,Y,S);
title("ICCS - Synergy");
xlabel("k"); ylabel("g");
fig = figure();
surf(X,Y,R);
xlabel("k"); ylabel("g");
title("ICCS - Redundancy");


%% GRADIENT PLOTS

% MMI
[DX_C,DY_C] = gradient(CSER(1:5:end,1:5:end),0.02*5,0.02*5);
[DX_S,DY_S] = gradient(S(1:5:end,1:5:end),0.02*5,0.02*5);
[DX_R,DY_R] = gradient(R(1:5:end,1:5:end),0.02*5,0.02*5);
fig = figure();
title("MMI");
quiver(k(1:5:end),g(1:5:end),DX_C,DY_C,'LineWidth',1);
hold on
quiver(k(1:5:end),g(1:5:end),DX_S,DY_S,'LineWidth',1);
hold on
quiver(k(1:5:end),g(1:5:end),DX_R,DY_R,'LineWidth',1);
axis equal
xlabel("k"); ylabel("g");


%% REAL ANALYSES (BASED ON MONKEY KIN2)- V sweep with A_W and A_S

folder = "Results/Minimal_Example_Real_Kin2_Vsweep/";

Cov_W = load(folder+"Cov_W", 'Cov_W').Cov_W;
Cov_S = load(folder+"Cov_S", 'Cov_S').Cov_S;

A_W = load(folder+"A_W", 'A_W').A_W;
A_S = load(folder+"A_S", 'A_S').A_S;
A = {A_S A_W};

L = 100;
lambdas = 0:1/L:1;

for a = 1:2
    if a==1, path=folder+"S/";
    else,   path=folder+"W/";  end

    try
        Unique_X = load(path+"Unique_X", "Unique_X").Unique_X;
        Unique_Y = load(path+"Unique_Y", "Unique_Y").Unique_Y;
        Synergy = load(path+"Synergy", "Synergy").Synergy;
        Redundancy = load(path+"Redundancy", "Redundancy").Redundancy;
        MI = load(path+"MI", "MI").MI;
        CSER = load(path+"CSER", "CSER").CSER;
        corr = load(path+"Correlation", "corr").corr;
        Oinfo = load(path+'Oinfo','Oinfo').Oinfo;
        Sinfo = load(path+'Sinfo','Sinfo').Sinfo;
    
    catch
        Unique_X = zeros(3,L+1);
        Unique_Y = zeros(3,L+1);
        Synergy = zeros(3,L+1);
        Redundancy = zeros(3,L+1);
    
        MI = zeros(1,L+1);
        CSER = zeros(1,L+1);
        corr = zeros(1,L+1);
        Oinfo = zeros(1,L+1);
        Sinfo = zeros(1,L+1);
        oCalc = infodynamics.measures.continuous.gaussian.OInfoCalculatorGaussian();
        sCalc = infodynamics.measures.continuous.gaussian.SInfoCalculatorGaussian();
       
        % going from S to W, varying lambda
        for l = 1:length(lambdas)
            V = lambdas(l)*Cov_W + (1-lambdas(l))*Cov_S;
        
            PID_atoms = PIDcalc_Mult(1,{A{a}},V,1,1,"joint",['y', 'y', 'y']);
            
            Unique_X(:,l) = PID_atoms(:,1);
            Unique_Y(:,l) = PID_atoms(:,2);
            Synergy(:,l) = PID_atoms(:,4);
            Redundancy(:,l) = PID_atoms(:,3);
            MI(l) = sum(PID_atoms(1,:));
            CSER(l) = .5*log(det(2*pi*exp(1)*V));
            temp = diag(diag(V).^(-1/2))*V*diag(diag(V).^(-1/2));
            corr(l) = temp(1,2);
            
            oCalc.initialise(2);     oCalc.setCovariance(V);
            Oinfo(l) = oCalc.computeAverageLocalOfObservations();
            sCalc.initialise(2);     sCalc.setCovariance(V);
            Sinfo(l) = sCalc.computeAverageLocalOfObservations();
        
        end
        
        save(path+'Unique_X.mat','Unique_X');
        save(path+'Unique_Y.mat','Unique_Y');
        save(path+'Synergy.mat','Synergy');
        save(path+'Redundancy.mat','Redundancy');
        save(path+'MI.mat','MI');
        save(path+'CSER.mat','CSER');
        save(path+'Correlation.mat','corr');
        save(path+'Oinfo.mat','Oinfo');
        save(path+'Sinfo.mat','Sinfo');
    end
    
    
    % draw plots (those are 2D, only 1 parameter lambda!)
    if a==1, A_spec=" with A from Sleep/Sedated ";
    else,   A_spec=" with A from Awake ";  end
    
    x = lambdas;
    
    % CSER and co.
    fig = figure();
    plot(x,CSER,x,corr,x,MI,'LineWidth',2);
    set(gca,'FontName','CMU serif','FontSize',15);
    legend("CSER", "Correlation", "MI",'Interpreter','latex','Location','southeast');
%     title("CSER, correlation and MI "+V_spec);
    ylabel("Conscious - Unconscious",'FontSize',15,'interpreter','latex');
    xlabel("$\lambda$",'FontSize',15,'interpreter','latex');
    ha=get(gcf,'children');
    set(ha(1),'position',[.7 .5 .1 .1]);
    save = path+"CSERandCo.png";
    exportgraphics(fig,save,'Resolution',300)
    % MMI
    fig = figure();
    plot(x,Unique_X(1,:),x,Unique_Y(1,:),x,Synergy(1,:),x,Redundancy(1,:),'LineWidth',2);
    set(gca,'FontName','CMU serif','FontSize',15);
%     title("MMI"+V_spec);
    ylabel("Conscious - Unconscious",'FontSize',15,'interpreter','latex');
    xlabel("$\lambda$",'FontSize',15,'interpreter','latex');
    legend("Unique X", "Unique Y", "Synergy", "Redundancy",'Interpreter','latex', ...
        'Location','southeast');
    ha=get(gcf,'children');
    set(ha(1),'position',[.7 .5 .1 .1]);
    save = path+"MMI.png";
    exportgraphics(fig,save,'Resolution',300)
%     % Idep
%     figure();
%     plot(x,Unique_X(2,:),x,Unique_Y(2,:),x,Synergy(2,:),x,Redundancy(2,:));
%     title("Idep"+A_spec);
%     legend("Unique X", "Unique Y", "Synergy", "Redundancy");
%     % Iccs
%     figure();
%     plot(x,Unique_X(3,:),x,Unique_Y(3,:),x,Synergy(3,:),x,Redundancy(3,:));
%     title("Iccs"+A_spec);
%     legend("Unique X", "Unique Y", "Synergy", "Redundancy");
end



%% REAL ANALYSES (BASED ON MONKEY KIN2)- A sweep with fixed V_W or V_S

folder = "Results/Minimal_Example_Real_Kin2_Asweep/";

Cov_W = load(folder+"Cov_W", 'Cov_W').Cov_W;
Cov_S = load(folder+"Cov_S", 'Cov_S').Cov_S;
V = {Cov_W, Cov_S};

A_W = load(folder+"A_W", 'A_W').A_W;
A_S = load(folder+"A_S", 'A_S').A_S;

L = 100;
lambdas = 0:1/L:1;

for v = 1:2
    if v==1, path=folder+"S/";
    else,   path=folder+"W/";  end

    try
        Unique_X = load(path+"Unique_X", "Unique_X").Unique_X;
        Unique_Y = load(path+"Unique_Y", "Unique_Y").Unique_Y;
        Synergy = load(path+"Synergy", "Synergy").Synergy;
        Redundancy = load(path+"Redundancy", "Redundancy").Redundancy;
        MI = load(path+"MI", "MI").MI;
        CSER = load(path+"CSER", "CSER").CSER;
        corr = load(path+"Correlation", "corr").corr;
        Oinfo = load(path+'Oinfo','Oinfo').Oinfo;
        Sinfo = load(path+'Sinfo','Sinfo').Sinfo;
    
    catch
        Unique_X = zeros(3,L+1);
        Unique_Y = zeros(3,L+1);
        Synergy = zeros(3,L+1);
        Redundancy = zeros(3,L+1);
    
        MI = zeros(1,L+1);
        CSER = zeros(1,L+1);
        corr = zeros(1,L+1);
        Oinfo = zeros(1,L+1);
        Sinfo = zeros(1,L+1);
        oCalc = infodynamics.measures.continuous.gaussian.OInfoCalculatorGaussian();
        sCalc = infodynamics.measures.continuous.gaussian.SInfoCalculatorGaussian();
       
        % going from S to W, varying lambda
        for l = 1:length(lambdas)
            A = lambdas(l)*A_W + (1-lambdas(l))*A_S;
        
            PID_atoms = PIDcalc_Mult(1,{A},V{v},1,1,"joint",['y', 'y', 'y']);
            
            Unique_X(:,l) = PID_atoms(:,1);
            Unique_Y(:,l) = PID_atoms(:,2);
            Synergy(:,l) = PID_atoms(:,4);
            Redundancy(:,l) = PID_atoms(:,3);
            MI(l) = sum(PID_atoms(1,:));
            CSER(l) = .5*log(det(2*pi*exp(1)*V{v}));
            temp = diag(diag(V{v}).^(-1/2))*V{v}*diag(diag(V{v}).^(-1/2));
            corr(l) = temp(1,2);
            
            oCalc.initialise(2);     oCalc.setCovariance(V{v});
            Oinfo(l) = oCalc.computeAverageLocalOfObservations();
            sCalc.initialise(2);     sCalc.setCovariance(V{v});
            Sinfo(l) = sCalc.computeAverageLocalOfObservations();
        
        end
        
        save(path+'Unique_X.mat','Unique_X');
        save(path+'Unique_Y.mat','Unique_Y');
        save(path+'Synergy.mat','Synergy');
        save(path+'Redundancy.mat','Redundancy');
        save(path+'MI.mat','MI');
        save(path+'CSER.mat','CSER');
        save(path+'Correlation.mat','corr');
        save(path+'Oinfo.mat','Oinfo');
        save(path+'Sinfo.mat','Sinfo');
    end
    
    
    % draw plots (those are 2D, only 1 parameter lambda!)
    if v==1, V_spec=" with V from Sleep/Sedated ";
    else,   V_spec=" with V from Awake ";  end
    
    x = lambdas;
    
    % CSER and co.
    fig = figure();
    plot(x,CSER,x,corr,x,MI,'LineWidth',2);
    set(gca,'FontName','CMU serif','FontSize',15);
    legend("CSER", "Correlation", "MI",'Interpreter','latex');
%     title("CSER, correlation and MI "+V_spec);
    ylabel("Conscious - Unconscious",'FontSize',15,'interpreter','latex');
    xlabel("$\lambda$",'FontSize',15,'interpreter','latex');
    save = path+"CSERandCo.png";
    exportgraphics(fig,save,'Resolution',300)
    % MMI
    fig = figure();
    plot(x,Unique_X(1,:),x,Unique_Y(1,:),x,Synergy(1,:),x,Redundancy(1,:),'LineWidth',2);
    set(gca,'FontName','CMU serif','FontSize',15);
%     title("MMI"+V_spec);
    ylabel("Conscious - Unconscious",'FontSize',15,'interpreter','latex');
    xlabel("$\lambda$",'FontSize',15,'interpreter','latex');
    legend("Unique X", "Unique Y", "Synergy", "Redundancy",'Interpreter','latex');
    save = path+"MMI.png";
    exportgraphics(fig,save,'Resolution',300)
%     % Idep
%     figure();
%     plot(x,Unique_X(2,:),x,Unique_Y(2,:),x,Synergy(2,:),x,Redundancy(2,:));
%     title("Idep"+V_spec);
%     legend("Unique X", "Unique Y", "Synergy", "Redundancy");
%     % Iccs
%     figure();
%     plot(x,Unique_X(3,:),x,Unique_Y(3,:),x,Synergy(3,:),x,Redundancy(3,:));
%     title("Iccs"+V_spec);
%     legend("Unique X", "Unique Y", "Synergy", "Redundancy");
end

function [] = heatmaps(k,g,Unique_X,Unique_Y,Synergy,Redundancy,type,CSER,corr,MI)
    path = "Results/Minimal_Example_Theor_temp/";

    kkk = [k(1), k(5:5:end-4),k(end)];
    xlabel = string(kkk);
    kind = 1:4:length(kkk);
    xlabelfin = strings(1,length(kkk));
    xlabelfin(kind) = xlabel(kind);
    dp = compose("%.2f", xlabelfin);
    xlabelfin(kind) = dp(kind);

    ggg = flip(g(1:5:end));
    ylabel = string(ggg);
    gind = 1:4:length(ggg);
    ylabelfin = strings(1,length(ggg));
    ylabelfin(gind) = ylabel(gind);
    dp = compose("%.2f", ylabelfin);
    ylabelfin(gind) = dp(gind);


    if exist('CSER','var')
        % CSER, Correlation and MI
        fig = figure('Position', [100 100 550 500]);
        a = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
                    flip([CSER(1:5:end,1), CSER(1:5:end,5:5:end-4), ...
                    CSER(1:5:end,end)]), 'Colormap',jet);
%         a.Title = 'CSER';
        grid off;
        a.XLabel = 'k';
        a.YLabel = 'g';
        save = path+'CSER.pdf';        
        a.XDisplayLabels = xlabelfin;        
        a.YDisplayLabels = ylabelfin;
        set(gca,'FontName','CMU serif','FontSize',15);
        h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
        h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
        exportgraphics(fig,save,'Resolution',300);

        fig = figure('Position', [100 100 550 500]);
        d = heatmap(kkk,ggg, ... 
                    flip([corr(1:5:end,1), corr(1:5:end,5:5:end-4), ...
                    corr(1:5:end,end)]), 'Colormap',jet);
%         d.Title = 'Correlation';
        grid off;
        d.XLabel = 'k';
        d.YLabel = 'g';
        save = path+'Corr.pdf';
        d.XDisplayLabels = xlabelfin;        
        d.YDisplayLabels = ylabelfin;
        set(gca,'FontName','CMU serif','FontSize',15);
        h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
        h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
        exportgraphics(fig,save,'Resolution',300);

        fig = figure('Position', [100 100 550 500]);
        d = heatmap(kkk,ggg, ... 
                    flip([MI(1,1:5:end,1)', squeeze(MI(1,1:5:end,5:5:end-4)), ...
                    MI(1,1:5:end,end)']), 'Colormap',jet);
%         d.Title = 'Mutual Information';
        grid off;
        d.XLabel = 'k';
        d.YLabel = 'g';
        save = path+'MI.pdf';
        d.XDisplayLabels = xlabelfin;        
        d.YDisplayLabels = ylabelfin;
        set(gca,'FontName','CMU serif','FontSize',15);
        h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
        h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
        exportgraphics(fig,save,'Resolution',300);
    end

    if type=="abs"
        Syn_lab = "Synergy ";
        Red_lab = "Redundancy ";
        Un_x_lab = "Unique X ";
        Un_y_lab = "Unique Y ";
    elseif type=="norm"
        Syn_lab = "Normalised Synergy ";
        Red_lab = "Normalised Redundancy ";
        Un_x_lab = "Normalised Unique X ";
        Un_y_lab = "Normalised Unique Y ";
    end
    
    % MMI
    fig = figure('Position', [100 100 550 500]);
    b = heatmap(kkk,ggg, ... 
                flip([Unique_X(1,1:5:end,1)', squeeze(Unique_X(1,1:5:end,5:5:end-4)), ...
                Unique_X(1,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Un_x_lab+'MMI';
    grid off;
    b.XLabel = 'k';
    b.YLabel = 'g';
    save = path+Un_x_lab+'MMI.pdf';
    b.XDisplayLabels = xlabelfin;        
    b.YDisplayLabels = ylabelfin;
    set(gca,'FontName','CMU serif','FontSize',15);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
    exportgraphics(fig,save,'Resolution',300);

    fig = figure('Position', [100 100 550 500]);
    b = heatmap(kkk,ggg, ... 
                flip([Unique_Y(1,1:5:end,1)', squeeze(Unique_Y(1,1:5:end,5:5:end-4)), ...
                Unique_Y(1,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Un_y_lab+'MMI';
    grid off;
    b.XLabel = 'k';
    b.YLabel = 'g';
    save = path+Un_x_lab+'MMI.pdf';
    b.XDisplayLabels = xlabelfin;        
    b.YDisplayLabels = ylabelfin;
    set(gca,'FontName','CMU serif','FontSize',15);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
    exportgraphics(fig,save,'Resolution',300);

    fig = figure('Position', [100 100 550 500]);
    b = heatmap(kkk,ggg, ... 
                flip([Synergy(1,1:5:end,1)', squeeze(Synergy(1,1:5:end,5:5:end-4)), ...
                Synergy(1,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Syn_lab+'MMI';
    grid off;
    b.XLabel = 'k';
    b.YLabel = 'g';
    save = path+Syn_lab+'MMI.pdf';
    b.XDisplayLabels = xlabelfin;        
    b.YDisplayLabels = ylabelfin;
    set(gca,'FontName','CMU serif','FontSize',15);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
    exportgraphics(fig,save,'Resolution',300);
    
    fig = figure('Position', [100 100 550 500]);
    c = heatmap(kkk,ggg, ... 
                flip([Redundancy(1,1:5:end,1)', squeeze(Redundancy(1,1:5:end,5:5:end-4)), ...
                Redundancy(1,1:5:end,end)']), 'Colormap',jet);
%     c.Title = Red_lab+ 'MMI';
    grid off;
    c.XLabel = 'k';
    c.YLabel = 'g';
    save = path+Red_lab+'MMI.pdf';
    c.XDisplayLabels = xlabelfin;        
    c.YDisplayLabels = ylabelfin;
    set(gca,'FontName','CMU serif','FontSize',15);
    h.NodeChildren(3).XAxis.Label.Interpreter = 'latex';
    h.NodeChildren(3).YAxis.Label.Interpreter = 'latex';
    exportgraphics(fig,save,'Resolution',300);
    
    
%     % Idep
%     fig = figure('Position', [100 100 550 500]);
%     b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
%                 flip([Unique_X(2,1:5:end,1)', squeeze(Unique_X(2,1:5:end,5:5:end-4)), ...
%                 Unique_X(2,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Un_x_lab+'MMI';
%     b.XLabel = 'k';
%     b.YLabel = 'g';
%     fig = figure('Position', [100 100 550 500]);
%     b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
%                 flip([Unique_Y(2,1:5:end,1)', squeeze(Unique_Y(2,1:5:end,5:5:end-4)), ...
%                 Unique_Y(2,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Un_y_lab+'MMI';
%     b.XLabel = 'k';
%     b.YLabel = 'g';
%     fig = figure('Position', [100 100 550 500]);
%     b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
%                 flip([Synergy(2,1:5:end,1)', squeeze(Synergy(2,1:5:end,5:5:end-4)), ...
%                 Synergy(2,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Syn_lab+'Idep';
%     b.XLabel = 'k';
%     b.YLabel = 'g';
%     fig = figure('Position', [100 100 550 500]);
%     c = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
%                 flip([Redundancy(2,1:5:end,1)', squeeze(Redundancy(2,1:5:end,5:5:end-4)), ...
%                 Redundancy(2,1:5:end,end)']), 'Colormap',jet);
%     c.Title = Red_lab+ 'Idep';
%     c.XLabel = 'k';
%     c.YLabel = 'g';
%     
%     
%     % Iccs
%     fig = figure('Position', [100 100 550 500]);
%     b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
%                 flip([Unique_X(3,1:5:end,1)', squeeze(Unique_X(3,1:5:end,5:5:end-4)), ...
%                 Unique_X(3,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Un_x_lab+'MMI';
%     b.XLabel = 'k';
%     b.YLabel = 'g';
%     fig = figure('Position', [100 100 550 500]);
%     b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
%                 flip([Unique_Y(3,1:5:end,1)', squeeze(Unique_Y(3,1:5:end,5:5:end-4)), ...
%                 Unique_Y(3,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Un_y_lab+'MMI';
%     b.XLabel = 'k';
%     b.YLabel = 'g';
%     fig = figure('Position', [100 100 550 500]);
%     b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
%                 flip([Synergy(3,1:5:end,1)', squeeze(Synergy(3,1:5:end,5:5:end-4)), ...
%                 Synergy(3,1:5:end,end)']), 'Colormap',jet);
%     b.Title = Syn_lab+'Iccs';
%     b.XLabel = 'k';
%     b.YLabel = 'g';
%     fig = figure('Position', [100 100 550 500]);
%     c = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
%                 flip([Redundancy(3,1:5:end,1)', squeeze(Redundancy(3,1:5:end,5:5:end-4)), ...
%                 Redundancy(3,1:5:end,end)']), 'Colormap',jet);
%     c.Title = Red_lab+ 'Iccs';
%     c.XLabel = 'k';
%     c.YLabel = 'g';
end


function [] = heatmaps_dyn(k,g,Oinfo, Sinfo)

    % MMI
    fig = figure();
    b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
                flip([Oinfo(1:5:end,1), Oinfo(1:5:end,5:5:end-4), ...
                Oinfo(1:5:end,end)]), 'Colormap',jet);
    b.Title = ' O Information';
    b.XLabel = 'k';
    b.YLabel = 'g';
    fig = figure();
    b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
                flip([Sinfo(1:5:end,1), Sinfo(1:5:end,5:5:end-4), ...
                Sinfo(1:5:end,end)]), 'Colormap',jet);
    b.Title = ' S Information';
    b.XLabel = 'k';
    b.YLabel = 'g';

    Oinfo_norm = Oinfo ./ Sinfo;

    fig = figure();
    b = heatmap([k(1), k(5:5:end-4),k(end)],flip(g(1:5:end)), ... 
                flip([Oinfo_norm(1:5:end,1), Oinfo_norm(1:5:end,5:5:end-4), ...
                Oinfo_norm(1:5:end,end)]), 'Colormap',jet);
    b.Title = ' O Information normalised';
    b.XLabel = 'k';
    b.YLabel = 'g';

end