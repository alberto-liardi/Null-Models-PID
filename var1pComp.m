%% VAR(1)-VAR(p) comparison of the data stored in PATH

% path = "Results/George_c20/Awake/MMI/";
path = "Results/S4_c2/Sleep/Iccs/";
subject = "S4";

% read from file all the necessary info
[Un_X, Un_Y, Red, Syn, MI] = MonkeyPIDLoader(path + "data.txt");

% comparison between Var(1) and Var(p)
DeltaR = Red{2}-Red{1};
DeltaS = Syn{2} - Syn{1};
DeltaUx = Un_X{2} - Un_X{1};
DeltaUy = Un_Y{2} - Un_Y{1};
DeltaMI = MI{2} - MI{1};

ViolinPlot(DeltaUx, DeltaUy, DeltaR, DeltaS, DeltaMI, subject+" - Var(p)-Var(1) PID atoms", ...
           path + "Violins.png");


% normalise w.r.t. MI

[Un_X_norm, Un_Y_norm, Red_norm, Syn_norm, MI_norm] = NormPID(Un_X, Un_Y, Red, Syn, MI);

DeltaX_norm = Un_X_norm{2} - Un_X_norm{1};
DeltaY_norm = Un_Y_norm{2} - Un_Y_norm{1};
DeltaS_norm = Syn_norm{2} - Syn_norm{1};
DeltaR_norm = Red_norm{2} - Red_norm{1};
DeltaMI_norm = MI_norm{2} - MI_norm{1};

ViolinPlot(DeltaX_norm, DeltaY_norm, DeltaR_norm, DeltaS_norm, DeltaMI_norm, ...
           subject+" - Var(p)-Var(1) PID atoms Normalised", path + "Violins_norm.png");

