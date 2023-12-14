%% PROGRAM TO COMPARE THE PID ATOMS FROM WAKE AND REST STATES (OR PLACEBO AND DRUG)

% aw = "Results/Su_c2/Awake/MMI/";
% as = "Results/Su_c2/Sedated/MMI/";

% set(0,'defaulttextinterpreter','latex');
% set(0,'DefaultTextFontname', 'CMU Serif');

% aw = "Results/S15_c2/Awake/MMI/";
% as = "Results/S15_c2/Sleep/MMI/";
% subject = "S15";

chan = '10';
drug = "KET";

for ss = 1
aw = "Results_psychedelics/"+drug+"/"+ss+"_c"+chan+"/"+drug+"/MMI/";
as = "Results_psychedelics/"+drug+"/"+ss+"_c"+chan+"/Placebo/MMI/";
subject = ss;

[Un_X_W, Un_Y_W, Red_W, Syn_W, MI_W] = MonkeyPIDLoader(aw + "data.txt");
[Un_X_S, Un_Y_S, Red_S, Syn_S, MI_S] = MonkeyPIDLoader(as + "data.txt");

DeltaS_1 = Syn_W{1} - Syn_S{1};
DeltaMI_1 = MI_W{1} - MI_S{1};
DeltaR_1 = Red_W{1} - Red_S{1};
DeltaX_1 = Un_X_W{1} - Un_X_S{1};
DeltaY_1 = Un_Y_W{1} - Un_Y_S{1};
    
DeltaS_p = Syn_W{2} - Syn_S{2};
DeltaMI_p = MI_W{2} - MI_S{2};
DeltaR_p = Red_W{2} - Red_S{2};
DeltaX_p = Un_X_W{2} - Un_X_S{2};
DeltaY_p = Un_Y_W{2} - Un_Y_S{2};

ViolinPlot(DeltaX_1, DeltaY_1, DeltaR_1, DeltaS_1, DeltaMI_1, ...
                 subject+" - Var(1) comparison W-S", aw + "../../comp_violins_Var(1).png");
ViolinPlot(DeltaX_p, DeltaY_p, DeltaR_p, DeltaS_p, DeltaMI_p, ...
                 subject+" - Var(p) comparison W-S", aw + "../../comp_violins_Var(p).png");

% normalise w.r.t. MI

[Un_X_W_norm, Un_Y_W_norm, Red_W_norm, Syn_W_norm, MI_W_norm] = NormPID(Un_X_W, Un_Y_W, Red_W, Syn_W, MI_W);
[Un_X_S_norm, Un_Y_S_norm, Red_S_norm, Syn_S_norm, MI_S_norm] = NormPID(Un_X_S, Un_Y_S, Red_S, Syn_S, MI_S);

DeltaS_1_norm = Syn_W_norm{1} - Syn_S_norm{1};
DeltaMI_1_norm = MI_W_norm{1} - MI_S_norm{1};
DeltaR_1_norm = Red_W_norm{1} - Red_S_norm{1};
DeltaX_1_norm = Un_X_W_norm{1} - Un_X_S_norm{1};
DeltaY_1_norm = Un_Y_W_norm{1} - Un_Y_S_norm{1};

DeltaS_p_norm = Syn_W_norm{2} - Syn_S_norm{2};
DeltaMI_p_norm = MI_W_norm{2} - MI_S_norm{2};
DeltaR_p_norm = Red_W_norm{2} - Red_S_norm{2};
DeltaX_p_norm = Un_X_W_norm{2} - Un_X_S_norm{2};
DeltaY_p_norm = Un_Y_W_norm{2} - Un_Y_S_norm{2};

ViolinPlot(DeltaX_1_norm, DeltaY_1_norm, DeltaR_1_norm, DeltaS_1_norm, DeltaMI_1_norm, ...
                 subject+" - Var(1) comparison W-S Normalised", aw + "../../comp_violins_Var(1)_norm.png");
ViolinPlot(DeltaX_p_norm, DeltaY_1_norm, DeltaR_p_norm, DeltaS_p_norm, DeltaMI_p_norm, ...
                 subject+" - Var(p) comparison W-S Normalised", aw + "../../comp_violins_Var(p)_norm.png");
end

%% Use null distribution and quantiles for comparison

% rng('default') % for reproducibility of random variables
% subject = '1';
% aw = "Results_psychedelics/PSIL/"+subject+"_c2/PSIL/MMI/";
% as = "Results_psychedelics/PSIL/"+subject+"_c2/Placebo/MMI/";
ks = 2;


InfoDyn = readmatrix(aw+"../InfoDyn.txt");
ps = InfoDyn(:,7:8);
model = struct('name',"VAR",'S',str2double(chan),'n',100);

for p = ks
    fprintf("doing p= "+sprintf('%01d',p)+"\n");
    model.p = ps(:,p);

    [Un_X_W_quant{p}, Un_Y_W_quant{p}, Red_W_quant{p}, Syn_W_quant{p}] = ...
        QuantPID(Un_X_W{p}, Un_Y_W{p}, Red_W{p}, Syn_W{p}, MI_W{p},model);
    [Un_X_S_quant{p}, Un_Y_S_quant{p}, Red_S_quant{p}, Syn_S_quant{p}] = ...
        QuantPID(Un_X_S{p}, Un_Y_S{p}, Red_S{p}, Syn_S{p}, MI_S{p},model);

    % saving the quantiles for future studies
    SaveQuantPID(Un_X_W_quant{p}, Un_Y_W_quant{p}, Red_W_quant{p}, ...
        yn_W_quant{p}, MI_W{p}, aw, p);
    SaveQuantPID(Un_X_S_quant{p}, Un_Y_S_quant{p}, Red_S_quant{p}, ...
        Syn_S_quant{p}, MI_S{p}, as, p);

end

DeltaS_1_quant = Syn_W_quant{1} - Syn_S_quant{1};
DeltaMI_1_quant = zeros(size(DeltaS_1_quant));
DeltaR_1_quant = Red_W_quant{1} - Red_S_quant{1};
DeltaX_1_quant = Un_X_W_quant{1} - Un_X_S_quant{1};
DeltaY_1_quant = Un_Y_W_quant{1} - Un_Y_S_quant{1};

DeltaS_p_quant = Syn_W_quant{2} - Syn_S_quant{2};
DeltaMI_p_quant = zeros(size(DeltaS_p_quant));
DeltaR_p_quant = Red_W_quant{2} - Red_S_quant{2};
DeltaX_p_quant = Un_X_W_quant{2} - Un_X_S_quant{2};
DeltaY_p_quant = Un_Y_W_quant{2} - Un_Y_S_quant{2};

ViolinPlot(DeltaX_1_quant, DeltaY_1_quant, DeltaR_1_quant, DeltaS_1_quant, DeltaMI_1_quant, ...
                 subject+" - Var(1) comparison W-S with quantiles", aw + "../../comp_violins_Var(1)_quant_temp.png");
ViolinPlot(DeltaX_p_quant, DeltaY_p_quant, DeltaR_p_quant, DeltaS_p_quant, DeltaMI_p_quant, ...
                 subject+" - Var(p) comparison W-S with quantiles", aw + "../../comp_violins_Var(p)_quant_temp.png");
