c20 = "Results/S8_c20/Sleep/Iccs/";
c2 = "Results/S8_c2/Sleep/Iccs/";
subject = "S8";

[Un_X_2, Un_Y_2, Red_2, Syn_2, MI_2] = MonkeyPIDLoader(c2 + "data.txt");
[Un_X_20, Un_Y_20, Red_20, Syn_20, MI_20] = MonkeyPIDLoader(c20 + "data.txt");

DeltaS_1 = Syn_2{1} - Syn_20{1};
DeltaMI_1 = MI_2{1} - MI_20{1};
DeltaR_1 = Red_2{1} - Red_20{1};
DeltaUx_1 = Un_X_2{1} - Un_X_20{1};
DeltaUy_1 = Un_X_2{1} - Un_X_20{1};

DeltaS_p = Syn_2{2} - Syn_20{2};
DeltaMI_p = MI_2{2} - MI_20{2};
DeltaR_p = Red_2{2} - Red_20{2};
DeltaUx_p = Un_X_2{2} - Un_X_20{2};
DeltaUy_p = Un_X_2{2} - Un_X_20{2};

ViolinPlot(DeltaUx_1, DeltaUy_1, DeltaR_1, DeltaS_1, DeltaMI_1, ...
                 subject+" - Var(1) comparison c2-c20", c2 + "../../channels_comp_violins_Var(1).png");
ViolinPlot(DeltaUx_p, DeltaUy_p, DeltaR_p, DeltaS_p, DeltaMI_p, ...
                 subject+" - Var(p) comparison c2-c20", c20 + "../../channels_comp_violins_Var(p).png");

% normalise w.r.t. MI

[Un_X_2_norm, Un_Y_2_norm, Red_2_norm, Syn_2_norm, MI_2_norm] = NormPID(Un_X_2, Un_Y_2, Red_2, Syn_2, MI_2);
[Un_X_20_norm, Un_Y_20_norm, Red_20_norm, Syn_20_norm, MI_20_norm] = NormPID(Un_X_20, Un_Y_20, Red_20, Syn_20, MI_20);

DeltaS_1_norm = Syn_2_norm{1} - Syn_20_norm{1};
DeltaMI_1_norm = MI_2_norm{1} - MI_20_norm{1};
DeltaR_1_norm = Red_2_norm{1} - Red_20_norm{1};
DeltaUx_1_norm = Un_X_2_norm{1} - Un_X_20_norm{1};
DeltaUy_1_norm = Un_X_2_norm{1} - Un_X_20_norm{1};

DeltaS_p_norm = Syn_2_norm{2} - Syn_20_norm{2};
DeltaMI_p_norm = MI_2_norm{2} - MI_20_norm{2};
DeltaR_p_norm = Red_2_norm{2} - Red_20_norm{2};
DeltaUx_p_norm = Un_X_2_norm{2} - Un_X_20_norm{2};
DeltaUy_p_norm = Un_X_2_norm{2} - Un_X_20_norm{2};

ViolinPlot(DeltaUx_1_norm, DeltaUy_1_norm, DeltaR_1_norm, DeltaS_1_norm, DeltaMI_1_norm, ...
                 subject+" - Var(1) comparison c2-c20 Normalised", c2 + "../../channels_comp_violins_Var(1)_norm.png");
ViolinPlot(DeltaUx_p_norm, DeltaUy_p_norm, DeltaR_p_norm, DeltaS_p_norm, DeltaMI_p_norm, ...
                 subject+" - Var(p) comparison c2-c20 Normalised", c20 + "../../channels_comp_violins_Var(p)_norm.png");