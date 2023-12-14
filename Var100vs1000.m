as100 = "Results_old/small_Kin2_c20/Sleep/MMI/";
aw100 = "Results_old/small_Kin2_c20/Awake/MMI/";
as1000 = "Results_old/Kin2_c20/Sedated/MMI/";
aw1000 = "Results_old/Kin2_c20/Awake/MMI/";
subject = "Kin2";

[Un_X_W_1000, Un_Y_W_1000, Red_W_1000, Syn_W_1000, MI_W_1000] = MonkeyPIDLoader(aw1000 + "data.txt");
[Un_X_S_1000, Un_Y_S_1000, Red_S_1000, Syn_S_1000, MI_S_1000] = MonkeyPIDLoader(as1000 + "data.txt");
[Un_X_W_1000, Un_Y_W_1000, Red_W_1000, Syn_W_1000, MI_W_1000] = NormPID(Un_X_W_1000, Un_Y_W_1000, Red_W_1000, Syn_W_1000, MI_W_1000);
[Un_X_S_1000, Un_Y_S_1000, Red_S_1000, Syn_S_1000, MI_S_1000] = NormPID(Un_X_S_1000, Un_Y_S_1000, Red_S_1000, Syn_S_1000, MI_S_1000);

[Un_X_W_100, Un_Y_W_100, Red_W_100, Syn_W_100, MI_W_100] = MonkeyPIDLoader(aw100 + "data.txt");
[Un_X_S_100, Un_Y_S_100, Red_S_100, Syn_S_100, MI_S_100] = MonkeyPIDLoader(as100 + "data.txt");
[Un_X_W_100, Un_Y_W_100, Red_W_100, Syn_W_100, MI_W_100] = NormPID(Un_X_W_100, Un_Y_W_100, Red_W_100, Syn_W_100, MI_W_100);
[Un_X_S_100, Un_Y_S_100, Red_S_100, Syn_S_100, MI_S_100] = NormPID(Un_X_S_100, Un_Y_S_100, Red_S_100, Syn_S_100, MI_S_100);

DeltaS_1000 = Syn_W_1000{1} - Syn_S_1000{1};
DeltaMI_1000 = MI_W_1000{1} - MI_S_1000{1};
DeltaR_1000 = Red_W_1000{1} - Red_S_1000{1};
DeltaX_1000 = Un_X_W_1000{1} - Un_X_S_1000{1};
DeltaY_1000 = Un_Y_W_1000{1} - Un_Y_S_1000{1};

DeltaS_100 = Syn_W_100{1} - Syn_S_100{1};
DeltaMI_100 = MI_W_100{1} - MI_S_100{1};
DeltaR_100 = Red_W_100{1} - Red_S_100{1};
DeltaX_100 = Un_X_W_100{1} - Un_X_S_100{1};
DeltaY_100 = Un_Y_W_100{1} - Un_Y_S_100{1};

mean(DeltaX_1000), mean(DeltaX_100)
mean(DeltaY_1000), mean(DeltaY_100)
mean(DeltaS_1000), mean(DeltaS_100)
mean(DeltaR_1000), mean(DeltaR_100)
mean(DeltaMI_1000), mean(DeltaMI_100)
