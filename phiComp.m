
channels = 10;
aw = "Results_p1_phi_test/George_c"+sprintf('%01d',channels)+"/Awake/";
as = "Results_p1_phi_test/George_c"+sprintf('%01d',channels)+"/Sedated/";

% aw = "Results_p1_phi/S2_c"+sprintf('%01d',channels)+"/Awake/";
% as = "Results_p1_phi/S2_c"+sprintf('%01d',channels)+"/Sleep/";
subject = "George";

data_w = readmatrix(aw+"InfoDyn.txt");
data_s = readmatrix(as+"InfoDyn.txt");

phi_wms_w = data_w(:,13)';
phi_r_w = data_w(:,14)';
phi_wms_s = data_s(:,13)';
phi_r_s = data_s(:,14)';

red_w = phi_r_w - phi_wms_w;
red_s = phi_r_s - phi_wms_s;

delta_phi_wms = phi_wms_w - phi_wms_s;
delta_phi_r = phi_r_w - phi_r_s;
delta_r = red_w - red_s;

ViolinPlot_small(delta_phi_wms, delta_phi_r, delta_r, ...
                 aw + "../phi_comp_violins.png", ...
                 subject+" - Calculation of \Phi for "+sprintf('%01d',channels)+" channels", ...
                 ["\phi_{WMS}", "\phi_R", "Red (MMI)"]);


% Normalise Phi measures w.r.t. MI

% loading MI
[~,~,~,~,MI_w] = MonkeyPIDLoader(aw+"/MMI/data.txt");
[~,~,~,~,MI_s] = MonkeyPIDLoader(as+"/MMI/data.txt");

norm_phi_wms_w = phi_wms_w./MI_w{1};
norm_phi_r_w = phi_r_w./MI_w{1};
norm_phi_wms_s = phi_wms_s./MI_s{1};
norm_phi_r_s = phi_r_s./MI_s{1};

norm_red_w = red_w./MI_w{1};
norm_red_s = red_s./MI_s{1};

norm_delta_phi_wms = norm_phi_wms_w - norm_phi_wms_s;
norm_delta_phi_r = norm_phi_r_w - norm_phi_r_s;
norm_delta_r = norm_red_w - norm_red_s;

ViolinPlot_small(norm_delta_phi_wms, norm_delta_phi_r, norm_delta_r, ...
                 aw + "../Normalised_phi_comp_violins.png", ...
                 subject+" - Normalised calculation of \Phi for "+sprintf('%01d',channels)+" channels", ...
                 ["\phi_{WMS}", "\phi_R", "Red (MMI)"]);



%% Use null distribution and quantiles for comparison

model = struct('name',"VAR",'S',channels,'n',100);

[phi_wms_w_quant, phi_r_w_quant, red_w_quant] = ...
    QuantPhi(phi_wms_w, phi_r_w, red_w, MI_w{1}, model);
[phi_wms_s_quant, phi_r_s_quant, red_s_quant] = ...
    QuantPhi(phi_wms_s, phi_r_s, red_s, MI_s{1}, model);

Delta_phi_wms_quant = phi_wms_w_quant - phi_wms_s_quant;
Delta_phi_r_quant = phi_r_w_quant - phi_r_s_quant;
Delta_red_phi_quant = red_w_quant - red_s_quant;

ViolinPlot_small(Delta_phi_wms_quant, Delta_phi_r_quant, Delta_red_phi_quant, ...
                 aw + "../Phi_comp_violins_Var(1)_quant.png", ...
                 subject+" - Var(1) Phi comparison W-S with quantiles", ...
                 ["\phi_{WMS}", "\phi_R", "Red (MMI)"]);
