S = 2;
p=1;
N_runs = 200;

gs_cell = cell(1,p);
phi_WMS = NaN(1,N_runs); phi_R = NaN(1,N_runs); 
red_phi = NaN(1,N_runs); MI = NaN(1,N_runs);

disp("begin run");

for j = 1:N_runs

%     if mod(j,N_runs/10)==0, disp(j); end
    % sample the matrices
    % NB. for Wishart degrees of freedom must be > size - 1
    
    % random VAR matrix
    A = var_rand(S,p,rand(1));
    while any(abs(A(:))<1e-18)
        A = var_rand(S,p,rand(1));
    end
    
    % target conditional covariance
    V = wishrnd(eye(S),S+1);

    gs = var_to_autocov(A,V,p);
    MI(j) = (0.5*log(det(gs(:,:,1)))-0.5*log(det(V)))/log(2);

    for s = 1:p+1
        gs_cell{s} = gs(:,:,s);
    end

    [phi_WMS(j), phi_R(j)] = integrated_info_VAR(MI(j),gs_cell);
    red_phi(j) =  phi_R(j)-phi_WMS(j);

               
end

Phis = [phi_WMS; phi_R; red_phi];

disp("done!");

%% quantiles

disp("begin run2");

model = struct('name',"VAR",'S',S,'n',100);

h = size(phi_WMS,2)/2;

[phi_wms_1_quant, phi_r_1_quant, red_1_quant] = ...
    QuantPhi(phi_WMS(1:h), phi_R(1:h), red_phi(1:h), MI(1:h), model);
[phi_wms_2_quant, phi_r_2_quant, red_2_quant] = ...
    QuantPhi(phi_WMS(h+1:end), phi_R(h+1:end), red_phi(h+1:end), MI(h+1:end), model);

Delta_phi_wms_quant = phi_wms_1_quant - phi_wms_2_quant;
Delta_phi_r_quant = phi_r_1_quant - phi_r_2_quant;
Delta_red_phi_quant = red_1_quant - red_2_quant;

ViolinPlot_small(Delta_phi_wms_quant, Delta_phi_r_quant, Delta_red_phi_quant, ...
                 [], "Var(1) Phi comparison with quantiles - random models", ...
                 ["\phi_{WMS}", "\phi_R", "Red (MMI)"]);

disp("mean Delta phi_WMS"); disp(mean(Delta_phi_wms_quant,2));
disp("mean Delta phi_R"); disp(mean(Delta_phi_r_quant,2));
disp("mean Delta phi_red"); disp(mean(Delta_red_phi_quant,2));

disp("done!");
