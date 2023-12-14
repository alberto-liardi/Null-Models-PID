%% computing integrated information metrics with VAR model

function [phi_WMS, phi_R] = integrated_info_VAR(MI, Gamma)

    [MIxx, MIxy, MIyx, MIyy] = MI_parts(Gamma);

    phi_WMS = MI - MIxx - MIyy; 
    % NB: using MMI redundancy here
    phi_R = phi_WMS + min([MIxx, MIxy, MIyx, MIyy]);

    if ~isreal(phi_R)
        disp("well...");
    end

end

% calculation of the mutual information between the parts of a bipartite system

function [MIxx, MIxy, MIyx, MIyy] = MI_parts(g)

    assert(iscell(g), "g matrices must be in a 1x(p+1) cell");
    assert(length(g)==2, "Calculating MI between the parts is valid only for p=1!")

    S = size(g{1},1);

    % calculating the covariances with lags t and t-1
    Var_x = g{1}(1:S/2,1:S/2);
    Var_y = g{1}(S/2+1:end,S/2+1:end);

    Sigma_xx = [Var_x g{2}(1:S/2,1:S/2); g{2}(1:S/2,1:S/2)' Var_x];
    Sigma_xy = [Var_x g{2}(1:S/2,S/2+1:end); g{2}(1:S/2,S/2+1:end)' Var_y];
    Sigma_yx = [Var_y g{2}(S/2+1:end,1:S/2); g{2}(S/2+1:end,1:S/2)' Var_x];
    Sigma_yy = [Var_y g{2}(S/2+1:end,S/2+1:end); g{2}(S/2+1:end,S/2+1:end)' Var_y];

    MIxx = 2*entropy(Var_x) - entropy(Sigma_xx);
    MIxy = entropy(Var_x) + entropy(Var_y) - entropy(Sigma_xy);
    MIyx = entropy(Var_y) + entropy(Var_x) - entropy(Sigma_yx);
    MIyy = 2*entropy(Var_y) - entropy(Sigma_yy);

end

function H = entropy(Sigma)
    
    arg = det(2*pi*exp(1)*Sigma);    
    H = 0.5*log(arg)/log(2);
    % NB: log(2) since mutual information is in BITS

end