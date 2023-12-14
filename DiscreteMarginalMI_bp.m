function [MIx, MIy] = DiscreteMarginalMI(ps,pe,HT,gate)

    if gate=="XOR"
        [MIx, MIy] = XORgate(ps,pe,HT);
    elseif gate=="OR"
        [MIx, MIy] = ORgate(ps,pe,HT);
    end
end
    
    
function [MIx, MIy] = XORgate(ps,pe,HT)

    % unpack probabilities
    cps = num2cell(ps);
    [p00,p01,p10,p11] = cps{:};

    % single probabilities
    px0 = p00 + p01;
    px1 = 1-px0;
    
    py0 = p00 + p10;
    py1 = 1-py0;

    % conditional probabilities on Z
    pZ0_x0 = p00/px0;
    pZ1_x0 = p01/px0;
    pZ0_x1 = p11/px1;
    pZ1_x1 = p10/px1;
    
    % conditional probabilities on T
    pT1_x1 = py0*(1-pe) + py1*pe; % = pT0_x0
    pT1_x0 = py1*(1-pe) + py0*pe; % 1-pT1_x1; % = pT0_x1
    pT0_x1 = pT1_x0;
    pT0_x0 = pT1_x1;
    
    % conditional entropies
    HT_x0 = -px0*(pT0_x0*log(pT0_x0) + pT1_x0*log(pT1_x0));
    HT_x1 = -px1*(pT0_x1*log(pT0_x1) + pT1_x1*log(pT1_x1));
    HTcondX = HT_x0 + HT_x1;
    
    % ...do the same for y
    pT1_y1 = px0*(1-pe) + px1*pe; % = pT0_y0
    pT1_y0 = px1*(1-pe) + px0*pe; % = pT0_y1
    pT0_y1 = pT1_y0;
    pT0_y0 = pT1_y1;
    
    HT_y0 = -py0*(pT0_y0*log(pT0_y0) + pT1_y0*log(pT1_y0));
    HT_y1 = -py1*(pT0_y1*log(pT0_y1) + pT1_y1*log(pT1_y1));
    HTcondY = HT_y0 + HT_y1;

%     HTcondX = - (pT1_x1*log(pT1_x1) + pT1_x0*log(pT1_x0));
%     HTcondY = - (pT1_y1*log(pT1_y1) + pT1_y0*log(pT1_y0));
    
    % mutual information
    MIx = HT - HTcondX;
    MIy = HT - HTcondY;

    if MIx<0 || MIy<0
        disp("problems...");
    end
    
end

function [MIx, MIy] = ORgate(ps,pe,HT)
    % single probabilities
    px0 = ps(1)+ps(2);
    px1 = 1-px0;
    
    py0 = ps(1)+ps(3);
    py1 = 1-py0;
    
    % conditional probabilities
    pT0_x0 = py0*(1-pe) + py1*pe; 
    pT1_x0 = py1*(1-pe) + py0*pe; % 1-pT1_x1; 
    pT1_x1 = 1-pe;
    pT0_x1 = pe;

    % conditional entropies
    HT_x0 = -px0*(pT0_x0*log(pT0_x0) + pT1_x0*log(pT1_x0));
    HT_x1 = -px1*(pT0_x1*log(pT0_x1) + pT1_x1*log(pT1_x1));
    HTcondX = HT_x0 + HT_x1;
    
    % ...do the same for y
    pT0_y0 = px0*(1-pe) + px1*pe; 
    pT1_y0 = px1*(1-pe) + px0*pe; 
    pT1_y1 = 1-pe;
    pT0_y1 = pe;
    
    HT_y0 = -py0*(pT0_y0*log(pT0_y0) + pT1_y0*log(pT1_y0));
    HT_y1 = -py1*(pT0_y1*log(pT0_y1) + pT1_y1*log(pT1_y1));
    HTcondY = HT_y0 + HT_y1;
    
    % mutual information
    MIx = HT - HTcondX;
    MIy = HT - HTcondY;
end