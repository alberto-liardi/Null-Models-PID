%% PID calculator for Gaussian system with MMI redundancy definition
% Sigma - Full covariance of the system
% S - Number of Sources
% T - Number of Targets

function [UnX, UnY, Red, Syn, MI_S1, MI_S2, MI] = PID_MMI_Gaussian(Sigma,S,T)

    [~,pd] = chol(Sigma);
    assert(pd==0, "The covariance matrix provided is not positive definite!");
    assert(size(Sigma,1)==S+T, "Mismatch between number of sources/targets and covariance dimensions.");

    Sigma_t = Sigma(end-T+1:end,end-T+1:end);
    Sigma_s = Sigma(1:S,1:S);
    Sigma_cross = Sigma(1:S,S+1:end)';

    MI = entropy(Sigma_t)+entropy(Sigma_s)-entropy(Sigma);
    
    % Sigma_s, Sigma_cross, Sigma_t, Sigma_full
    
    % calulate single MI (considering 2 sources)
            
    % S1
    Sigma_TS1 = cell(2,2);
    Sigma_TS1{1,1} = Sigma_t;
    Sigma_TS1{1,2} = Sigma_cross(:,1:S/2);
    Sigma_TS1{2,1} = Sigma_TS1{1,2}';
    Sigma_TS1{2,2} = Sigma_s(1:S/2,1:S/2);
    Sigma_TS1 = cell2mat(Sigma_TS1);
    MI_S1 = entropy(Sigma_t)+entropy(Sigma_s(1:S/2,1:S/2)) ...
            - entropy(Sigma_TS1);
    
    % S2
    Sigma_TS2 = cell(2,2);
    Sigma_TS2{1,1} = Sigma_t;
    Sigma_TS2{1,2} = Sigma_cross(:,S/2+1:end);
    Sigma_TS2{2,1} = Sigma_TS2{1,2}';
    Sigma_TS2{2,2} = Sigma_s(S/2+1:end,S/2+1:end);
    Sigma_TS2 = cell2mat(Sigma_TS2);
    MI_S2 = entropy(Sigma_t)+entropy(Sigma_s(S/2+1:end,S/2+1:end)) ...
            - entropy(Sigma_TS2);
    
    
    % now calculate PID atoms
    
    % MMI definition
    Red = min(MI_S1, MI_S2);
    UnX = MI_S1 - Red;
    UnY = MI_S2 - Red;
    Syn = MI - UnY - UnX - Red;

end


function H = entropy(Sig)
    
    arg = det(2*pi*exp(1)*Sig);    
    H = 0.5*log(arg);

end