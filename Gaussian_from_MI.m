%% Get random Gaussian system with specific Mutual Information
% MI - total mutual information between sources and targets
% Sources - Number of sources
% Targets - Number of targets


function [Cov,flag,Coef,g] = Gaussian_from_MI(MI,S,T)

    if (~exist('S','var') || isempty(S)), S = 2; end
    if ~exist('T','var') || isempty(T), T = 1; end

    x0 = [1 2 0 5 10 20 50];
    opts = optimset('display', 'none');        
    flag=0;
    
    % sample the matrices
    % NB. for Wishart degrees of freedom must be > size - 1
    
    % source covariance
    Sigma_s = wishrnd(eye(S),S+1);
    % linear coefficients
    A = normrnd(0,1,T,S);
    % target conditional covariance
    Sigma_u = wishrnd(eye(T),T+1);
    
    % find g such that mutual information is MI_tot
    l = 1;
    [alpha, ~, code] = fzero(@(x) fun(x,Sigma_u,A,Sigma_s,MI), x0(l), opts);
    
    while( (code == -3 || alpha<0) && l+1 <= length(x0) )
        l = l+1;
        [alpha, ~, code] = fzero(@(x) fun(x,Sigma_u,A,Sigma_s,MI), x0(l), opts);
    end
    if( code == -3 || alpha<0 ), error("problem with MI"); end 
    
    g = alpha^(-1);
    
    % check
    % 0.5*log(det(2*pi*exp(1)*(A*Sigma_s*A'+g*Sigma_u))) ...
    % - 0.5*log(det(2*pi*exp(1)*(g*Sigma_u)))
    
    % get the covariance matrix
    Sigma_t = A * Sigma_s * A.' + g*Sigma_u;
    
    % calculate cross covariance terms
    Sigma_cross = A*Sigma_s;

    Sigma = [Sigma_s, Sigma_cross';
                  Sigma_cross, Sigma_t];

    % check the result of the optimisation
    I_check = entropy(Sigma_t)+entropy(Sigma_s)-entropy(Sigma);

    if(abs(MI-I_check)>1e-7)
        disp(abs(MI-I_check)); flag=1; error("problem with MI");
    end

    Cov = Sigma; Coef = A;

end


function H = entropy(Sigma)
    
    arg = det(2*pi*exp(1)*Sigma);    
    H = 0.5*log(arg);

end

function y = fun(x,Sigma_u,A,Sigma_s,I)

    T = length(Sigma_u);
    y = det(eye(T)+ x * (Sigma_u \ A * Sigma_s' * A.'))-exp(2*I);

end