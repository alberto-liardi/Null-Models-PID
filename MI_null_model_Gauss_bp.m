%% MI-PID ATOMS NULL MODEL FOR GAUSSIAN SYSTEMS

function PIDs = MI_null_model_Gauss_bp(MI,Sources,Targets,N_runs,red_fun)

    if (~exist('Sources','var') || isempty(Sources)), Sources = 2; end
    if ~exist('Targets','var') || isempty(Targets), Targets = 1; end
    if ~exist('N_runs','var') || isempty(N_runs), N_runs = 1e4; end
    if ~exist('red_fun','var') || isempty(red_fun), red_fun = "MMI"; end

    Syn = NaN(1,N_runs); Red = NaN(1,N_runs);
    UnX = NaN(1,N_runs); UnY = NaN(1,N_runs);
    
    x0 = [1 2 0 5 10 20 50];
    opts = optimset('display', 'none');
    cc = 0;
    S = Sources;
    T = Targets;
        
    for j = 1:N_runs
        
        if(mod(j,N_runs/10)==0), disp(j); end
        
        
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
        if( code == -3 || alpha<0 ), cc = cc+1; continue; end 
        
        g = alpha^(-1);
        
        % check
        % 0.5*log(det(2*pi*exp(1)*(A*Sigma_s*A'+g*Sigma_u))) ...
        % - 0.5*log(det(2*pi*exp(1)*(g*Sigma_u)))
        
        % get the covariance matrix
        Sigma_t = A * Sigma_s * A.' + g*Sigma_u;
        
        % calculate cross covariance terms
        Sigma_cross = A*Sigma_s;
    
        Sigma_full = [Sigma_s, Sigma_cross';
                      Sigma_cross, Sigma_t];
    
        I_check = entropy(Sigma_t)+entropy(Sigma_s)-entropy(Sigma_full);
        
        if(abs(MI-I_check)>1e-7)
            disp("problem with MI"); disp(abs(MI-I_check)); 
            continue
        end
        
        
        % now calculate PID atoms using {red_fun} definition
        if red_fun=="MMI"
            % MMI
            [UnX(1,j), UnY(1,j), Red(1,j), Syn(1,j)] = PID_MMI_Gaussian(Sigmal_full,S,T);
        elseif red_fun=="Idep"
            % Idep
            PI = num2cell(calc_pi_Idep_mvn(Sigma_full, [S/2 S/2 T]));
            [Red(1,j), UnX(1,j), UnY(1,j), Syn(1,j)] = PI{:};
        elseif red_fun=="Iccs"
            % Iccs 
            PI = num2cell(calc_pi_mvn(lattice2d(), Sigma_full, [S/2 S/2 T], @Iccs_mvn_P2).PI);
            [Red(1,j), UnX(1,j), UnY(1,j), Syn(1,j)] = PI{:};
        end
    
        try
            PIDcheck([UnX(1,j), UnY(1,j), Red(1,j), Syn(1,j)]);
        catch
            cc = cc+1;
            continue
        end
        assert(isreal(Syn(1,j)), "not real");
    
    end

    PIDs = [UnX; UnY; Red; Syn];
end

function H = entropy(Sigma)
    
    arg = det(2*pi*exp(1)*Sigma);    
    H = 0.5*log(arg);

end

function y = fun(x,Sigma_u,A,Sigma_s,I)

    T = length(Sigma_u);
    y = det(eye(T)+ x * (Sigma_u \ A * Sigma_s' * A.'))-exp(2*I);

end