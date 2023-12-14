%% MI-PID ATOMS NULL MODEL FOR GAUSSIAN SYSTEMS

function PIDs = MI_null_model_Gauss(MI,Sources,Targets,N_runs,red_fun)

    if (~exist('Sources','var') || isempty(Sources)), Sources = 2; end
    if ~exist('Targets','var') || isempty(Targets), Targets = 1; end
    if ~exist('N_runs','var') || isempty(N_runs), N_runs = 1e4; end
    if ~exist('red_fun','var') || isempty(red_fun), red_fun = "MMI"; end

    Syn = NaN(3,N_runs); Red = NaN(3,N_runs);
    UnX = NaN(3,N_runs); UnY = NaN(3,N_runs);
    PIDs = NaN(4,N_runs,3);
    
    cc = 0;
    S = Sources;
    T = Targets;
        
    for j = 1:N_runs
        
        if(mod(j,N_runs/10)==0), disp(j); end
        
        % generate random Gaussian system with specific MI
        [Sigma_full, check] = Gaussian_from_MI(MI,S,T);
    
        if(check==1)
            disp("problem with MI optimisation"); 
            cc = cc+1;
            continue
        end

        
        % now calculate PID atoms using {red_fun} definition
        if red_fun=="MMI" || red_fun=="all"
            % MMI
            [UnX(1,j), UnY(1,j), Red(1,j), Syn(1,j)] = PID_MMI_Gaussian(Sigma_full,S,T);
        end
        if red_fun=="Idep" || red_fun=="all"
            % Idep
            [UnX(2,j), UnY(2,j), Red(2,j), Syn(2,j)] = PID_Idep_Gaussian(Sigma_full,S,T);
        end
        if red_fun=="Iccs" || red_fun=="all"
            % Iccs 
            [UnX(3,j), UnY(3,j), Red(3,j), Syn(3,j)] = PID_Iccs_Gaussian(Sigma_full,S,T);
        end
    
        % check consistency of PID atoms (only MMI)
        try
            PIDcheck([UnX(1,j), UnY(1,j), Red(1,j), Syn(1,j)]);
        catch
            cc = cc+1;
            continue
        end
        assert(isreal(Syn(1,j)), "not real");
    
    end

    for n = 1:3
        PIDs(:,:,n) = [UnX(n,:); UnY(n,:); Red(n,:); Syn(n,:)];
    end
end
