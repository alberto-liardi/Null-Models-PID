%% MI-PID ATOMS NULL MODEL FOR DISCRETE SYSTEMS

function PIDs = MI_null_model_Discrete(MI,N_runs)

%     if ~exist('boolean','var') || isempty(gate), gate = "XOR"; end
    if ~exist('N_runs','var') || isempty(N_runs), N_runs = 7e4; end

    % N_runs is the number of runs of the OR gate.
    % NB: XOR gate must be run half times that value for consistency

    % set parameters for the study 
    x0 = [1 -1 2 -2 5 -5];
    opts = optimset('display', 'none');
    cc = 0;
    
    assert(MI<-log(0.5), ...
           "Value of Mutual Information is incompatible with the system!")
    
    Syn = NaN(1,N_runs); Red = NaN(1,N_runs);
    UnX = NaN(1,N_runs); UnY = NaN(1,N_runs);
    
    % sample the parameters for the dirichlet from a uniform distribution
    a = repmat(unifrnd(0.1, 3),1,4);

    runarr = 1:floor(N_runs/7);
    
    for gate = ["XOR","XOR2","XOR3","OR","OR2","OR3","OR4"]

        fprintf("doing gate %s!\n", gate);

        if gate == "XOR"
            boolfun = @XORgate;
            js = runarr;
        elseif gate == "XOR2"
            boolfun = @XOR2gate;
            js = floor(N_runs/7)+[runarr];
        elseif gate == "XOR3"
            boolfun = @XOR3gate;
            js = floor(2*N_runs/7)+[runarr];
        elseif gate == "OR"
            boolfun = @ORgate;
            js = floor(3*N_runs/7)+[runarr];
        elseif gate == "OR2"
            boolfun = @OR2gate;
            js = floor(4*N_runs/7)+[runarr];
        elseif gate == "OR3"
            boolfun = @OR3gate;
            js = floor(5*N_runs/7)+[runarr];
        elseif gate == "OR4"
            boolfun = @OR4gate;
            js = floor(6*N_runs/7)+[runarr];
        end
    
        for j = js

            if mod(j,N_runs/10)==0, disp(j), end

            % sample the probabilities of the Sources S from the dirichlet (only once)
            ps = drchrnd(1,a);
            
            % get the probability of Z=f(X,Y)
            [pz0,pz1] = boolfun(ps);
        
            % for some parameters of the Dirichlet the MI value may not be reached,
            % check with the maximum of MI (at x=1 or x=0) and skip the step if needed
            if -pz1*log(pz1)-pz0*log(pz0)<MI, disp("skipping step!"), continue; end
            
            % find the discrete noise epsilon to get the value of MI
            l = 1;
            [f, ~, code, ~] = fzero(@(x) discrfun(x, pz0, pz1, MI), x0(l), opts);
            
            while( code ~= 1 && l+1 <= length(x0) )
                l = l+1;
                [f, ~, code, ~] = fzero(@(x) discrfun(x, pz0, pz1, MI), x0(l), opts);
            end
            if( code ~= 1 ), cc = cc+1; fprintf("optimisation failed "+sprintf('%01d',code)+"\n"); 
                             continue; end 
            
            % now do one more calculation
            pe = Sigmoid(f,0.5);
            if ~isfinite(discrfun(f, pz0, pz1, MI))
                cc = cc+1;
                disp("Optimisation function returned NaN value!");
                continue;
            end
            pt0 = pz0*(1-pe) + pz1*pe;
            pt1 = pz1*(1-pe) + pz0*pe; % 1-pt0;
            
            HT = -pt0*log(pt0)-pt1*log(pt1);
            HTcondZ = -pe*log(pe)-(1-pe)*log(1-pe);
            MI_calc = HT - HTcondZ;
            if(abs(MI_calc-MI)>MI/1e4), cc=cc+1; 
                disp("Mismatch between true MI and MI calculated during the optimisation!")
                disp(abs(MI_calc-MI)); %disp(MI); 
                continue; end
            
            % now need to calculate the marginals MIx, MIy
            [MIx, MIy] = DiscreteMarginalMI(ps,pe,HT,gate);
        
            Red(j) = min(MIx, MIy);
            UnX(j) = MIx - Red(j);
            UnY(j) = MIy - Red(j);
            Syn(j) = MI_calc - UnX(j) - UnY(j) - Red(j);
    
    %         if Red(j)<0
    %             disp("problems...");
    %         end
    
    %         PIDcheck([UnX(j), UnY(j), Red(j), Syn(j)]);
    
        end
    end

    PIDs = [UnX; UnY; Red; Syn];

end

function y = discrfun(x, pz0, pz1, MI_value)

%     disp(x); disp("here");
    
    x = Sigmoid(x,0.5);
%     disp(x);
    % can calculate the entropy of the target T=Z+epsilon manually:
    pt0 = pz0*(1-x) + pz1*x;
    pt1 = pz1*(1-x) + pz0*x; % 1-pt0;
    
    HT = -pt0*log(pt0)-pt1*log(pt1);
    HTcondZ = -x*log(x)-(1-x)*log(1-x);
    I = HT - HTcondZ;
    
    y = MI_value - I;

end

function y = Sigmoid(x,M)
    
    y = M/(1+exp(-x));

end

function [pz0,pz1] = XORgate(ps)
    pz0 = ps(1)+ps(4);
    pz1 = 1-pz0;
end

function [pz0,pz1] = XOR2gate(ps)
    pz0 = ps(1)+ps(2);
    pz1 = 1-pz0;
end

function [pz0,pz1] = XOR3gate(ps)
    pz0 = ps(1)+ps(3);
    pz1 = 1-pz0;
end

function [pz0,pz1] = ORgate(ps)
    pz0 = ps(1);
    pz1 = 1-ps(1);
end

function [pz0,pz1] = OR2gate(ps)
    pz0 = ps(2);
    pz1 = 1-pz0;
end

function [pz0,pz1] = OR3gate(ps)
    pz0 = ps(3);
    pz1 = 1-pz0;
end

function [pz0,pz1] = OR4gate(ps)
    pz0 = ps(4);
    pz1 = 1-pz0;
end


function r = drchrnd(n,a)
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
end