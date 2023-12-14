%% MI-PID ATOMS NULL MODEL FOR AUTOREGRESSIVE SYSTEMS

function [PIDs, Phis] = MI_null_model_VAR(MI,S,p,N_runs,red_fun,metric)
    
    if ~exist('S','var') || isempty(S), S = 2; end
    if ~exist('p','var') || isempty(p), p = 1; end
    if ~exist('N_runs','var') || isempty(N_runs), N_runs = 1e4; end
    if ~exist('red_fun','var') || isempty(red_fun), red_fun = ["MMI"]; end
    if ~exist('metric','var') || isempty(metric), metric = ['y','n']; end

    if mod(S,2)~=0, error("Provide even number of sources!"); end

    PID_red = ['n','n','n']; rind=[];
    if ismember("MMI",red_fun), PID_red(1) = 'y'; rind(end+1)=1; end
    if ismember("Idep",red_fun), PID_red(2) = 'y'; rind(end+1)=2; end
    if ismember("Iccs",red_fun), PID_red(3) = 'y'; rind(end+1)=3; end
    if ismember("all",red_fun), PID_red(:) = 'y'; rind=1:3; end

    Syn = NaN(3,N_runs); Red = NaN(3,N_runs);
    UnX = NaN(3,N_runs); UnY = NaN(3,N_runs);
    PIDs = NaN(4,N_runs,3);

    phi_WMS = NaN(1,N_runs); phi_R = NaN(1,N_runs);
    red_phi = NaN(1,N_runs);

    cc=0;

%     if MI>43, fprintf("MI is "+sprintf('%02d',MI)+"\n"); end
    
    for j = 1:N_runs

        if mod(j,N_runs/10)==0, disp(j); end

        [A,V,check] = VAR_from_MI(MI,p,S);
        
        if(check==1)
            disp("problem with MI optimisation"); 
            cc = cc+1;
            continue
        end

        if metric(1)=='y'
            PID = PIDcalc_Mult(p,A,V,S/2,S/2,"joint",PID_red);
            
            try
                PIDcheck(PID(1,:));
                assert(isreal(PID), "Not real");
            catch 
                cc = cc+1;
                disp("PID values not real or not consistent!");
                continue;
            end

            for r = rind
                MI_calc = sum(PID(r,:));
    
                if(abs(MI_calc-MI)>MI/1e4), cc=cc+1; 
                    disp("Mismatch between true MI and MI calculated during the optimisation!")
                    disp(abs(MI_calc-MI)); %disp(MI); 
                    continue; 
                end
            end
                    
            UnX(:,j) = PID(:,1);    UnY(:,j) = PID(:,2);
            Red(:,j) = PID(:,3);    Syn(:,j) = PID(:,4);

        end
        
        if metric(2)=='y'
            Gs = var_to_autocov(B,V,p);
            MI_calc = (0.5*log(det(Gs(:,:,1)))-0.5*log(det(V)))/log(2);

            Gs_cell = cell(1,p+1);
            for s = 1:p+1
                Gs_cell{s} = Gs(:,:,s);
            end

            [phi_WMS(j), phi_R(j)] = integrated_info_VAR(MI_calc,Gs_cell);
            red_phi(j) =  phi_R(j)-phi_WMS(j);
        end
  
    end

    for r = 1:3
        PIDs(:,:,r) = [UnX(r,:); UnY(r,:); Red(r,:); Syn(r,:)];
    end
    Phis = [phi_WMS; phi_R; red_phi];
%     fprintf("done! excluded values were " + sprintf('%01d',cc)+"\n");
end
