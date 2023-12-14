%% MI-PID ATOMS NULL MODEL FOR AUTOREGRESSIVE SYSTEMS

function [PIDs, Phis] = MI_null_model_VAR(MI,Sources,p,N_runs,red_fun,metric)
    
    if ~exist('Sources','var') || isempty(Sources), S = 2; else S = Sources; end
    if ~exist('p','var') || isempty(p), p = 1; end
    if ~exist('N_runs','var') || isempty(N_runs), N_runs = 1e4; end
    if ~exist('red_fun','var') || isempty(red_fun), red_fun = "MMI"; end
    if ~exist('metric','var') || isempty(metric), metric = ['y','n']; end

    if red_fun == "MMI"
        PID_red = ['y', 'n', 'n'];
        r=1;
    elseif red_fun == "Idep"
        PID_red = ['n', 'y', 'n'];
        r=2;
    elseif red_fun == "Iccs"
        PID_red = ['n', 'n', 'y'];
        r=3;
    end
    
    opts = optimset('display', 'none');
    C = cell(1,p); gs_cell = cell(1,p);
    cc = 0;
    
    % initial conditions for the optimisation
    x0 = [0 -1 2 -2 5 -5];

    Syn = NaN(1,N_runs); Red = NaN(1,N_runs);
    UnX = NaN(1,N_runs); UnY = NaN(1,N_runs);

    phi_WMS = NaN(1,N_runs); phi_R = NaN(1,N_runs);
    red_phi = NaN(1,N_runs);

%     if MI>43, fprintf("MI is "+sprintf('%02d',MI)+"\n"); end
    
    for j = 1:N_runs

%         if mod(j,N_runs/10)==0, disp(j); end
        % sample the matrices
        % NB. for Wishart degrees of freedom must be > size - 1
        
        % random VAR matrix
        A = var_rand(S,p,rand(1));
        while any(abs(A(:))<1e-18)
            A = var_rand(S,p,rand(1));
        end
        
        % target conditional covariance
        V = wishrnd(eye(S),S+1);
%         A_old = A; V_old=V;
%         [A,V] = var_normalise(A,V);
        
        % find g such that mutual information is MI_tot
        l = 1;
%         t0 = tic();
        [g, ~, code, out] = fzero(@(x) fun(x,A,V,MI), x0(l), opts);
          
        while( code ~= 1 && l+1 <= length(x0) )
            l = l+1;
            [g, ~, code, out] = fzero(@(x) fun(x,A,V,MI), x0(l), opts);
        end
        if( code ~= 1 ), cc = cc+1; fprintf("optimisation failed "+sprintf('%01d',code)+"\n"); 
                         continue; end 

        % now do one more calculation
        gg = Sigmoid(g,1);
        [B,~] = specnorm(A,gg);
        for s = 1:p
            C{s} = B(:,:,s);
        end
%         C = num2cell(B,[1,2]);
%         C = squeeze(C);
        if ~isfinite(fun(g,B,V,MI))
            cc = cc+1;
            disp("Optimisation function returned NaN value!");
            continue;
        end

%         fprintf('%f\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t', MI, p, l, g, gg, out.iterations, ...
%            out.intervaliterations, out.funcCount, cc); toc(t0);

        if metric(1)=='y'
            PID = PIDcalc_Mult(p,C,V,S/2,S/2,"joint",PID_red);
            
            try
                PIDcheck(PID);
                assert(isreal(PID), "Not real");
            catch 
                cc = cc+1;
                disp("PID values not real or not consistent!");
                continue;
            end

            MI_calc = sum(PID(r,:));

            if(abs(MI_calc-MI)>MI/1e4), cc=cc+1; 
                disp("Mismatch between true MI and MI calculated during the optimisation!")
                disp(abs(MI_calc-MI)); %disp(MI); 
                continue; 
            end
                    
            Syn(j) = PID(r,4);    Red(j) = PID(r,3);
            UnX(j) = PID(r,1);    UnY(j) = PID(r,2);

        end
        
        if metric(2)=='y'
            gs = var_to_autocov(B,V,p);
            MI_calc = (0.5*log(det(gs(:,:,1)))-0.5*log(det(V)))/log(2);

            for s = 1:p+1
                gs_cell{s} = gs(:,:,s);
            end

            [phi_WMS(j), phi_R(j)] = integrated_info_VAR(MI_calc,gs_cell);
            red_phi(j) =  phi_R(j)-phi_WMS(j);
        end
  
    end

    PIDs = [UnX; UnY; Red; Syn];
    Phis = [phi_WMS; phi_R; red_phi];
%     fprintf("done! excluded values were " + sprintf('%01d',cc)+"\n");
end

function y = fun(x,A,V,I)

    x = Sigmoid(x,1);
    [BB,~] = specnorm(A,x);

    % sometimes for large p and specific (small) values of the spectral radius 
    % the Lyapunov equation does not have a unique solution
    try
        G = var_to_autocov(BB,V,0);
        MI_value = (0.5*log(det(G(:,:,1)))-0.5*log(det(V)))/log(2);
        if ~isreal(MI_value)
            % when the spectral radius is close to 1 the MI can become complex
            error("error in optimising the g!");
        end
        y = MI_value-I;

    catch 
        y=NaN;
        disp("NAN value encountered");

    end
%     disp(y);

end

function y = Sigmoid(x,M)
    
    y = M/(1+exp(-x));

end