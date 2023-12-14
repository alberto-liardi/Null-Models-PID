%% Human vs Monkey analyses: comparison between the phi-measures

% load data for every subject
subjects = {"Kin2", "George", "Su", "Chibi", 2, 4, 5, 8, 11, 12, 13, 15, 16};
species = ['M', 'M', 'M', 'M', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
conditions = ["W", "S"];
channels = [2,6,10];

for chan = 1:length(channels)
    
    phi_wms_w = []; phi_wms_s = [];
    phi_r_w = []; phi_r_s = [];
    red_w = []; red_s = [];
    for sub = 1:length(subjects)
        for cond = 1:length(conditions)
    
            % loading the data (phis)
            path = createPath(species(sub), subjects{sub}, conditions(cond), channels(chan));
            InfoDyn = readmatrix(path+"InfoDyn.txt");

%             phi_wms = mean(InfoDyn(:,13));
%             phi_r = mean(InfoDyn(:,14));
%             red = phi_r - phi_wms;

            [phi_wms,phi_r,red]=normalise_phi(InfoDyn(:,13)',InfoDyn(:,14)',path);
            phi_wms = mean(phi_wms,2);
            phi_r = mean(phi_r,2);
            red = mean(red,2);
    
            if cond==1
                phi_wms_w(end+1) = phi_wms;
                phi_r_w(end+1) = phi_r;
                red_w(end+1) = red;
            elseif cond==2
                phi_wms_s(end+1) = phi_wms;
                phi_r_s(end+1) = phi_r;
                red_s(end+1) = red;
            end
    
        end
    
    end

    delta_phi_wms = phi_wms_w - phi_wms_s;
    delta_phi_r = phi_r_w - phi_r_s;
    delta_r = red_w - red_s;
    
    % creating the plots
    HMViolinPlot_small(delta_phi_wms(1:4), delta_phi_r(1:4), delta_r(1:4), ...
                 delta_phi_wms(5:end), delta_phi_r(5:end), delta_r(5:end), ...
                 "Results_p1_phi/HMAnalyses/HM_phi_c"+sprintf('%01d',channels(chan))+"_norm.png", ...
                 "Normalised W-S $\phi$ Human and Monkeys - "+sprintf('%01d',channels(chan))+" channels", ...
                 ["\phi_{WMS}", "\phi_R", "Red (MMI)"]);

end

%% Human vs Monkey analyses: quantiles comparison

disp("begin run");

% load data
subjects = {"Kin2", "George", "Su", "Chibi", 2, 4, 5, 8, 11, 12, 13, 15, 16};
species = ['M', 'M', 'M', 'M', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
% subjects = {2, 4, 5, 8, 11, 12, 13, 15, 16};
% species = ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
conditions = ["W", "S"];
channels = [2,6,10];

for c = 1:length(channels)
    chan = channels(c);
    delta_phi_wms=[]; delta_phi_r=[]; delta_red_phi=[];

    for sub = 1:length(subjects)
        fprintf("Doing subject "+subjects(sub)+"\n"); 
        for cond = 1:length(conditions)
    
            if chan == 20
                if sub>4 && (subjects{sub}==4 || subjects{sub}==5 || subjects{sub}==13)
                    chan = 16;
                elseif sub>4 && subjects{sub}==16
                    chan = 12;
                else
                    chan = 20;
                end
            end
    
            % load data (PIDs and model order)
            path = createPath(species(sub), subjects{sub}, conditions(cond), chan);
            InfoDyn = readmatrix(path+"InfoDyn.txt");
            [~,~,~,~,MI] = MonkeyPIDLoader(path+"/MMI/data.txt");
    
            phi_wms = InfoDyn(:,13)';
            phi_r = InfoDyn(:,14)';
            red = phi_r - phi_wms;

            % Instead of normalising w.r.t. MI, use the Null distribution!
            model = struct('name',"VAR",'S',chan,'metric',['n','y']);
            model.n = 100;
            [phi_wms_q, phi_r_q, red_q] = ...
                QuantPhi(phi_wms, phi_r, red, MI{1}, model);
    
            if cond==1
                [phi_wms_w, phi_r_w, red_w] = meanPhi(phi_wms_q, phi_r_q, red_q);
            elseif cond==2
                [phi_wms_s, phi_r_s, red_s] = meanPhi(phi_wms_q, phi_r_q, red_q);
            end
    
        end
         
        delta_phi_wms(end+1) = phi_wms_w - phi_wms_s;
        delta_phi_r(end+1) = phi_r_w - phi_r_s;
        delta_red_phi(end+1) = red_w - red_s;
    
    end

    % creating the plots
    HMViolinPlot_small(delta_phi_wms(1:4), delta_phi_r(1:4), delta_red_phi(1:4), ...
                 delta_phi_wms(5:end), delta_phi_r(5:end), delta_red_phi(5:end), ...
                 "Results_p1_phi/HMAnalyses/HM_phi_c"+sprintf('%01d',chan)+"_quant.png", ...
                 "Quantiles W-S $\phi$ Human and Monkeys - "+sprintf('%01d',chan)+" channels", ...
                 ["\phi_{WMS}", "\phi_R", "Red (MMI)"]);
    
end


function [norm_phi_wms,norm_phi_r,norm_red] = normalise_phi(phi_wms,phi_r,dir)

    dir = dir+"MMI/data.txt";
    [~,~,~,~,MI] = MonkeyPIDLoader(dir);

    norm_phi_wms = phi_wms./MI{1};
    norm_phi_r = phi_r./MI{1};
    norm_red = (phi_r-phi_wms)./MI{1};


end 

function [Y1, Y2, Y3] = meanPhi(X1, X2, X3)

    Y1 = mean(X1,2);
    Y2 = mean(X2,2);
    Y3 = mean(X3,2);

end

function path = createPath(spec, subj, con, chan)

    % creating the name of the path
    if spec == 'H'
        pref = 'S';
        name = sprintf('%01d',subj);
    else
        pref = '';  
        name = subj;
    end
    
    if con == 'S' && spec == 'M', state = 'Sedated';
    elseif con ~= 'W', state = 'Sleep';
    else,    state = 'Awake'; 
    end
    
    path = strcat('Results_p1_phi/', pref, name, '_c', ...
                  sprintf('%01d',chan), '/', state, '/');

end

