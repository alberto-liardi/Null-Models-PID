%% Human analysis: Drug vs Placebo

drug = "PSIL";
Norm_method = "None";
chan = 10;
ps = 1; % only useful for the null model
red_fun = ["MMI", "Idep", "Iccs"];
% red_fun = ["MMI"];

disp("begin run");

% load data
if drug=="KET"
    subjects = 1:19;
elseif drug=="LSD"
    subjects = 1:15;
elseif drug=="PSIL"
    subjects = 1:14;
end

subjects=3:14;

if Norm_method=="None"
    Norm_fun = @identity;
    save_app = "";
    tit_app = "";
elseif Norm_method=="MI"
    Norm_fun = @MI_norm;
    save_app = "_MI_norm";
    tit_app = "MI norm";
elseif Norm_method=="Null"
    Norm_fun = @Quant_norm;
    save_app = "_Null_norm";
    tit_app = "Null norm";
else
    error("Invalid normalisation method!");
end
    
conditions = [drug, "Placebo"];

for r = 1:length(red_fun)
    fprintf("\nUsing %s definition \n", red_fun(r));
    for sub = 1:length(subjects)
        fprintf("\ndoing subject number %d \n", sub)
        for cond = 1:length(conditions)
            fprintf("\ndoing condition %s \n", conditions(cond))

            % loading data
            path = createPath(subjects(sub), drug, conditions(cond), red_fun(r), chan);
            [UnX, UnY, Red, Syn, MI] = MonkeyPIDLoader(path + "data.txt");
            
            % using normalisation method (if any)
            [UnX, UnY, Red, Syn, MI] = Norm_fun(path, chan, UnX, UnY, Red, ...
                                                Syn, MI, ps, red_fun(r));        

            if cond==1
                UnX_W = UnX{1}; UnY_W = UnY{1};
                Red_W = Red{1}; Syn_W = Syn{1};
                MI_W = MI{1};
            elseif cond==2
                UnX_S = UnX{1}; UnY_S = UnY{1};
                Red_S = Red{1}; Syn_S = Syn{1};
                MI_S = MI{1};
            end
        end

        % finally do the plots
%         HViolinPlot(UnX_W, UnY_W, Red_W, Syn_W, ...
%                      MI_W, red_fun(r) + ": Drug Var(1) " + ...
%                      sprintf('%01d',chan)+" channels "+tit_app, "Results_psychedelics_all_redfun/"+drug+...
%                      "/H_"+red_fun(r)+"_Var(1)_c"+...
%                      sprintf('%01d',chan)+save_app+"_drug.png")
%         
%         HViolinPlot(UnX_S, UnY_S, Red_S, Syn_S, ...
%                      MI_S, red_fun(r) + ": Placebo Var(1) " + ...
%                      sprintf('%01d',chan)+" channels "+tit_app, "Results_psychedelics_all_redfun/"+drug+...
%                      "/H_"+red_fun(r)+"_Var(1)_c"+...
%                      sprintf('%01d',chan)+save_app+"_placebo.png")

        PIDs = [MI_W', UnX_W', UnY_W', Red_W', Syn_W'];
        T = array2table(PIDs);
        T.Properties.VariableNames(1:5) = {'MI','Unique X','Unique Y','Redundancy','Synergy'};
        out_folder = "../Null_model_figures/"+drug+"_all_redfun/S"+sprintf('%01d',subjects(sub))+"_D/";
        if ~exist(out_folder, 'dir')
            mkdir(out_folder);
        end
        writetable(T,out_folder+red_fun(r)+"_VAR1_C10"+save_app+".csv");
    
        PIDs = [MI_S', UnX_S', UnY_S', Red_S', Syn_S'];
        T = array2table(PIDs);
        T.Properties.VariableNames(1:5) = {'MI','Unique X','Unique Y','Redundancy','Synergy'};
        out_folder = "../Null_model_figures/"+drug+"_all_redfun/S"+sprintf('%01d',subjects(sub))+"_P/";
        if ~exist(out_folder, 'dir')
            mkdir(out_folder);
        end
        writetable(T,out_folder+red_fun(r)+"_VAR1_C10"+save_app+".csv");
    
    end
end

disp("Done execution!");

function path = createPath(subj, drug, state, red, chan)

    % creating the name of the path
    name = sprintf('%01d',subj);
    path = strcat('Results_psychedelics_all_redfun/', drug, '/', name, '_c', ...
                  sprintf('%01d',chan), '/', state, '/', red, '/');

end


function [UnX, UnY, Red, Syn, MI] = identity(~, ~, UnX, UnY, Red, Syn, MI, ~, ~)

    % doing nothing
end 

function [UnX, UnY, Red, Syn, MI] = MI_norm(~, ~, UnX, UnY, Red, Syn, MI, ~, ~)
    
    [UnX, UnY, Red, Syn, MI] = NormPID(UnX, UnY, Red, Syn, MI);

end

function [UnX, UnY, Red, Syn, MI] = Quant_norm(path, chan, UnX, UnY, Red, Syn, MI, ks, red_fun)
    % Instead of normalising w.r.t. MI, use the Null distribution!

    app = ["_1","_p"];
    for k = ks
        name=path+"Quantiles"+app(k)+".mat";
        % if quantiles have been previosuly computed, load them 
        if isfile(name)
            [UnX{k}, UnY{k}, Red{k}, Syn{k}, MI{k}] = LoadQuantPID(name);
%         
%         % otherwise compute the quantiles with the null distribution
        else
            % load model orders
            InfoDyn = readmatrix(path+"../InfoDyn.txt");
                
            if k==1
                p = InfoDyn(:,7); % array of ones...
            elseif k==2
                p = InfoDyn(:,8);
            end

            %t0 = tic();
            model = struct('name',"VAR",'S',chan,'p',p,'red_fun',red_fun);
            model.n = 100;
            [UnX{k}, UnY{k}, Red{k}, Syn{k}] = ...
                QuantPID(UnX{k}, UnY{k}, Red{k}, Syn{k}, MI{k}, model);
            %fprintf("%d \t", ps); toc(t0);
        
            % save quantiles on file
            SaveQuantPID(UnX{k}, UnY{k}, Red{k}, Syn{k}, MI{k}, path, k);
        end

        % if k is only one model order, set the other one's PID quantities to zero
        % TODO: maybe it's better to set it to NaN? need to check if it throws errors later
        if k==1
            [UnX{2}, UnY{2}, Red{2}, Syn{2}, MI{2}] = deal(zeros(1,length(UnX{2})));
        elseif k==2
            [UnX{1}, UnY{1}, Red{1}, Syn{1}, MI{1}] = deal(zeros(1,length(UnX{1})));
        end

        % set the mutual information to zero as they do not matter anymore
        [MI{1}, MI{2}] = deal(zeros(1,length(MI{1})));
    
    end
    
end