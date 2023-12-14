%% Human vs Monkey analyses with the PID null model

% rng('default') % for reproducibility of random variables
disp("begin run");

% load data
% subjects = {"Kin2", "George", "Su", "Chibi", 2, 4, 5, 8, 11, 12, 13, 15, 16};
% species = ['M', 'M', 'M', 'M', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
subjects = {2, 4, 5, 8, 11, 12, 13, 15, 16};
species = ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
red_fun = ["MMI"];
conditions = ["W", "S"];
channels = [10];
ps=1;   


DeltaX = cell(1,2);     DeltaY = cell(1,2);
DeltaS = cell(1,2);     DeltaR = cell(1,2);
DeltaMI = cell(1,2);

for r = 1;length(red_fun)
    for c = 1:length(channels)
        chan = channels(c);
        for sub = 1:length(subjects)
            fprintf("Doing subject "+subjects(sub)+"\n"); 
            for cond = 1:length(conditions)
                fprintf("\ndoing condition %s \n", conditions(cond))
        
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
                path = createPath(species(sub), subjects{sub}, conditions(cond), red_fun(r), chan);
                [UnX, UnY, Red, Syn, MI] = MonkeyPIDLoader(path + "data.txt");
        
                % Instead of normalising w.r.t. MI, use the Null distribution!
                [UnX, UnY, Red, Syn, MI] = Quant_norm(path, chan, UnX, UnY, Red, ...
                                                    Syn, MI, ps, red_fun(r));
                
                if cond==1
                    [UnX_W, UnY_W, Red_W, Syn_W, MI_W] = meanPID(UnX, UnY, Red, Syn, MI);
                elseif cond==2
                    [UnX_S, UnY_S, Red_S, Syn_S, MI_S] = meanPID(UnX, UnY, Red, Syn, MI);
                end
        
            end
            
    
            % taking the difference between W and S
            for p = 1:2
        
                DeltaX{p}(sub) = UnX_W(p) - UnX_S(p);
                DeltaY{p}(sub) = UnY_W(p) - UnY_S(p);
                DeltaS{p}(sub) = Syn_W(p) - Syn_S(p);
                DeltaR{p}(sub) = Red_W(p) - Red_S(p);
                DeltaMI{p}(sub) = 0;
        
            end
        
        end
    end

    disp("Done execution!");
    
    % finally do the plots
    if ps==1
%         HMViolinPlot(DeltaX{1}(1:4), DeltaY{1}(1:4), DeltaR{1}(1:4), DeltaS{1}(1:4), ...
%                      DeltaMI{1}(1:4), DeltaX{1}(5:end), DeltaY{1}(5:end), DeltaR{1}(5:end), ...
%                      DeltaS{1}(5:end), DeltaMI{1}(5:end), red_fun(r) + ": W-S Var(1) Human and Monkeys - " + ...
%                      sprintf('%01d',chan)+" channels", "Results_p10/HMAnalyses/HM_"+red_fun(r)+"_Var(1)_c" + ...
%                      sprintf('%01d',chan)+"_quant.png");

          HViolinPlot(DeltaX{1}, DeltaY{1}, DeltaR{1}, DeltaS{1}, ...
                 DeltaMI{1}, red_fun(r) + ": W-S Var(1) " + ...
                 sprintf('%01d',chan)+" channels Null norm", "Results_p10/null_model"+...
                 "/H_"+red_fun(r)+"_Var(1)_c"+...
                 sprintf('%01d',chan)+"_null.png")
    elseif ps==2
%         HMViolinPlot(DeltaX{p}(1:4), DeltaY{p}(1:4), DeltaR{p}(1:4), DeltaS{p}(1:4), ...
%                      DeltaMI{p}(1:4), DeltaX{p}(5:end), DeltaY{p}(5:end), DeltaR{p}(5:end), ...
%                      DeltaS{p}(5:end), DeltaMI{p}(5:end), red_fun(r) + ": W-S Var(p) Human and Monkeys - " + ...
%                      sprintf('%01d',chan)+" channels", "Results_p10/HMAnalyses/HM_"+red_fun(r)+"_Var(p)_c" + ...
%                      sprintf('%01d',chan)+"_quant.png");
        HViolinPlot(DeltaX{p}, DeltaY{p}, DeltaR{p}, DeltaS{p}, ...
                 DeltaMI{p}, red_fun(r) + ": W-S Var(p) " + ...
                 sprintf('%01d',chan)+" channels Null norm ", "Results_p10/null_model"+...
                 "/H_"+red_fun(r)+"_Var(p)_c"+...
                 sprintf('%01d',chan)+save_app+"_null.png")
    end
end



function [UnX_C, UnY_C, Red_C, Syn_C, MI_C] = meanPID(UnX, UnY, Red, Syn, MI)
    
    UnX_C = zeros(1,2);     UnY_C = zeros(1,2);
    Syn_C = zeros(1,2);     Red_C = zeros(1,2);
    MI_C = zeros(1,2);

    for p = 1:2
        UnX_C(p) =  mean(UnX{p});
        UnY_C(p) =  mean(UnY{p});
        Syn_C(p) =  mean(Syn{p});
        Red_C(p) =  mean(Red{p});
        MI_C(p) = mean(MI{p});
    end

end


function path = createPath(spec, subj, con, red, chan)

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
    
    path = strcat('Results_p10/', pref, name, '_c', ...
                  sprintf('%01d',chan), '/', state, '/', red, '/');

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

