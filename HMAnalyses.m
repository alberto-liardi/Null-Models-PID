%% Human vs Monkey analyses 

% load data
% subjects = {"Kin2", "George", "Su", "Chibi", 2, 4, 5, 8, 11, 12, 13, 15, 16};
% species = ['M', 'M', 'M', 'M', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
red_fun = ["MMI"];%, "Idep", "Iccs"];      
conditions = ["W", "S"];
chan = 20;

subjects = {2, 4, 5, 8, 11, 12, 13, 15, 16};
species = ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];

DeltaX = cell(1,2);     DeltaY = cell(1,2);
DeltaS = cell(1,2);     DeltaR = cell(1,2);
DeltaMI = cell(1,2);

r = 1;
% for r = 1:length(red_fun)
    for sub = 1:length(subjects)

        if sub>4 && (subjects{sub}==4 || subjects{sub}==5 || subjects{sub}==13)
%             channels = 2;
            channels = 16;
        elseif sub>4 && subjects{sub}==16
%             channels = 2;
            channels = 12;
        else
%             channels = 2;
            channels = 20;
        end

        for cond = 1:length(conditions)

            path = createPath(species(sub), subjects{sub}, conditions(cond), red_fun(r), channels);
            [UnX, UnY, Red, Syn, MI] = MonkeyPIDLoader(path + "data.txt");
%             [UnX, UnY, Red, Syn, MI] = NormPID(UnX, UnYq, Red, Syn, MI);

            if cond==1
                [UnX_W, UnY_W, Red_W, Syn_W, MI_W] = meanPID(UnX, UnY, Red, Syn, MI);
            elseif cond==2
                [UnX_S, UnY_S, Red_S, Syn_S, MI_S] = meanPID(UnX, UnY, Red, Syn, MI);
            end

        end

        % taking the difference between W and S
        for p = 1:2

%             DeltaX{p}(sub) = (UnX_W(p) - UnX_S(p) )/ ((UnX_W(p) + UnX_S(p))/2)*100;
%             DeltaY{p}(sub) = (UnY_W(p) - UnY_S(p) )/ ((UnY_W(p) + UnY_S(p))/2)*100;
%             DeltaS{p}(sub) = (Syn_W(p) - Syn_S(p) )/ ((Syn_W(p) + Syn_S(p))/2)*100;
%             DeltaR{p}(sub) = (Red_W(p) - Red_S(p) )/ ((Red_W(p) + Red_S(p))/2)*100;
%             DeltaMI{p}(sub) = (MI_W(p) - MI_S(p) )/ ((MI_W(p) + MI_S(p))/2)*100;

            DeltaX{p}(sub) = UnX_W(p) - UnX_S(p);
            DeltaY{p}(sub) = UnY_W(p) - UnY_S(p);
            DeltaS{p}(sub) = Syn_W(p) - Syn_S(p);
            DeltaR{p}(sub) = Red_W(p) - Red_S(p);
            DeltaMI{p}(sub) = MI_W(p) - MI_S(p);


        end

    end

    % finally do the plots

    HMViolinPlot(DeltaX{1}(1:4), DeltaY{1}(1:4), DeltaR{1}(1:4), DeltaS{1}(1:4), ...
                 DeltaMI{1}(1:4), DeltaX{1}(5:end), DeltaY{1}(5:end), DeltaR{1}(5:end), ...
                 DeltaS{1}(5:end), DeltaMI{1}(5:end), red_fun(r) + ": W-S Var(1) Human and Monkeys - " + ...
                 sprintf('%01d',chan)+" channels", "Results/HMAnalyses/HM_"+red_fun(r)+"_Var(1)_c"+...
                 sprintf('%01d',chan)+".png")

   HMViolinPlot(DeltaX{p}(1:4), DeltaY{p}(1:4), DeltaR{p}(1:4), DeltaS{p}(1:4), ...
                 DeltaMI{p}(1:4), DeltaX{p}(5:end), DeltaY{p}(5:end), DeltaR{p}(5:end), ...
                 DeltaS{p}(5:end), DeltaMI{p}(5:end), red_fun(r) + ": W-S Var(p) Human and Monkeys - " + ...
                 sprintf('%01d',chan)+" channels", "Results/HMAnalyses/HM_"+red_fun(r)+"_Var(p)_c"+...
                 sprintf('%01d',chan)+".png")

% end

%% now adding Monkey Sleep

subjects = {"George", "Chibi"};
species = ['M', 'M'];
red_fun = ["MMI"];%, "Idep", "Iccs"];      
conditions = ["W", "Sleep"];
channels = 20;

% for r = 1:length(red_fun)
for sub = 1:length(subjects)

    for cond = 1:length(conditions)

        path = createPath(species(sub), subjects{sub}, conditions(cond), red_fun(r), channels);
        [UnX, UnY, Red, Syn, MI] = MonkeyPIDLoader(path + "data.txt");
        [UnX, UnY, Red, Syn, MI] = NormPID(UnX, UnY, Red, Syn, MI);

        if cond==1
            [UnX_W, UnY_W, Red_W, Syn_W, MI_W] = meanPID(UnX, UnY, Red, Syn, MI);
        elseif cond==2
            [UnX_S, UnY_S, Red_S, Syn_S, MI_S] = meanPID(UnX, UnY, Red, Syn, MI);
        end

    end

    for p = 1:2
    
        DeltaX{p}(end+1) = UnX_W(p) - UnX_S(p);
        DeltaY{p}(end+1) = UnY_W(p) - UnY_S(p);
        DeltaS{p}(end+1) = Syn_W(p) - Syn_S(p);
        DeltaR{p}(end+1) = Red_W(p) - Red_S(p);
        DeltaMI{p}(end+1) = MI_W(p) - MI_S(p);
    
    end
end
% end

HMViolinPlot_big(DeltaX{1}(1:4), DeltaY{1}(1:4), DeltaR{1}(1:4), DeltaS{1}(1:4), ...
                 DeltaMI{1}(1:4), DeltaX{1}(5:end-2), DeltaY{1}(5:end-2), DeltaR{1}(5:end-2), ...
                 DeltaS{1}(5:end-2), DeltaMI{1}(5:end-2), DeltaX{1}(end-1:end), DeltaY{1}(end-1:end), ... 
                 DeltaR{1}(end-1:end), DeltaS{1}(end-1:end), DeltaMI{1}(end-1:end), ... 
                 red_fun(r) + ": Norm W-S Var(1) Human and Monkeys - " + ...
                 sprintf('%01d',20)+" channels", "Results/HMAnalyses/HM_"+red_fun(r)+"_Var(1)_c20_norm.png")

HMViolinPlot_big(DeltaX{2}(1:4), DeltaY{2}(1:4), DeltaR{2}(1:4), DeltaS{2}(1:4), ...
                 DeltaMI{2}(1:4), DeltaX{2}(5:end-2), DeltaY{2}(5:end-2), DeltaR{2}(5:end-2), ...
                 DeltaS{2}(5:end-2), DeltaMI{2}(5:end-2), DeltaX{2}(end-1:end), DeltaY{2}(end-1:end), ... 
                 DeltaR{2}(end-1:end), DeltaS{2}(end-1:end), DeltaMI{2}(end-1:end), ... 
                 red_fun(r) + ": Norm W-S Var(p) Human and Monkeys - " + ...
                 sprintf('%01d',20)+" channels", "Results/HMAnalyses/HM_"+red_fun(r)+"_Var(p)_c20_norm.png")





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
    
    path = strcat('Results/', pref, name, '_c', ...
                  sprintf('%01d',chan), '/', state, '/', red, '/');

end

