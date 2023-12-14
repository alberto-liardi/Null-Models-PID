%% Human vs Monkey analyses 

% load data
subjects = {"Kin2", "George", "Su", "Chibi", 2, 4, 5, 8, 11, 12, 13, 15, 16};
species = ['M', 'M', 'M', 'M', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
red_fun = ["MMI", "Idep", "Iccs"];      conditions = ["W", "S"];
channels = 2;

corr_W = cell(1,2);     corr_S = cell(1,2);

for sub = 1:length(subjects)
    for cond = 1:length(conditions)

        path = createPath(species(sub), subjects{sub}, conditions(cond), channels);
        
        InfoDyn = readmatrix(path+"InfoDyn.txt");
        corr = mean(InfoDyn(:,11:12));

        if cond==1    
            for p = 1:2
                corr_W{p}(sub) = corr(p);
            end
        elseif cond==2
            for p = 1:2
                corr_S{p}(sub) = corr(p);
            end
        end

    end
end

Delta_corr = cell(1,2);

for p = 1:2
    Delta_corr{p} = corr_W{p} - corr_S{p};
end

% finally do the plots

% var(1) awake / sleep-sedated / awake-sleep/sedated
HMViolinPlot_small(corr_W{1}(1:4), corr_S{1}(1:4), Delta_corr{1}(1:4), ...
                   corr_W{1}(5:end), corr_S{1}(5:end), Delta_corr{1}(5:end), ...
                   "Results/HMAnalyses/Corr_var(1).png", ...
                   "Humans and monkeys Residual Correlation - Var(1)", ...
                   ["Awake", "Sleep", "W-S"]);

% var(1) awake / sleep-sedated / awake-sleep/sedated
HMViolinPlot_small(corr_W{2}(1:4), corr_S{2}(1:4), Delta_corr{2}(1:4), ...
                   corr_W{2}(5:end), corr_S{2}(5:end), Delta_corr{2}(5:end), ...
                   "Results/HMAnalyses/Corr_var(p).png", ...
                   "Humans and monkeys Residual Correlation - Var(p)", ...
                   ["Awake", "Sleep", "W-S"]);

function path = createPath(spec, subj, con, chan)

    % creating the name of the path
    if spec == 'H'
        pref = 'S';
        name = sprintf('%01d',subj);
    else
        pref = '';  
        name = subj;
    end
    
    if con == 'S', state = 'Sleep';
    else,    state = 'Awake'; end
    
    path = strcat('Results/', pref, name, '_c', ...
                  sprintf('%01d',chan), '/', state, '/');

end