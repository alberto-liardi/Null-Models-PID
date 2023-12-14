rng('default') % for reproducibility of random variables
javaaddpath("infodynamics.jar");

disp(" "); disp("---- Start of compilation ----"); disp(" ");

subject_info.species = 'H';
subject_info.drug = "PSIL";

if subject_info.drug=="KET"
    subjects = 1:19;
elseif subject_info.drug=="LSD"
    subjects = 1:15;
elseif subject_info.drug=="PSIL"
    subjects = 1:14;
end

subjects = 3:14;
channels = [10];     conditions = [subject_info.drug, "Placebo"];

for sub = 1:length(subjects)
    subject_info.name = subjects(sub);

    for cond = 1:length(conditions)
        subject_info.condition = conditions(cond);

        for chan = 1:length(channels)
            parameters = struct('channels', channels(chan), 'runs', 100, ...
                'mmorder', 1, 'Red_fun', ['y', 'y', 'y']);
            
            fprintf("doing subject %d condition %s and channels %d", ...
                    subject_info.name,subject_info.condition,channels(chan));
            wrapTimeSeries(subject_info, parameters);

        end
    end
end

disp(" "); disp("---- End of compilation ----"); disp(" ");