%rng('default') % for reproducibility of random variables
javaaddpath("infodynamics.jar");

disp(" "); disp("---- Start of compilation ----"); disp(" ");

% subject_info.species = 'H';
% subject_info.condition = 'S';
% subject_info.name = 4;

% subjects = {"Kin2", "George", "Su", "Chibi", 2, 4, 5, 8, 11, 12, 13, 15, 16};
% species = ['M', 'M', 'M', 'M', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
% subjects = {"Kin2", "Su", "George", "Chibi"}; species = ['M', 'M', 'M', 'M']; 
% subjects = {"George", "Su", "Chibi", 2, 4, 5, 8, 11, 12, 13, 15, 16};
% species = ['M', 'M', 'M', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];
channels = [10];     conditions = ["Sed", "Sleep", "W"];

subjects = {8, 11, 12, 13, 15, 16};
species = ['H', 'H', 'H', 'H', 'H', 'H'];

% subjects = {2, 4, 5, 8, 11, 12, 13, 15, 16};
% species = ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'];

% subjects = {4}; conditions = ["S", "W"]; species = ['H'];

assert( length(subjects) == length(species), ...
    "check consistency between subjects and species");

for sub = 1:length(subjects)
    subject_info.species = species(sub);
    subject_info.name = subjects{sub};

    for cond = 1:length(conditions)
        subject_info.condition = conditions(cond);
        if (subject_info.species == "H" && subject_info.condition == "Sed")
            continue
        elseif subject_info.species == "M"
            if (subject_info.name == "Kin2" || subject_info.name == "Su") && ...
            subject_info.condition == "Sleep"
                continue
            end
        end

        for chan = 1:length(channels)
            parameters = struct('channels', channels(chan), 'runs', 100, ...
                'mmorder', 1, 'Red_fun', ['y', 'n', 'n']);
            
            try
                wrapTimeSeries(subject_info, parameters);
            catch
                parameters.channels = parameters.channels-4;
                wrapTimeSeries(subject_info, parameters);
            end
        end
    end
end


disp(" "); disp("---- End of compilation ----"); disp(" ");