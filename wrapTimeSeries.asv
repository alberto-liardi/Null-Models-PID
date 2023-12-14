function [] = wrapTimeSeries(subject_info, parameters)

% loading data (if necessary)
% if exist('data',     'var') == 0
    disp("Loading data... "+subject_info.name+" "+subject_info.condition+...
         " "+sprintf('%01d',parameters.channels)+" channels");
    
    if subject_info.species == 'M'

        if subject_info.condition == "Sleep"
        % !!! IMPORTANT: Monkey Sleep data must be in a folder called "Monkeys_data_Sleep"
            data = SleepLoader(subject_info.name, "S", "Monkeys_data_Sleep");
        elseif subject_info.condition == "Sed"
            % !!! IMPORTANT: Monkey Sedated data must be in a folder called "Monkeys_data"
            data = KTMDLoader(subject_info.name, "S", "Monkeys_data");
        elseif subject_info.condition == "W"
            % !!! IMPORTANT: Monkey Sedated data must be in a folder called "Monkeys_data"
            data = KTMDLoader(subject_info.name, "W", "Monkeys_data"); 
        end

        subject_info.fs = 1000;

    elseif subject_info.species == 'H'

        if isfield(subject_info, 'drug')
            % !!! IMPORTANT: Human data must be in a folder called "../psychedelics_MEG_600Hz"
            [data, subject_info] = PsychMEGLoad(subject_info);
        else
            % !!! IMPORTANT: Human data must be in a folder called "ValdasPreprocessed"
            [data, subject_info] = ValdasLoad(subject_info);
        end

    end
    disp("... done.");
% end

% start the timer
t0 = tic();

% running analyses
Information = TimeSeriesAnalyses(data, subject_info, parameters);

% stop the timer and print elapsed time 
disp(" "); toc(t0);

% saving everything on file
saveTimeSeriesAnalyses(Information);

% clean the workspace
%clear