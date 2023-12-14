function [time_series, subject_info] = ValdasLoad(subject_info)

condition = subject_info.condition;
subject = subject_info.name;

if condition == "Sleep"
    csv = readtable('ValdasPreprocessed/MAIN data file.csv');
    for x = 1:size(csv,1)
        if csv{x,7}==subject && csv{x,5}==2
            ID = csv{x,1};
            Nights = csv{x,8};
            Awakenings = csv{x,9};
            break
        end
    end

    filename = strcat('ID_', sprintf('%03d',ID), '_S', sprintf('%01d',subject), ...
           '_N', sprintf('%01d',Nights), '_A', sprintf('%01d',Awakenings), '.set');

    subject_info.ID = ID;
    subject_info.Nights = Nights;
    subject_info.Awakenings = Awakenings;

elseif condition == 'W'
    if subject == 2 
        filename = 'baseline 1.set';
    elseif subject == 4
        filename = 'baseline 2.set';
    elseif subject == 5
        filename = 'baseline 3.set';
    elseif subject == 8
        filename = 'baseline 4.set';
    elseif subject == 11
        filename = 'baseline 5.set';
    elseif subject == 12
        filename = 'baseline 6.set';
    elseif subject == 13
        filename = 'baseline 7.set';
    elseif subject == 15
        filename = 'baseline 8.set';
    elseif subject == 16
        filename = 'baseline 9.set';
    end
end

data_folder = 'ValdasPreprocessed';

cfg = [];
cfg.dataset = [data_folder, '/', filename];
cfg.channel = {'all', '-REF'};
data = ft_preprocessing(cfg);

cfg = [];
cfg.length = 5;
data = ft_redefinetrial(cfg, data);

subject_info.fs = data.fsample;

time_series = data.trial;
for n = 1:length(time_series)
    time_series{n} = time_series{n}';
end

