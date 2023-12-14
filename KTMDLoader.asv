function [ data ] = KTMDLoader(subject, condition, data_folder)

sampling_rate = 1000;

if strcmp(subject, 'Chibi')
  if strcmp(condition, 'W')
    start_sample = 1478.51*sampling_rate;
    end_sample = 2685.72*sampling_rate;
    sess = 'Session1';
  elseif strcmp(condition, 'S')
    start_sample = 699.73*sampling_rate;
    end_sample = 2498.80*sampling_rate;    
    sess = 'Session2';
  else
    error(sprintf('Unknown condition %s for anaesthesia dataset', condition));
  end

elseif strcmp(subject, 'George')
  if strcmp(condition, 'W')
    start_sample = 695.39*sampling_rate;
    end_sample = 1296.33*sampling_rate;
    sess = 'Session1';
  elseif strcmp(condition, 'S')
    start_sample = 1798.73*sampling_rate;
    end_sample = 3051.68*sampling_rate;   
    sess = 'Session1';
  else
    error(sprintf('Unknown condition %s for anaesthesia dataset', condition));
  end

elseif strcmp(subject, 'Kin2')
  if strcmp(condition, 'W')
    start_sample = 23.76*sampling_rate;
    end_sample = 1239.97*sampling_rate;
    sess = 'Session2';
  elseif strcmp(condition, 'S')
    start_sample = 585.62*sampling_rate;
    end_sample = 2396.60*sampling_rate;   
    sess = 'Session3';
  else
    error(sprintf('Unknown condition %s for anaesthesia dataset', condition));
  end

elseif strcmp(subject, 'Su')
  if strcmp(condition, 'W')
    start_sample = 1601.33*sampling_rate;
    end_sample = 2922.45*sampling_rate;
    sess = 'Session1';
  elseif strcmp(condition, 'S')
    start_sample = 925.06*sampling_rate;
    end_sample = 2862.40*sampling_rate;   
    sess = 'Session2';
  else
    error(sprintf('Unknown condition %s for anaesthesia dataset', condition));
  end

else
  error(sprintf('Unknown subject %s for anaesthesia dataset', subject));
end

data_path = sprintf('%s/%s/%s', data_folder, subject, sess);
fl = dir([data_path, '/ECoG_ch*']);
nb_channels = length(fl);
f = load([data_path, '/ECoG_ch1.mat']);
nb_samples  = length(f.ECoGData_ch1);
M = zeros([nb_samples, nb_channels]);

for i=1:length(fl)
  f = load(sprintf('%s/%s.mat', data_path, ['ECoG_ch', num2str(i)]));
  M(:,i) = f.(['ECoGData_ch', num2str(i)]);
end

M = M(start_sample:end_sample-1, :);


ws = 1200;
i = 1;
data = cell([1, floor(length(M)/ws)]);
while i*ws < length(M)
  data{i} = M((i-1)*ws+1:i*ws, :);
  i = i+1;
end


