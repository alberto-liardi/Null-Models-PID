function [ data ] = SleepLoader(subject, condition, data_folder)

sampling_rate = 1000;

if strcmp(subject, 'Chibi')
  if strcmp(condition, 'W')
    start_sample = 2124.09*sampling_rate;
    end_sample = 2747.71*sampling_rate;
    sess = 'Session2';
  elseif strcmp(condition, 'S')
    start_sample = 75.44*sampling_rate;
    end_sample = 3000.50*sampling_rate;    
    sess = 'Session1';
  else
    error(sprintf('Unknown condition %s for anaesthesia dataset', condition));
  end

elseif strcmp(subject, 'George')
  if strcmp(condition, 'W')
    start_sample = 3382.25*sampling_rate;
    end_sample = 3982.85*sampling_rate;
    sess = 'Session2';
  elseif strcmp(condition, 'S')
    start_sample = 242.64*sampling_rate;
    end_sample = 2819.22*sampling_rate;   
    sess = 'Session1';
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


