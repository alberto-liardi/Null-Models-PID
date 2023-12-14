subjects = {2}; species = ['H']; 
conditions = ['S'];

subject_info.name = subjects{1};
subject_info.condition = conditions;
subject_info.species = species;

[data, subject_info] = ValdasLoad(subject_info);
% data = KTMDLoader(subject_info.name, subject_info.condition, "Monkeys_data");
% subject_info.fs = 1000;

tot_channels = size(data{1},2);
c = randperm(tot_channels, 2);

tot_epochs = length(data);
e = 1:tot_epochs-1;
% e = tot_epochs;

tot_period = size(data{1},1);

% building time series
Z = zeros(2, tot_period, tot_epochs-1);
for j = 1:2
    for n = 1:tot_epochs-1
        Z(j,:,n) = data{e(n)}(:,c(j))';
    end
end

% filtering
down_s = subject_info.fs/200;
Y = zeros(2, tot_period/down_s, tot_epochs-1);
for j = 1:2
    for n = 1:tot_epochs-1
        Y(j,:,n) = LowFilter(Z(j,:,n), down_s, subject_info.fs);
    end
end
Z = Y;

[A,p,V,~] = EvolCoefFit(Z, 30, false);

    