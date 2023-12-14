function [time_series, subject_info] = PsychMEGLoad(subject_info)

if subject_info.drug == "KET"
    % subject ids for KETAMINE
    IDs=["240107_6", "310713_51" ,"270407_1", "210813_51", "210813_52", ...
         "280813_51", "280813_52", "110913_51", "021013_51" ,"301112_51", ...
         "091013_51", "161013_51" ,"161013_52", "231013_51", "221112_50", ...
         "171212_1" ,"181113_1", "040806_1", "041213_1"];
elseif subject_info.drug == "PSIL"
    % subject ids for PSILOCYBIN
    IDs=["231109_1","160312_2","190311_1","020311_50","081211_51","050312_1", ...
         "020312_50","010312_1","240107_6","241111_2","020312_51","271108_3", ...
         "021211_51","030507_1"];
elseif subject_info.drug == "LSD"
    % subject ids for LSD
    IDs=["260614_1","010514_1","310714_1","040914_1","270613_1","010813_4", ...
         "040914_2","290514_2","070814_2","090714_2","140514_1","140514_2", ...
         "200814_1","230911_1","041213_1"];
end

file = "../psychedelics_MEG_600Hz/"+subject_info.drug+"/"+IDs(subject_info.name);

if subject_info.condition == "Placebo"
    file = file + "_PLA.mat";   
else 
    file = file + "_" + subject_info.drug + ".mat";  
end

subject_info.fs = 600;
subject_info.ID = IDs(subject_info.name);

data = load(file).dat;
time_series = cell(1,length(data));
for n = 1:length(data)
    time_series{n} = data{n}';
end