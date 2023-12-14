% rng('default') % for reproducibility of random variables
rng(82);

% MI_list = [0.01,0.05,0.1,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0];
MI_start = 3; MI_end = 6; n_MI=100;
% MI_list = logspace(log10(MI_start),log10(MI_end),n_MI);
% MI_list = linspace(MI_start,MI_end,n_MI);
MI_list = unifrnd(MI_start,MI_end,1,n_MI);
% max_Syn = 0.1125*MI_list+0.15 +0.6;
% max_Syn(max_Syn>0.98) = 0.98;
PIDs = NaN(4,length(MI_list),3);
quantPID = NaN(4,length(MI_list),3);
S = 10; p=1;
null_PIDs = cell(1,length(MI_list));
dir_out = "../Null_model_figures/Different_PID_study_VAR/high_syn_NMI_rand_S"+...
          sprintf('%d',S)+"p1/"; %high_syn_Null_lin_

for i = 1:length(MI_list)

%     if(mod(i,10)==0), disp(i); end
    disp(i);

    MI = MI_list(i);
    disp(MI);
    
    % construct the null models and get the quantiles
    null_PIDs{i} = MI_null_model_VAR(MI,[],[],1000,"all");

    n=0;
    norm_syn=0;
    while norm_syn<0.75
        % generate a random Gaussian system
        [A,V,check] = VAR_from_MI(MI,p,S);

        if check==1, continue; end
        
        % calculate PID atoms
        PIDs(:,i,:) = PIDcalc_Mult(p,A,V,S/2,S/2,'joint',['y','y','y'])';

        n=n+1;
        if(mod(n,10)==0), disp(n); end

%         norm_syn = comp_quantile(null_PIDs{i}(:,:,1), PIDs(:,i,1));
%         norm_syn = norm_syn(4);
        Syn = PIDs(4,i,1);
        norm_syn = Syn/MI;
%         norm_syn = 100;

        if n>1000, MI=unifrnd(MI_start,MI_end,1); MI_list(i)=MI; fprintf("NEW MI %0.3f",MI); n=0; end 

    end
    disp(norm_syn);
   
    for n = 1:3
%         PIDcheck(PIDs(:,i,n));
        quantPID(:,i,n) = comp_quantile(null_PIDs{i}(:,:,n), PIDs(:,i,n));
    end

end

% save everything on file 

if ~exist(dir_out, 'dir')
       mkdir(dir_out)
end

% raw PID
save22csv([MI_list; PIDs(:,:,1)],dir_out, "PIDs_MMI");
save22csv([MI_list; PIDs(:,:,2)],dir_out, "PIDs_Idep");
save22csv([MI_list; PIDs(:,:,3)],dir_out, "PIDs_Iccs");

% NMI
save22csv([MI_list; PIDs(:,:,1)./MI_list],dir_out, "PIDs_MMI_MI_norm");
save22csv([MI_list; PIDs(:,:,2)./MI_list],dir_out, "PIDs_Idep_MI_norm");
save22csv([MI_list; PIDs(:,:,3)./MI_list],dir_out, "PIDs_Iccs_MI_norm");

% quantiles
save22csv([MI_list; quantPID(:,:,1)],dir_out, "PIDs_MMI_Null_norm");
save22csv([MI_list; quantPID(:,:,2)],dir_out, "PIDs_Idep_Null_norm");
save22csv([MI_list; quantPID(:,:,3)],dir_out, "PIDs_Iccs_Null_norm");



%% null model results comparison for different definitions

% MMI
fig=figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle("Null model PID atoms")
subplot(2,2,1);
scatter(MI_list,quantPID(1,:,1));
hold on
scatter(MI_list,quantPID(1,:,2));
hold on
scatter(MI_list,quantPID(1,:,3));
legend("MMI", "Idep", "Iccs")
title("Unique X")

subplot(2,2,2);
scatter(MI_list,quantPID(2,:,1));
hold on
scatter(MI_list,quantPID(2,:,2));
hold on
scatter(MI_list,quantPID(2,:,3));
legend("MMI", "Idep", "Iccs")
title("Unique Y")

subplot(2,2,3);
scatter(MI_list,quantPID(3,:,1));
hold on
scatter(MI_list,quantPID(3,:,2));
hold on
scatter(MI_list,quantPID(3,:,3));
legend("MMI", "Idep", "Iccs")
title("Redundancy")

subplot(2,2,4);
scatter(MI_list,quantPID(4,:,1));
hold on
scatter(MI_list,quantPID(4,:,2));
hold on
scatter(MI_list,quantPID(4,:,3));
legend("MMI", "Idep", "Iccs")
title("Synergy")

%% NMI results comparison for different definitions


% MMI
fig=figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle("NMI PID atoms")
subplot(2,2,1);
scatter(MI_list,PIDs(1,:,1)/MI_list);
hold on
scatter(MI_list,PIDs(1,:,2)/MI_list);
hold on
scatter(MI_list,PIDs(1,:,3)/MI_list);
legend("MMI", "Idep", "Iccs")
title("Unique X")

subplot(2,2,2);
scatter(MI_list,PIDs(2,:,1)/MI_list);
hold on
scatter(MI_list,PIDs(2,:,2)/MI_list);
hold on
scatter(MI_list,PIDs(2,:,3)/MI_list);
legend("MMI", "Idep", "Iccs")
title("Unique Y")

subplot(2,2,3);
scatter(MI_list,PIDs(3,:,1)/MI_list);
hold on
scatter(MI_list,PIDs(3,:,2)/MI_list);
hold on
scatter(MI_list,PIDs(3,:,3)/MI_list);
legend("MMI", "Idep", "Iccs")
title("Redundancy")

subplot(2,2,4);
scatter(MI_list,PIDs(4,:,1)/MI_list);
hold on
scatter(MI_list,PIDs(4,:,2)/MI_list);
hold on
scatter(MI_list,PIDs(4,:,3)/MI_list);
legend("MMI", "Idep", "Iccs")
title("Synergy")

%% raw PID atoms results comparison

% MMI
fig=figure('units','normalized','outerposition',[0 0 1 1]);
sgtitle("raw PID atoms")
subplot(2,2,1);
scatter(MI_list,PIDs(1,:,1));
hold on
scatter(MI_list,PIDs(1,:,2));
hold on
scatter(MI_list,PIDs(1,:,3));
legend("MMI", "Idep", "Iccs")
title("Unique X")

subplot(2,2,2);
scatter(MI_list,PIDs(2,:,1));
hold on
scatter(MI_list,PIDs(2,:,2));
hold on
scatter(MI_list,PIDs(2,:,3));
legend("MMI", "Idep", "Iccs")
title("Unique Y")

subplot(2,2,3);
scatter(MI_list,PIDs(3,:,1));
hold on
scatter(MI_list,PIDs(3,:,2));
hold on
scatter(MI_list,PIDs(3,:,3));
legend("MMI", "Idep", "Iccs")
title("Redundancy")

subplot(2,2,4);
scatter(MI_list,PIDs(4,:,1));
hold on
scatter(MI_list,PIDs(4,:,2));
hold on
scatter(MI_list,PIDs(4,:,3));
legend("MMI", "Idep", "Iccs")
title("Synergy")

function [] = save22csv(save, dir, name)

    save = save';
    data = array2table(save);
    data.Properties.VariableNames(1:5) = {'MI','Unique X', 'Unique Y', 'Redundancy', 'Synergy'};
    writetable(data,dir+name+".csv");

end
