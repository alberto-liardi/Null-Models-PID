rng('default') % for reproducibility of random variables
disp("begin run");

% choose task to create the appropriate graph
% task = "single";
task = "average";

if task=="average"
    MI_start = 1e-4;    
    MI_end = 4;
    n_MI = 100;
    n_sample = 1e4;
elseif task=="single"
    MI_start = 1e-1;
    MI_end = 3;
    n_MI = 6;
    n_sample = 1e5;
else, error("Input a correct task"); 
end

% Log-spacing Mutual Information in MI_list
MI_list = logspace(log10(MI_start),log10(MI_end),n_MI);
% MI_list = 0.1:0.2:1;

% vectors containing number of sources S, targets T 
Sources = [4 8]; 
Targets = [2 2];
% Sources = [8];
% Targets = [4];
% Sources = [2]; 
% Targets = [1];

Synergies = NaN(length(MI_list),n_sample);
Redundancies = NaN(length(MI_list),n_sample);
UniqueX = NaN(length(MI_list),n_sample);
UniqueY = NaN(length(MI_list),n_sample);

num_bin = zeros(1,3);
edges = cell(1,3);
x0 = [1 2 0 5 10 20 50];
opts = optimset('display', 'none');
cc = 0;

for k = 1:length(Sources)
    for i = 1:length(MI_list)

        if(mod(i,10)==0), disp(i); end
    
        I = MI_list(i);
        S = Sources(k);
        T = Targets(k);
    
        for j = 1:n_sample
            % sample the matrices
            % NB. for Wishart degrees of freedom must be > size - 1
            
            % source covariance
            Sigma_s = wishrnd(eye(S),S+1);
            % linear coefficients
            A = normrnd(0,1,T,S);
            % target conditional covariance
            Sigma_u = wishrnd(eye(T),T+1);
            
            % find g such that mutual information is MI_tot
            l = 1;
            [alpha, ~, code] = fzero(@(x) fun(x,Sigma_u,A,Sigma_s,I), x0(l), opts);
            
            while( (code == -3 || alpha<0) && l+1 <= length(x0) )
                l = l+1;
                [alpha, ~, code] = fzero(@(x) fun(x,Sigma_u,A,Sigma_s,I), x0(l), opts);
            end
            if( code == -3 || alpha<0 ), cc = cc+1; continue; end 
            
            g = alpha^(-1);

            % check
            % 0.5*log(det(2*pi*exp(1)*(A*Sigma_s*A'+g*Sigma_u))) ...
            % - 0.5*log(det(2*pi*exp(1)*(g*Sigma_u)))
            
            % get the covariance matrix
            Sigma_t = A * Sigma_s * A.' + g*Sigma_u;
            
            % calculate cross covariance terms
            Sigma_cross = A*Sigma_s;

            Sigma_full = [Sigma_s, Sigma_cross';
                          Sigma_cross, Sigma_t];

            I_check = entropy(Sigma_t)+entropy(Sigma_s)-entropy(Sigma_full);
            
%             disp(abs(I-I_check));
%             assert(abs(I-I_check)<1e-6, "problem with MI");
            if(abs(I-I_check)/I>1e-4)
                disp("problem with MI"); disp(abs(I-I_check)); 
                MI_list(i)=NaN;
                continue
            end
            
            % calulate marginal MI (considering 2 sources)
            
            % S1
            Sigma_TS1 = cell(2,2);
            Sigma_TS1{1,1} = Sigma_t;
            Sigma_TS1{1,2} = Sigma_cross(:,1:S/2);
            Sigma_TS1{2,1} = Sigma_TS1{1,2}';
            Sigma_TS1{2,2} = Sigma_s(1:S/2,1:S/2);
            Sigma_TS1 = cell2mat(Sigma_TS1);
            MI_S1 = entropy(Sigma_t)+entropy(Sigma_s(1:S/2,1:S/2)) ...
                    - entropy(Sigma_TS1);
            
            % S2
            Sigma_TS2 = cell(2,2);
            Sigma_TS2{1,1} = Sigma_t;
            Sigma_TS2{1,2} = Sigma_cross(:,S/2+1:end);
            Sigma_TS2{2,1} = Sigma_TS2{1,2}';
            Sigma_TS2{2,2} = Sigma_s(S/2+1:end,S/2+1:end);
            Sigma_TS2 = cell2mat(Sigma_TS2);
            MI_S2 = entropy(Sigma_t)+entropy(Sigma_s(S/2+1:end,S/2+1:end)) ...
                    - entropy(Sigma_TS2);
            
            
            % now calculate PID atoms
            
            % MMI definition
            Red = min(MI_S1, MI_S2);
            Un_X = MI_S1 - Red;
            Un_Y = MI_S2 - Red;
            Syn = I - Un_Y - Un_X - Red;
            
%             fprintf("MI_S1 = " + MI_S1 + "\t MI_S2 = " + MI_S2 + "\n" ...
%                     + "MI = " + I + "\n\n");

            try
                PIDcheck([Un_X, Un_Y, Red, Syn]);
            catch
                cc = cc+1;
                continue
            end
            assert(isreal(Syn), "not real");
            
            % print = "Barrett Redundancy is: " + Red + ...
            %     "\nUnique information of X is: " + Un_X + ...
            %     "\nUnique information of Y is: " + Un_Y + ...
            %     "\nSynergy is: " + Syn + "\n\n";
            % fprintf(print); 
            
            Synergies(i,j) = Syn;
            Redundancies(i,j) = Red;
            UniqueX(i,j) = Un_X;
            UniqueY(i,j) = Un_Y;
        end

        PIDs = {UniqueX(i,:)+UniqueY(i,:); Synergies(i,:); Redundancies(i,:)};
        for pid = 1:3
            [~,edges{pid}] = histcounts(PIDs{pid});
            if length(edges{pid})>num_bin(pid), num_bin(pid)=length(edges{pid}); end
        end
    
    end

    PIDs = {UniqueX+UniqueY, Synergies, Redundancies};

    ID = "S" + sprintf('%01d',S) + "T" + sprintf('%01d',T);
    fprintf("done "+ ID+"\n");

    if task=="single", single_wat(MI_list, num_bin, PIDs, ID);
    elseif task=="average", aver_plot(MI_list, PIDs, ID);
    end

    file_write(ID, n_sample, MI_list, PIDs, task);
     
end

disp("done all!");


%% debugging

plot(linspace(-30,30),arrayfun(@(a) fun(a,Sigma_u,A,Sigma_s,I), linspace(-30,30)))

[alpha, ~, code] = fzero(@(x) fun(x,Sigma_u,A,Sigma_s,I), 5, opts);
g = alpha^-1;
disp(g);
disp(alpha);

%%
% get the covariance matrix
Sigma_t = A * Sigma_s * A.' + g*Sigma_u;

% calculate cross covariance terms
Sigma_cross = A*Sigma_s;

Sigma_full = [Sigma_s, Sigma_cross';
              Sigma_cross, Sigma_t];

% calulate single MI (considering 2 sources)


%% read and compare data

% Sources = [2 4 8 20 2 4 8 20 8 20 20]; 
% Targets = [1 1 1 1 2 2 2 2 4 4 8];
% Sources = [20 20 20 20]; 
% Targets = [1 2 4 8];
Sources = [2 4 8 20];
Targets = [1 1 1 1];
MI_list = logspace(log10(1e-4),log10(4),100);

cmap = colormap(jet);

% IDs = ["S2T1", "S4T1", "S8T1", "S20T1", "S2T2", "S4T2", "S8T2", "S20T2", "S8T4", "S20T4", "S20T8"];
% IDs = ["S20T1", "S20T2", "S20T4", "S20T8"];
IDs = ["S2T1", "S4T1", "S8T1", "S20T1"];
type = ["synergy", "redundancy", "unique"];
for p = 1:3
    pid = cell(1,length(Sources));
    averages = cell(1,length(Sources));
    for i = 1:length(Sources)
        pid{i} = readmatrix("Results_null_model/Gaussian/"+IDs(i)+"/data_"+type(p)+"_av.txt");
        averages{i} = mean(pid{i},2,'omitnan');
    end
    
    fig = figure();
    for i = 1:length(Sources)
        plot(MI_list, averages{i}, 'LineWidth', 1);
        hold on
    end
    title("Null model "+type(p)+" averages comparison",'FontSize',15,'interpreter','latex');
    ylabel(type(p)+" averages");
    xlabel("Mutual Information");
    legend(IDs, 'Location','northwest');
    save = "Results_null_model/Gaussian//All_averages_"+type(p)+".png";
    exportgraphics(fig,save,'Resolution',300);
end

%% PLOT ALL PID ATOMS IN THE SAME GRAPH (only for a system)
Sources = 2;
Targets = 1;
MI_list = logspace(log10(1e-4),log10(4),100);

fig = figure();
cmap = colormap(jet);

type = ["synergy", "redundancy", "unique"];
for p = 1:3
    pid = readmatrix("Results_null_model/Gaussian/S2T1/data_"+type(p)+"_av.txt");
    averages = mean(pid,2,'omitnan');
        
    plot(MI_list, averages'./MI_list, 'LineWidth', 1);
    hold on
end
% plot(MI_list, MI_list, 'LineWidth', 1);


title("NMI PID atoms averages S2T1",'FontSize',15,'interpreter','latex');
ylabel("NMI PID atoms");
xlabel("Mutual Information");
legend([type], 'Location','northwest');
save = "Results_null_model/Gaussian//S2T1_all_PID_NMI.pdf";
exportgraphics(fig,save,'Resolution',300);


function [] = single_wat(MI_list, num_bin, PIDatoms, ID)

    type = ["Unique", "Synergy", "Redundancy"];
    for k = 1:3

        X = zeros(length(MI_list),num_bin(k));
        Y = zeros(length(MI_list),num_bin(k));
        for j = 1:length(MI_list)
            Y(j,:) = MI_list(j);
        end
        Z = zeros(length(MI_list),num_bin(k));
    
        for i = 1:length(MI_list)
            pid = PIDatoms{k};
            
            [NN,edges] = histcounts(pid(i,:), num_bin(k));
            NN = NN./sum(NN);              % normalise w.r.t. N_tot
            centers = movmean(edges,2);
            X(i,:) = centers(2:end);
            Z(i,:) = NN;
        
        end
    
        % waterfall
        fig = figure();
        waterfall(X,Y,Z);
        xlabel(type(k));
        ylabel("Mutual Information");
        zlabel("Counts");
        title("Wishart "+type(k)+" waterfall "+ID,'FontSize',15,'interpreter','latex');
        save = "Results_null_model/Gaussian/Single_waterfall"+ID+"_"+type(k)+".png";
        exportgraphics(fig,save,'Resolution',300);
    end

end

function [] = aver_plot(MI_list, PIDatoms, ID)
    
    type = ["Unique", "Synergy", "Redundancy"];
    for k = 1:3
        % calculate and plot means
        averages = zeros(1,length(MI_list));
        for i = 1:length(MI_list)
            pid = PIDatoms{k};
            averages(i) = mean(pid(i,:),'omitnan');
        end
    
        fig = figure();
    
        plot(MI_list, averages, 'LineWidth',1);
        hold on
        b = polyfit(MI_list,averages,1);
        plot(MI_list, polyval(b,MI_list));
        hold on 
        c = polyfit(MI_list,averages,2);    
        plot(MI_list, polyval(c,MI_list));
    
        title("Wishart "+type(k)+" averages "+ID,'FontSize',15,'interpreter','latex');
        ylabel(type(k)+" averages");
        xlabel("Mutual Information");
        legend("average","linear fit","quadratic fit");
        save = "Results_null_model/Gaussian//Averages_"+ID+"_"+type(k)+".png";
        exportgraphics(fig,save,'Resolution',300);
    end
end

function [] = file_write(ID, n_sample, MI_list, PIDatoms, task)

    if task=="single", tsk = "sin";
    elseif task=="average", tsk = "av";
    end

    path = "Results_null_model/Gaussian/"+ID;
    if not(isfolder(path)), mkdir(path); end

    name = path+"/sum_"+tsk+".txt";
    fileID = fopen(name,'w');
    fprintf(fileID,'%s \t %s \t %s \t %s \t %s \n',...
            'Initial MI', 'Final MI', 'Number of MI', 'Number of samplings', 'ID');
    sum = [MI_list(1), MI_list(end), length(MI_list), n_sample];
    fprintf(fileID,'%d \t %d \t %d \t %d \t %s \n \n', sum, ID);

    type = ["unique", "synergy", "redundancy"];
    for k = 1:3
        name = path+"/data_"+type(k)+"_"+tsk+".txt";
        writematrix(PIDatoms{k},name);
    end

end

function H = entropy(Sigma)
    
    arg = det(2*pi*exp(1)*Sigma);    
    H = 0.5*log(arg);

end

function y = fun(x,Sigma_u,A,Sigma_s,I)

    T = length(Sigma_u);
    y = det(eye(T)+ x * (Sigma_u \ A * Sigma_s * A.'))-exp(2*I);

end