rng('default') % for reproducibility of random variables
javaaddpath("infodynamics_new.jar");

disp(" "); disp("---- Start of compilation ----"); disp(" ");

subject = 'H';
condition = 'S';
subject_num = 4;

% loading data (if necessary)
if exist('data',     'var') == 0
    disp("Loading data...");

    if subject == 'M'
        % !!! IMPORTANT: Monkey data must be in a folder called "Monkeys_data"
        data = KTMDLoader("Kin2", condition, "Monkeys_data");
    
    elseif subject == 'H'
        % !!! IMPORTANT: Human data must be in a folder called "ValdasPreprocessed"
        [data, fs] = ValdasLoad(subject_num,condition);

    end
    disp("... done.");
end

% setting parameters of the data
tot_epochs = length(data);
tot_channels = size(data{1},2);
tot_period = size(data{1},1);

% setting parameter for the analyses
channels = 20;
epochs = 30; if tot_epochs<epochs, epochs=tot_epochs; end
mmorder = 30;
target = "joint";
morders = [1,mmorder];
tot_N = 2;

% choosing redundancy function
fun = ["MMI", "Idep", "Iccs"];
Red_fun = ['y', 'y', 'y'];
rn = 0;
for r = 1:length(Red_fun)
    if Red_fun(r) == 'y', rn = rn+1; end
end
Red_num = rn;

% saving results
Redundancy = cell(length(fun),length(morders));
Synergy = cell(length(fun),length(morders));
Unique_X = cell(length(fun),length(morders));
Unique_Y = cell(length(fun),length(morders));

path = "Results/S4_c20/Sleep/";
if not(isfolder(path)), mkdir(path); end

k = cell(1,length(morders));
morder_hist = cell(1,length(morders));
CSER = cell(1,length(morders));
tc = cell(1,length(morders));
oIn = cell(1,length(morders));
sIn = cell(1,length(morders));
tcCalc = infodynamics.measures.continuous.gaussian.MultiInfoCalculatorGaussian();
oCalc = infodynamics.measures.continuous.gaussian.OInfoCalculatorGaussian();
sCalc = infodynamics.measures.continuous.gaussian.SInfoCalculatorGaussian();

% start the timer
t0 = tic();

% V_hist = {};

% repeating the simulation tot_N times
for N = 1:1:tot_N

    % choosing #channels different random channels and some random epochs
    c = randperm(tot_channels, channels);
    e = randperm(tot_epochs, epochs);

    % building time series
    Z = zeros(channels, tot_period, epochs);
    for j = 1:channels
        for n = 1:epochs
            Z(j,:,n) = data{e(n)}(:,c(j))';
        end
    end

    % filtering
    if subject == 'M'
        ds = 5; fs = 1000;
    elseif subject == 'H'
        ds = fs/200;
    end
        
    Y = zeros(channels, tot_period/ds, epochs);
    for j = 1:channels
        for n = 1:epochs
            Y(j,:,n) = LowFilter(Z(j,:,n), ds, fs);
        end
    end
    Z = Y; 

    for l = 1:length(morders)
        
        % fitting the evolution coefficient, model order and (optional) residuals
        mord = morders(l);
        [A,p,V,E] = EvolCoefFit(Z, mord, false);

        % computing some quantities (CSER, Oinfo, total correlation)
        morder_hist{l}(end+1) = p;
        CSER{l}(end+1) = .5*log(det(2*pi*exp(1)*V));
        tcCalc.initialise(channels);        tcCalc.setCovariance(V);    
        tc{l}(end+1) = tcCalc.computeAverageLocalOfObservations();
        oCalc.initialise(channels);     oCalc.setCovariance(V);
        oIn{l}(end+1) = oCalc.computeAverageLocalOfObservations();
        sCalc.initialise(channels);     sCalc.setCovariance(V);
        sIn{l}(end+1) = sCalc.computeAverageLocalOfObservations();
%         V_hist{end+1}=V;

        % computing residual correlation (only if 2 channels are considered!)
        if(channels==2)
            corr = diag(diag(V).^(-1/2))*V*diag(diag(V).^(-1/2));
            assert(abs(corr(1,2)-corr(2,1))<1e-7, "correlation not symmetric");
            assert(abs(corr(1,1)-corr(2,2))<1e-7 && abs(corr(1,1)-1)<1e-7, ...
                   "correlation has not ones on diagonal");
            k{l}(end+1) = corr(1,2);
        end

        % computing PID atoms
        PID_atoms= PIDcalc_Mult(p,A,V,channels/2,channels/2,... 
                                target, Red_fun);

        for r = 1:length(Red_fun)
            if (Red_fun(r) == 'y')
                
                % check PID results
                PIDcheck(PID_atoms(r,:));

                % store results
                Unique_X{r,l}(end+1) = PID_atoms(r,1);
                Unique_Y{r,l}(end+1) = PID_atoms(r,2);
                Redundancy{r,l}(end+1) = PID_atoms(r,3);
                Synergy{r,l}(end+1) = PID_atoms(r,4);

            end
        end
    end

end

% stop the timer and print elapsed time 
disp(" "); toc(t0);

% saving everything on file

for r = 1:length(Red_fun)
    if (Red_fun(r) == 'y')
        final_path = path+fun(r)+'/';
        if not(isfolder(final_path)), mkdir(final_path); end

        T = table(Redundancy{r,1}', Redundancy{r,2}', Synergy{r,1}', Synergy{r,2}', ...
                  Unique_X{r,1}', Unique_X{r,2}', Unique_Y{r,1}', Unique_Y{r,2}', ...
                  'VariableNames', ["Redundancy VAR(1)", "Redundancy VAR(p)", ...
                  "Synergy VAR(1)", "Synergy VAR(p)", "Unique information X VAR(1)", ...
                  "Unique information X VAR(p)", "Unique information Y VAR(1)", ...
                  "Unique information Y VAR(p)"] );
        
        writetable(T, final_path+'data.txt','Delimiter','\t');
        
        % save some parameters
        fileID = fopen(final_path + 'summary.txt','w');
        fprintf(fileID,'%s \t %s \t %s \t %s \t %s \t %s \n',...
                'channels', 'epochs', 'max model order', 'N. runs', 'target', 'Redundancy');
        sum = [channels, epochs, mmorder, tot_N];
        fprintf(fileID,'%d \t %d \t %d \t %d \t %s \t %s \n', sum, target, fun(r));
        fclose(fileID);

        if subject == 'H' && condition=='S'
            fileID = fopen(final_path + 'ID.txt','w');
            fprintf(fileID,'%s \t %s \t %s \t %s \t \n','ID', 'Subject', 'Nights', 'Awakening');
            sum = [ID, subject_num, Nights, Awakenings];
            fprintf(fileID,'%d \t %d \t %d \t %d \n', sum);
            fclose(fileID);
        end

    end
end

if channels~=2
    k{1} = zeros(1,tot_N);
    k{2} = zeros(1,tot_N);
end

% save CSER, Total Correlation, O-Information, S-information, Residual Correlation
T = table(CSER{1}', CSER{2}', tc{1}', tc{2}', oIn{1}', oIn{2}', morder_hist{1}', ...
          morder_hist{2}', sIn{1}', sIn{2}', k{1}', k{2}', ...
          'VariableNames', ["CSER Var(1)", "CSER Var(p)", ...
          "Tot Corr Var(1)", "Tot Corr Var(p)", "O info Var(1)", ...
          "O info Var(p)", "Model Order Var(1)", "Model Order Var(p)", ...
          "S info Var(1)", "S info Var(p)", "Residual Corr Var(1)", "Residual Corr Var(p)"] );

writetable(T, path+'InfoDyn.txt','Delimiter','\t');

tc_n = cell(1,length(morders));
oIn_n = cell(1,length(morders));
tc_n{1} = tc{1}./sIn{1}; tc_n{2} = tc{2}./sIn{2};
oIn_n{1} = oIn{1}./sIn{1}; oIn_n{2} = oIn{2}./sIn{2};

T = table(tc_n{1}', tc_n{2}', oIn_n{1}', oIn_n{2}', ...
          'VariableNames', ["Tot Corr Norm Var(1)", "Tot Corr Norm Var(p)", ...
          "O info Norm Var(1)", "O info Norm Var(p)"] );

writetable(T, path+'InfoDyn_Norm.txt','Delimiter','\t');


disp(" "); disp("---- End of compilation ----"); disp(" ");