function Information = TimeSeriesAnalyses(data, subject_info, par)

javaaddpath("infodynamics.jar");

% setting parameters of the data
tot_epochs = length(data);
tot_channels = size(data{1},2);
tot_period = size(data{1},1);

% setting parameter for the analyses
if ~isfield(par, 'channels'), par.channels = 2; end
if ~isfield(par, 'epochs'), par.epochs = 30; end, if tot_epochs<par.epochs, par.epochs=tot_epochs; end
if ~isfield(par, 'mmorder'), par.mmorder = 30; end
target = "joint";
morders = [1,par.mmorder];
if ~isfield(par, 'runs'), par.runs = 100; end

% choosing redundancy function
fun = ["MMI", "Idep", "Iccs"];
if ~isfield(par, 'Red_fun'), par.Red_fun = ['y', 'y', 'y']; end

% storing results
Redundancy = cell(length(fun),length(morders));
Synergy = cell(length(fun),length(morders));
Unique_X = cell(length(fun),length(morders));
Unique_Y = cell(length(fun),length(morders));

k = cell(1,length(morders));
morder_hist = cell(1,length(morders));
CSER = cell(1,length(morders));
tc = cell(1,length(morders));
oIn = cell(1,length(morders));
sIn = cell(1,length(morders));
tcCalc = infodynamics.measures.continuous.gaussian.MultiInfoCalculatorGaussian();
oCalc = infodynamics.measures.continuous.gaussian.OInfoCalculatorGaussian();
sCalc = infodynamics.measures.continuous.gaussian.SInfoCalculatorGaussian();

V_hist = cell(1,par.runs);
A_hist = cell(1,par.runs);

phi_WMS = zeros(1,par.runs);
phi_R = zeros(1,par.runs);

% repeating the simulation tot_N times
for N = 1:1:par.runs

    % choosing #channels different random channels and some random epochs
    c = randperm(tot_channels, par.channels);
    e = randperm(tot_epochs, par.epochs);

    % building time series
    Z = zeros(par.channels, tot_period, par.epochs);
    for j = 1:par.channels
        for n = 1:par.epochs
            Z(j,:,n) = data{e(n)}(:,c(j))';
        end
    end

    % filtering
    down_s = subject_info.fs/200;
    Y = zeros(par.channels, tot_period/down_s, par.epochs);
    for j = 1:par.channels
        for n = 1:par.epochs
            Y(j,:,n) = LowFilter(Z(j,:,n), down_s, subject_info.fs);
        end
% TODO - avoid one for 1-cycle... need to check with butterworth filter
%Y(j,:,1:par.epochs) = LowFilter(Z(j,:,1:par.epochs), down_s, subject_info.fs);
    end
    Z = Y;

    for l = 1:length(morders)
        
        % fitting the evolution coefficient, model order and (optional) residuals
        mord = morders(l);
        [A,p,V,~] = EvolCoefFit(Z, mord, false);
                
        % computing some quantities (CSER, Oinfo, Sinfo, Total Correlation)
        morder_hist{l}(end+1) = p;
        CSER{l}(end+1) = .5*log(det(2*pi*exp(1)*V));
        tcCalc.initialise(par.channels);        tcCalc.setCovariance(V);    
        tc{l}(end+1) = tcCalc.computeAverageLocalOfObservations();
        oCalc.initialise(par.channels);     oCalc.setCovariance(V);
        oIn{l}(end+1) = oCalc.computeAverageLocalOfObservations();
        sCalc.initialise(par.channels);     sCalc.setCovariance(V);
        sIn{l}(end+1) = sCalc.computeAverageLocalOfObservations();
        
        if mord == 1,   V_hist{N}=V;    A_hist{N}=A{1};    end

        % computing residual correlation (only if 2 channels are considered!)
        if(par.channels==2)
            corr = diag(diag(V).^(-1/2))*V*diag(diag(V).^(-1/2));
            assert(abs(corr(1,2)-corr(2,1))<1e-7, "correlation not symmetric");
            assert(abs(corr(1,1)-corr(2,2))<1e-7 && abs(corr(1,1)-1)<1e-7, ...
                   "correlation has not ones on diagonal");
            k{l}(end+1) = corr(1,2);
        end

        % computing PID atoms
        [PID_atoms, Gammas] = PIDcalc_Mult(p,A,V,par.channels/2,par.channels/2,... 
                                target, par.Red_fun);

        if(mord==1 && l==1)
            [phi_WMS(N), phi_R(N)] = integrated_info_VAR(sum(PID_atoms(1,:)),Gammas);
        end

        for r = 1:length(par.Red_fun)
            if (par.Red_fun(r) == 'y')
                
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

% saving total correlation of residuals (only for 2 channels)
if par.channels~=2
    k{1} = zeros(1,par.runs);
    k{2} = zeros(1,par.runs);
end

% Normalise information dynamics with respect to the Sigma-Information
tc_n = cell(1,length(morders));
oIn_n = cell(1,length(morders));
tc_n{1} = tc{1}./sIn{1}; tc_n{2} = tc{2}./sIn{2};
oIn_n{1} = oIn{1}./sIn{1}; oIn_n{2} = oIn{2}./sIn{2};


% save everything we calculated in the Information struct
Information.subject_info = subject_info;

Information.parameters = par;
Information.parameters.fun = fun;

Information.UnX = Unique_X;
Information.UnY = Unique_Y;
Information.Red = Redundancy;
Information.Syn = Synergy;

Information.CSER = CSER;
Information.TC = tc;
Information.Sinfo = sIn;
Information.Oinfo = oIn;
Information.mord_hist = morder_hist;
Information.corr = k;

Information.phi_WMS = phi_WMS;
Information.phi_R = phi_R;

Information.TC_norm = tc_n;
Information.Oinfo_norm = oIn_n;

Information.V = V_hist;
Information.A = A_hist;

end