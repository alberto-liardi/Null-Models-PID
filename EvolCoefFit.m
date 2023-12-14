%% ------ CODE ADAPTED FROM MVGC2/demo/mvgc_demo_var.m ------

%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ AA,morder,V,E ] = EvolCoefFit(X, momax, plotm)

    % VAR model parameters
    
    if ~exist('moact',     'var'), moact     = 6;       end % model order
    
    % VAR model order estimation
    
    if ~exist('moregmode', 'var'), moregmode = 'LWR';   end % VAR model estimation regression mode ('OLS' or 'LWR')
    if ~exist('mosel',     'var'), mosel     = 'LRT';   end % model order selection ('ACT', 'AIC', 'BIC', 'HQC', 'LRT', or supplied numerical value)
    if ~exist('momax',     'var'), momax     = 15; end % maximum model order for model order selection
    
    % VAR model parameter estimation
    
    if ~exist('regmode',   'var'), regmode   = 'LWR';   end % VAR model estimation regression mode ('OLS' or 'LWR')
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~exist('seed',      'var'), seed      = 0;       end % random seed (0 for unseeded)
    if (~exist('plotm',     'var')  || ~plotm),      plotm = [];
    else, plotm = 0; end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%% VAR modelling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Remove temporal mean and normalise by temporal variance.
    % Not strictly necessary, but may help numerical stability
    % if data has very large or very small values.
    
    X = demean(X,true);
    
    % Model order estimation
    
    % Calculate and plot VAR model order estimation criteria up to specified maximum model order.
    
    ptic('\n*** tsdata_to_varmo... ');
    if isnumeric(plotm), plotm = plotm+1; end
    [moaic,mobic,mohqc,molrt] = tsdata_to_varmo(X,momax,moregmode,[],[],plotm);
    ptoc;
    
    % Select and report VAR model order.
    
    morder = moselect(sprintf('VAR model order selection (max = %d)',momax), mosel,'ACT',moact,'AIC',moaic,'BIC',mobic,'HQC',mohqc,'LRT',molrt);
    assert(morder > 0,'selected zero model order! GCs will all be zero!');
    if morder >= momax, fprintf(2,'*** WARNING: selected maximum model order (may have been set too low)\n'); end
    
    % VAR model estimation
    
    % Estimate VAR model of selected order from data.
    
    ptic('\n*** tsdata_to_var... ');
    [A,V,E] = tsdata_to_var(X,morder,regmode);
    ptoc;
    
    % Check for failed regression
    
    assert(~isbad(A),'VAR estimation failed - bailing out');
    
    % Report information on the estimated VAR, and check for errors.
    %
    % IMPORTANT: We check the VAR model for stability and symmetric positive-
    % definite residuals covariance matrix. THIS CHECK IS ESSENTIAL - subsequent
    % routines may fail if there are errors here. If there are problems with the
    % data (e.g. non-stationarity, colinearity, etc.) there's also a good chance
    % they'll show up at this point - and the diagnostics may supply useful
    % information as to what went wrong.
    
    info = var_info(A,V);
    assert(~info.error,'VAR error(s) found - bailing out');
    
    % saving coefficient matrices in cells
    AA = cell(1,morder);
    for l = 1:1:morder
        AA{l} = A(:,:,l);
    end

end