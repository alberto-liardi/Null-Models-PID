%% function that returns quantile of Phi measures with the null model

% model: struct with fields:
% name - "Gauss" or "VAR"
% S sources - T targets - n Number of runs - p array of model orders

% NB: only insert ONE model order, UnX and other PID atoms must be arrays,
% not cells of arrays.

function [qPhi_wms, qPhi_R, qPhi_red] = QuantPhi(Phi_wms, Phi_R, Phi_red, MI, model)

if isfield(model, 'p')
    % change length in size(,1)
    assert(length(MI)==length(model.p));
else 
    model.p=ones(length(MI));
end
if ~isfield(model, 'n'), model.n=[]; end
if ~isfield(model, 'S'), model.S=2; end
if ~isfield(model, 'T'), model.T=1; end

nPhis = cell(1,length(MI));

if model.name == "VAR"
    for l = 1:length(MI)
        [~, nPhis{l}] = MI_null_model_VAR(MI(l),model.S,model.p(l),model.n,['n','y']);
    end
else 
    error("Not a valid model inserted, exiting."); 
end

qPhis = zeros(3,length(MI));
for m = 1:length(MI) 
    values = [Phi_wms(m), Phi_R(m), Phi_red(m)]';
    qPhis(1:3,m) = comp_quantile(nPhis{m},values);
end

qPhi_wms = qPhis(1,:);  qPhi_R = qPhis(2,:);  qPhi_red = qPhis(3,:);

% TODO/TOTRY: PIDs{1:length(MI)} = MI_null_model_VAR(MI(l),model.S,model.p(l),model.n);
    