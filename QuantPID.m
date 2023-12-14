%% function that returns quantile of a set of PID atoms with the null model

% model: struct with fields:
% name - "Gauss" or "VAR"
% S sources - T targets - n Number of runs - p array of model orders

% NB: only insert ONE model order, UnX and other PID atoms must be arrays,
% not cells of arrays.

function [qUnX, qUnY, qRed, qSyn] = QuantPID(UnX, UnY, Red, Syn, MI, model)

if isfield(model, 'p')
    % change length in size(,1)
    assert(length(MI)==length(model.p));
else 
    model.p=ones(length(MI));
end
if ~isfield(model, 'n'), model.n=[]; end
if ~isfield(model, 'S'), model.S=2; end
if ~isfield(model, 'T'), model.T=1; end
if ~isfield(model, 'red_fun'), model.red_fun="MMI"; end

if model.red_fun=="MMI", r=1;
elseif model.red_fun=="Idep", r=2;
elseif model.red_fun=="Iccs", r=3;
end

nPIDs = cell(1,length(MI));
% TODO -- keep using nPIDs{l} but then change the for loop down there so as
% to work with quantiles in 3D vectors

if model.name == "Gauss"
    for l = 1:length(MI)
        all_PIDs = MI_null_model_Gauss(MI(l),model.S,model.T,model.n);   
        nPIDs{l} = all_PIDs(:,:,r);
    end

elseif model.name == "VAR"
    for l = 1:length(MI) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         t0 = tic();
        all_PIDs = MI_null_model_VAR(MI(l), model.S, model.p(l), model.n, ...
                                     [model.red_fun], ['y','n']);
        nPIDs{l} = all_PIDs(:,:,r);
%         fprintf("%d\t%d\t", l, model.p(l)); toc(t0); 
%         if mod(l,10)==0, fprintf("doing MI number %d\n",l); end 
    end
else 
    error("Not a valid model inserted, exiting."); 
end

qPIDs = zeros(4,length(MI));
for m = 1:length(MI) 
    values = [UnX(m), UnY(m), Red(m), Syn(m)]';
    qPIDs(1:4,m) = comp_quantile(nPIDs{m},values);
end

qUnX = qPIDs(1,:);  qUnY = qPIDs(2,:);
qRed = qPIDs(3,:);  qSyn = qPIDs(4,:);

% TODO/TOTRY: PIDs{1:length(MI)} = MI_null_model_VAR(MI(l),model.S,model.p(l),model.n);
    