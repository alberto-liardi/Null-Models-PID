rng('default') % for reproducibility of random variables
disp("start of the program");

% set parameters for the study 
x0 = [1 -1 2 -2 5 -5];
MI_list = [0.4];
N_runs = 100;
gate = "XOR";

opts = optimset('display', 'none');
cc = 0;

if gate == "XOR"
    boolfun = @XORgate;
elseif gate == "OR"
    boolfun = @ORgate;
end

assert(max(MI_list)<-log(0.5), ...
       "Value of Mutual Information is incompatible with the system!")

Syn = NaN(1,N_runs); Red = NaN(1,N_runs);
UnX = NaN(1,N_runs); UnY = NaN(1,N_runs);

% sample the parameters for the dirichlet from a uniform distribution
% a = [0.1 0.1, 0.1, 0.1];
% a = repmat(5,1,4);
a = repmat(unifrnd(0.1, 3),1,4);

for MI = MI_list
    for j = 1:N_runs
    % sample the probabilities of the Sources S from the dirichlet (only once)
    ps = drchrnd(1,a);
    
    % get the probability of Z=f(X,Y)
    [pz0,pz1] = boolfun(ps);

    % for some parameters of the Dirichlet the MI value may not be reached,
    % check with the maximum of MI (at x=1 or x=0) and skip the step if needed
    if -pz1*log(pz1)-pz0*log(pz0)<MI, disp("skipping step!"), continue; end
    
    % find the discrete noise epsilon to get the value of MI
    l = 1;
    [f, ~, code, out] = fzero(@(x) discrfun(x, pz0, pz1, MI), x0(l), opts);
    
    while( code ~= 1 && l+1 <= length(x0) )
        l = l+1;
        [f, ~, code, out] = fzero(@(x) discrfun(x, pz0, pz1, MI), x0(l), opts);
    end
    % if( code ~= 1 ), cc = cc+1; fprintf("optimisation failed "+sprintf('%01d',code)+"\n"); 
    %                  continue; end 
    
    % now do one more calculation
    pe = Sigmoid(f,0.5);
    % if ~isfinite(discrfun(f, pz0, pz1, MI))
    %     cc = cc+1;
    %     disp("Optimisation function returned NaN value!");
    %     continue;
    % end
    pt0 = pz0*(1-pe) + pz1*pe;
    pt1 = pz1*(1-pe) + pz0*pe; % 1-pt0;
    
    HT = -pt0*log(pt0)-pt1*log(pt1);
    HTcondZ = -pe*log(pe)-(1-pe)*log(1-pe);
    MI_calc = HT - HTcondZ;
    % if(abs(MI_calc-MI)>MI/1e4), cc=cc+1; 
    %     disp("Mismatch between true MI and MI calculated during the optimisation!")
    %     disp(abs(MI_calc-MI)); %disp(MI); 
    %     continue; end
    
    % now need to calculate the marginals MIx, MIy
    [MIx, MIy] = DiscreteMarginalMI(ps,pe,HT,gate)

    Red(j) = min(MIx, MIy);
    UnX(j) = MIx - Red(j);
    UnY(j) = MIy - Red(j);
    Syn(j) = MI - UnX(j) - UnY(j) - Red(j);

%     Red, Syn, UnX, UnY
%     pe

    end
end

disp("end of the program");

function y = discrfun(x, pz0, pz1, MI_value)

%     disp(x); disp("here");
    
    x = Sigmoid(x,0.5);
%     disp(x);
    % can calculate the entropy of the target T=Z+epsilon manually:
    pt0 = pz0*(1-x) + pz1*x;
    pt1 = pz1*(1-x) + pz0*x; % 1-pt0;
    
    HT = -pt0*log(pt0)-pt1*log(pt1);
    HTcondZ = -x*log(x)-(1-x)*log(1-x);
    I = HT - HTcondZ;
    
    y = MI_value - I;

end

function y = Sigmoid(x,M)
    
    y = M/(1+exp(-x));

end

function [pz0,pz1] = XORgate(ps)
    pz0 = ps(1)+ps(4);
    pz1 = ps(2)+ps(3);
end

function [pz0,pz1] = ORgate(ps)
    pz0 = ps(1);
    pz1 = 1-ps(1);
end

function r = drchrnd(n,a)
    p = length(a);
    r = gamrnd(repmat(a,n,1),1,n,p);
    r = r ./ repmat(sum(r,2),1,p);
end

