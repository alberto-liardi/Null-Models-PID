%% MAXIMUM REDUNDANCY SYSTEM
% X=Y+small white noise
% Therefore Unique = 0; Redundancy = MIx = MIy, Syn = MI-Red

rng('default') % for reproducibility of random variables
disp("begin run");

% MAXIMUM UNIQUE INFORMATION
% gate = "XOR2";

% MAXIMUM SYNERGY - knowing just one variable does not tell you anything
gate = "XOR";

% setting probabilities
ps = [0.25,0.25,0.25,0.25];
pe = 0.4;

% MAXIMUM REDUNDANCY
% gate = "OR2";
% eps = 0.01;
% ps = [eps/2,(1-eps)/2,(1-eps)/2,eps/2];
% pe = 0.5-eps;

% random system
% gate = "OR";
% ps = [0.1,0.2,0.4,0.3];
% pe = 0.4;

if gate == "XOR"
    boolfun = @XORgate;
elseif gate == "XOR2"
    boolfun = @XOR2gate;
elseif gate == "XOR3"
    boolfun = @XOR3gate;
elseif gate == "OR"
    boolfun = @ORgate;
elseif gate == "OR2"
    boolfun = @OR2gate;
elseif gate == "OR3"
    boolfun = @OR3gate;
elseif gate == "OR4"
    boolfun = @OR4gate;
end 

[pz0,pz1] = boolfun(ps);

pt0 = pz0*(1-pe) + pz1*pe;
pt1 = pz1*(1-pe) + pz0*pe; % 1-pt0;

HT = -pt0*log(pt0)-pt1*log(pt1);
HTcondZ = -pe*log(pe)-(1-pe)*log(1-pe);
MI = HT - HTcondZ;

[MIx, MIy] = DiscreteMarginalMI(ps,pe,HT,gate);

% now calculate PID atoms
Red = min(MIx, MIy);
UnX = MIx - Red;
UnY = MIy - Red;
Syn = MI - UnX - UnY - Red;

UnX, UnY, Red, Syn, return
% constructing and running the null model
null_PIDs = MI_null_model_Discrete(MI, 7e4);
null_Reds = null_PIDs(3,:);
null_Syns = null_PIDs(4,:);
null_Uns = null_PIDs(1,:) + null_PIDs(2,:);

figure();
histogram(null_Reds);
hold on
xline(Red,'--', 'LineWidth', 2, 'Color', 'black');
xlabel('Redundancies');
ylabel('Counts');
title("Redundancy distribution for MI = "+num2str(MI));
legend('Red', "MAX = "+sprintf('%.2f',Red));

figure();
histogram(null_Syns);
hold on
xline(Syn,'--', 'LineWidth', 2, 'Color', 'black');
xlabel('Synergies');
ylabel('Counts');
title("Synergy distribution for MI = "+num2str(MI))
legend('Syn', "MAX = "+sprintf('%.2f',Syn));

figure();
histogram(null_Uns);
hold on
xline(UnX+UnY,'--', 'LineWidth', 2, 'Color', 'black');
xlabel('Unique Information');
ylabel('Counts');
title("Unique distribution for MI = "+num2str(MI))
legend('Un', "MAX = "+sprintf('%.2f',UnX+UnY));

quantPID = comp_quantile(null_PIDs, [UnX,UnY,Red,Syn]');
disp(quantPID);

disp("end of execution");


function [pz0,pz1] = XORgate(ps)
    pz0 = ps(1)+ps(4);
    pz1 = 1-pz0;
end

function [pz0,pz1] = XOR2gate(ps)
    pz0 = ps(1)+ps(2);
    pz1 = 1-pz0;
end

function [pz0,pz1] = XOR3gate(ps)
    pz0 = ps(1)+ps(3);
    pz1 = 1-pz0;
end

function [pz0,pz1] = ORgate(ps)
    pz0 = ps(1);
    pz1 = 1-ps(1);
end

function [pz0,pz1] = OR2gate(ps)
    pz0 = ps(2);
    pz1 = 1-pz0;
end

function [pz0,pz1] = OR3gate(ps)
    pz0 = ps(3);
    pz1 = 1-pz0;
end

function [pz0,pz1] = OR4gate(ps)
    pz0 = ps(4);
    pz1 = 1-pz0;
end

