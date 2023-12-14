%% MAXIMUM REDUNDANCY SYSTEM
% X=Y+small white noise
% Therefore Unique = 0; Redundancy = MIx = MIy, Syn = MI-Red

rng('default') % for reproducibility of random variables
disp("begin run");

max_mord = 1;

x = 1:100;
X = 0.5*x + normrnd(0,1,1,100);
Y = X+normrnd(0,0.001,1,100);
Z = [X;Y];

% fitness procedure
[A,p,V,~] = EvolCoefFit(Z, max_mord, false);

PID= PIDcalc_Mult(p,A,V,1,1,"joint",['y','n','n']);
Red = PID(1,3);   MI = sum(PID(1,:));

null_PIDs = MI_null_model_VAR(MI);
null_Reds = null_PIDs(3,:);

figure();
histogram(null_Reds);
hold on
xline(Red,'--', 'LineWidth', 2, 'Color', 'black');

xlabel('Redundancies');
ylabel('Counts');
title("Redundancy distribution for MI = "+num2str(MI));
legend('Red', "MAX = "+sprintf('%.2f',Red));

disp("end of execution");


%% MAXIMUM UNIQUE SYSTEM

rng('default') % for reproducibility of random variables
disp("begin run");

max_mord = 1;

x = 1:100;
X = 0.5*x + normrnd(0,1,1,100);
Y = normrnd(0,1,1,100);

Z = [X;Y];

% fitness procedure
[A,p,V,~] = EvolCoefFit(Z, max_mord, false);

PID= PIDcalc_Mult(p,A,V,1,1,"joint",['y','n','n']);
Un = max(PID(1,1), PID(1,2));   MI = sum(PID(1,:));

null_PIDs = MI_null_model_VAR(MI);
null_Uns = null_PIDs(1,:)+null_PIDs(2,:);

figure();
histogram(null_Uns);
hold on
xline(Un,'--', 'LineWidth', 2, 'Color', 'black');

xlabel('Unique');
ylabel('Counts');
title("Unique distribution for MI = "+num2str(MI));
legend('Un', "MAX = "+sprintf('%.2f',Un));

disp("end of execution");

%% MAXIMAL SYNERGISTIC SYSTEM

rng('default') % for reproducibility of random variables
disp("begin run");

max_mord = 1;

a = 0.5+rand(1)*20;
A = [-a a+0.1; -a a+0.1];

V = eye(2);
PID= PIDcalc_Mult(p,{A},V,1,1,"joint",['y','n','n']);
Syn = PID(1,4);   MI = sum(PID(1,:));

null_PIDs = MI_null_model_VAR(MI);
null_Syns = null_PIDs(4,:);

figure();
histogram(null_Syns);
hold on
xline(Syn,'--', 'LineWidth', 2, 'Color', 'black');

xlabel('Synergy');
ylabel('Counts');
title("Synergy distribution for MI = "+num2str(MI));
legend('Syn', "MAX = "+sprintf('%.2f',Syn));

disp("end of execution");
