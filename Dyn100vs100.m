% compare 100 vs 1000 results Kin2

path = "Results/Kin2_c20/Awake/InfoDyn.txt";

info = readmatrix(path);
O_100 = info(:,5:6);
tc_100 = info(:,3:4);

O_100_n = O_100 ./ (2*tc_100-O_100);
O_100_n_e_W = 2*sqrt(var(O_100_n))/10;
O_100_norm_W = mean(O_100_n);


path = "Results_old/Kin2_c20/Awake/InfoDyn.txt";

info = readmatrix(path);
O_1000 = info(:,5:6);
tc_1000 = info(:,3:4);

O_1000_n = O_1000 ./ (2*tc_1000-O_1000);
O_1000_n_e_W = 2*sqrt(var(O_1000_n))/sqrt(1000);
O_1000_norm_W = mean(O_1000_n);

% W = (O_100_norm_W - O_1000_norm_W)./(O_100_norm_W + O_1000_norm_W) * 200;


% now SEDATED condition


path = "Results/Kin2_c20/Sleep/InfoDyn.txt";

info = readmatrix(path);
O_100 = info(:,5:6);
tc_100 = info(:,3:4);

O_100_n = O_100 ./ (2*tc_100-O_100);
O_100_n_e_S = 2*sqrt(var(O_100_n))/sqrt(100);
O_100_norm_S = mean(O_100_n);


path = "Results_old/Kin2_c20/Sedated/InfoDyn.txt";

info = readmatrix(path);
O_1000 = info(:,5:6);
tc_1000 = info(:,3:4);

O_1000_n = O_1000 ./ (2*tc_1000-O_1000);
O_1000_n_e_S = 2*sqrt(var(O_1000_n))/sqrt(1000);
O_1000_norm_S = mean(O_1000_n);

% S = (O_100_norm_S - O_1000_norm_S)./(O_100_norm_S + O_1000_norm_S) * 200;

%
%
D100 = O_100_norm_W - O_100_norm_S;
D1000 = O_1000_norm_W - O_1000_norm_S;

% [h,p] = ttest2(D100,D1000,'Vartype','unequal');
[h,p] = ttest(D100-D1000);

disp("done");