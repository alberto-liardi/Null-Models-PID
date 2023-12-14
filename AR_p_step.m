rng('default') % for reproducibility of random variables

disp(" "); disp("---- Start of compilation ----"); disp(" ");

p=1;

% using newly-made function:
A = RandEvolMatrix(p);

% first need a set of all evolution matrices
% A1 = [0.01 0.0; 0.10 0.10];
% A2 = [0.1 0.; 0.1 0.1];
% A3 = [0. 0.; 0. 0.]+0.1;
% A = { A1, A2, A3 };


% define vectors (need to be array of size at least p)
X = 1:1:3;
Y = 10:-1:8;

% using a matrix which contains both time series
Z = vertcat( X, Y ); 

% evolving the system for some steps following the above AR
Z = SystEvol(Z,5,p,A);

% now creating the matrices for the Lyapunov equation

% coefficient matrix
A_L = horzcat( eye(2*(p-1)), zeros(2*(p-1),2) ); 
A_L = vertcat( cell2mat(A) , A_L);
disp("Coefficient evolution matrix: "); disp(A_L);

% covariance of random noise variables
Lambda_L = horzcat( eye(2), zeros(2, 2*(p-1)) ); 
Lambda_L = vertcat( Lambda_L , zeros(2*(p-1), 2*p) );    
disp("Noise covariance matrix: "); disp(Lambda_L);

% solving Lyapunov equation for Gamma: 
Gamma = dlyap(A_L,Lambda_L);
disp("Gamma matrices: "); disp(Gamma);

% deconstructing and reordering gamma matrices
split = zeros(1,p)+2;
Gamma = Gamma([1 2],:); % only the first line is linear independent
disp(Gamma);
Gamma = mat2cell(Gamma, 2, split);

% now find last gamma matrix (Gamma_p) from YW equation:
Gamma_new = 0;
for l = 1:1:p
    Gamma_new = Gamma_new + A{l}*Gamma{p-l+1};
end
Gamma{end+1} = Gamma_new; 


% can now construct Sigma matrices

% covariance matrices
sigma_cov_XX = zeros(p,p);
sigma_cov_XY = zeros(p,p);
sigma_cov_YY = zeros(p,p);
for s = 1:1:p
    for t = s:1:p
        sigma_cov_XX(s,t) = Gamma{t-s+1}(1,1);
        sigma_cov_XY(s,t) = Gamma{t-s+1}(1,2);
        sigma_cov_YY(s,t) = Gamma{t-s+1}(2,2);
    end
end
% can be written more efficiently just but tayloring the first row
% adequately

% concatenating matrices to obtain the full covariance matrix (_tot)
sigma_cov_XX = sigma_cov_XX' + triu(sigma_cov_XX,1);
sigma_cov_YY = sigma_cov_YY' + triu(sigma_cov_YY,1);
sigma_cov_XY = sigma_cov_XY' + triu(sigma_cov_XY,1);
sigma_cov_YX = sigma_cov_XY';

% now constructing the total covariance matrix with the proper ordering
sigma_cov = cell(p,p);
sigma_cov(:) = {[0 0; 0 0]};
for s = 1:1:p
    for t = s:1:p
        sigma_cov{s,t} = Gamma{t-s+1};
    end
end
sigma_cov = cell2mat(sigma_cov);
sigma_cov = sigma_cov - tril(sigma_cov,-1);
sigma_cov_tot = triu(sigma_cov,1) + sigma_cov';

%cross covariance matrices - take second row of all Gamma_k matrices (k!=0)
sigma_cross_YX = zeros(1,p);
sigma_cross_YY = zeros(1,p);
for s = 1:1:p
    sigma_cross_YX(s) = Gamma{s+1}(2,1);
    sigma_cross_YY(s) = Gamma{s+1}(2,2);
end

% total cross-covariance matrix with the proper ordering
sigma_cross_tot = horzcat( Gamma{2:end} );
sigma_cross_tot = sigma_cross_tot(2,:);

% variance of the target
varY = Gamma{1}(2,2); %==Gamma_0(2,2)
% disp(varY);

% computing partial variances
partial_var_tot = varY - sigma_cross_tot / sigma_cov_tot * sigma_cross_tot';
partial_var_YY = varY - sigma_cross_YY / sigma_cov_YY * sigma_cross_YY';
partial_var_YX = varY - sigma_cross_YX / sigma_cov_XX * sigma_cross_YX';

print= "partial variance sigma(Y;X^-,Y^-): " + partial_var_tot + "\n" + ...
       "partial variance sigma(Y,X^-): " + partial_var_YX + "\n" + ... 
       "partial variance sigma(Y,Y^-): " + partial_var_YY + " \n";
fprintf(print);

% can now compute physical quantities:

% predictive information PE (==mutual information MI)
PE = 0.5*log( varY / partial_var_tot);
% self entropy SE
SE = 0.5*log( varY / partial_var_YY);
% trasfer entropy TE
TE = 0.5*log( partial_var_YY / partial_var_tot);
% cross entropy CE
CE = 0.5*log( varY / partial_var_YX);
% conditional self entropy cSE
cSE = 0.5*log( partial_var_YX / partial_var_tot);

print= "\n" + "predictive information: " + PE + " \n" + ...
        "self entropy: " + SE + " \n" + "transfer entropy: " + TE + ...
        " \n" + "cross entropy: " + CE + " \n" + ...
        "conditional self entropy: " + cSE + " \n" ; 
fprintf(print);

% ------------------------------------------------------------------ %
% now want to compute physical quantities with definition of entropy %
% ------------------------------------------------------------------ %

fprintf("\nNow calculating mutual information with entropies.");

sigma_cov_TS = cell(p+1,p+1);
sigma_cov_TS(:) = {[0 0; 0 0]};
for s = 1:1:p+1
    for t = s:1:p+1
        sigma_cov_TS{s,t} = Gamma{t-s+1};
    end
end
sigma_cov_TS = cell2mat(sigma_cov_TS);
sigma_cov_TS = sigma_cov_TS - tril(sigma_cov_TS,-1);
sigma_cov_tot_TS= triu(sigma_cov_TS,1) + sigma_cov_TS';

sigma_cov_S = cell(p,p);
sigma_cov_S(:) = {[0 0; 0 0]};
for s = 1:1:p
    for t = s:1:p
        sigma_cov_S{s,t} = Gamma{t-s+1};
    end
end
sigma_cov_S = cell2mat(sigma_cov_S);
sigma_cov_S = sigma_cov_S - tril(sigma_cov_S,-1);
sigma_cov_tot_S = triu(sigma_cov_S,1) + sigma_cov_S';

% now have to add first row and column
row = horzcat( Gamma{2:end} );
row = row(1,:);

column = vertcat( Gamma{1} );
column = vertcat( column(1,1), row' );

sigma_cov_tot_S = vertcat( row, sigma_cov_tot_S );
sigma_cov_tot_S = horzcat( column, sigma_cov_tot_S );

% variance of the target
sigma_cov_tot_T = Gamma{1}(2,2); %==Gamma_0(2,2)

% calculating entropies
H_T = .5*log(det(2*pi*exp(1) * sigma_cov_tot_T));
H_S = .5*log(det(2*pi*exp(1) * sigma_cov_tot_S));
H_TS = .5*log(det(2*pi*exp(1) * sigma_cov_tot_TS));

% mutual information
MI = H_T + H_S - H_TS;
fprintf("\nMutual Information computed this way is:"); disp(MI);



% Now to do PID some matrices that have been calculated before are needed
% again (e.g. XX, YY, ...)


% covariance matrices
sigma_X = zeros(p+1,p+1);
sigma_Y = zeros(p+1,p+1);
for s = 1:1:p+1
    for t = s:1:p+1
        sigma_X(s,t) = Gamma{t-s+1}(1,1);
        sigma_Y(s,t) = Gamma{t-s+1}(2,2);
    end
end
% can be written more efficiently just but tayloring the first row
% adequately

% concatenating matrices to obtain the full covariance matrix (_tot)
sigma_X_tot = sigma_X' + triu(sigma_X,1);
sigma_Y_tot = sigma_Y' + triu(sigma_Y,1);

row = [];
for l = 1:1:p+1
    row(end+1) = Gamma{l}(2,1);
end
sigma_YX_tot = vertcat( row, sigma_X_tot );
column = vertcat( Gamma{1}(2,2), row' );
sigma_YX_tot = horzcat( column, sigma_YX_tot);

sigma_Y_past = sigma_Y_tot(2:end,2:end); 

% calculating I(Y,X_tot):
H_T = .5*log(det(2*pi*exp(1) * sigma_cov_tot_T));
H_S = .5*log(det(2*pi*exp(1) * sigma_X_tot));
H_TS = .5*log(det(2*pi*exp(1) * sigma_YX_tot));

MI_YX = H_T + H_S - H_TS;
fprintf("Mutual Information brought by X_tot is: " + MI_YX);


% calculating I(Y,Y_past):
H_T = .5*log(det(2*pi*exp(1) * sigma_cov_tot_T));
H_S = .5*log(det(2*pi*exp(1) * sigma_Y_past));
H_TS = .5*log(det(2*pi*exp(1) * sigma_Y_tot));

MI_YY = H_T + H_S - H_TS;

fprintf("\nMutual Information brought by Y_past is:"); disp(MI_YY);

% problems with imaginary numbers and with positivity of the MI obtained

% now calculating PID quantities
Red = min(abs(MI_YX), abs(MI_YY)); % <- NB I put the absolute value 
Uni_X = MI_YX - Red;
Uni_Y = MI_YY - Red;
Syn = MI - Uni_Y - Uni_X - Red;

print = "Barrett Redundancy is: " + Red + ...
        "\nUnique information of X is: " + Uni_X + ...
        "\nUnique information of Y is: " + Uni_Y + ...
        "\nSynergy is: " + Syn + "\n";
fprintf(print); 

disp(" "); disp("---- End of compilation ----"); disp(" ");