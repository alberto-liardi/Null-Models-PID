
function [ PID, Gamma ] = PIDcalc_Mult(p,A,V,L1,L2,target,Red)

% ---------------------------------------------------------------------- %
% Compute PID quantities with the relation between Entropy - MI  %
% ---------------------------------------------------------------------- %

% MULTIVARIATE GENERALISATION

% L1 length of the first source
% L2 length of the second source
L = L1 + L2; % length of the target

if ~exist('V','var'), V = eye(L); end
assert(iscell(A), "Coefficient matrices must be in a 1xp cell");
assert(0<p, "Enter a valid model order"); 
assert(length(A)==p, "Model order and number of evolution matrices are not the same");
assert(size(A{1},1) == L, "Sources and evolution matrix dimensions are not compatible");

% creating the matrices for the Lyapunov equation

% coefficient matrix
% A_L = horzcat( eye(L*(p-1)), zeros(L*(p-1),L) ); 
% A_L = vertcat( cell2mat(A) , A_L);
% % disp("Coefficient evolution matrix: "); disp(A_L);
% 
% % covariance of random noise variables matrix
% Lambda_L = horzcat( V, zeros(L, L*(p-1)) ); 
% Lambda_L = vertcat( Lambda_L , zeros(L*(p-1), L*p) );    
% % disp("Noise covariance matrix: "); disp(Lambda_L);
% 
% % solving Lyapunov equation for Gamma: 
% Gamma = dlyap(A_L,Lambda_L);
% 
% % deconstructing and reordering gamma matrices
% split = zeros(1,p)+L;
% Gamma = Gamma(1:L,:); % only the first line is linear independent
% Gamma = mat2cell(Gamma, L, split);
% 
% 
% % now find last gamma matrix (Gamma_p) from YW equation:
% Gamma_new = 0;
% for l = 1:1:p
%     Gamma_new = Gamma_new + A{l}*Gamma{p-l+1};
% end
% Gamma{end+1} = Gamma_new; 
% % disp("Gamma matrices: "); 
% % for g = 1:1:p+1
% %     print="Gamma"+g+"\n"; fprintf(print);
% %     disp(Gamma{g});
% % end

B = zeros(L,L,p);
for s = 1:p
    B(:,:,s) = A{s};
end
% disp("START lyap");
G = var_to_autocov(B,V,p);
Gamma = cell(1,p+1);
for s = 1:p+1
    Gamma{s} = G(:,:,s);
end
% disp("END lyap");


% can now construct Sigma matrices

% ----- single future generalisation missing -----

if( target == "joint")

    % T = (T1, T2) is the target
    % S = (S1, S2) is the source
    % TS all of them (target(s) + sources)

    sigma_cov_TS = cell(p+1,p+1);
    sigma_cov_TS(:) = {zeros(L,L)};
    for s = 1:1:p+1
        for t = s:1:p+1
            sigma_cov_TS{s,t} = Gamma{t-s+1};
        end
    end
    % going from the upper diagonal to the full matrix 
    sigma_cov_TS = cell2mat(sigma_cov_TS);
    sigma_cov_TS = sigma_cov_TS - tril(sigma_cov_TS,-1);
    sigma_cov_tot_TS= triu(sigma_cov_TS,1) + sigma_cov_TS';
    
    sigma_cov_tot_S = sigma_cov_tot_TS(1:end-L, 1:end-L);
%     disp(sigma_cov_tot_S);
    
    sigma_cov_tot_T = Gamma{1};
%     disp(sigma_cov_tot_T);
    
    % calculating entropies [log(det(A)) = Tr(log(A))]
%     H_T = .5*(size(sigma_cov_tot_T,1)*log(2*pi*exp(1)) + ...
%               trace(logm(sigma_cov_tot_T)));
%     H_S = .5*(size(sigma_cov_tot_S,1)*log(2*pi*exp(1)) + ...
%               trace(logm(sigma_cov_tot_S)));  
%     H_TS = .5*(size(sigma_cov_tot_TS,1)*log(2*pi*exp(1)) + ...
%               trace(logm(sigma_cov_tot_TS)));

    H_T = 0.5*logdet(sigma_cov_tot_T);
    H_S = 0.5*logdet(sigma_cov_tot_S);
    H_TS = 0.5*logdet(sigma_cov_tot_TS);
    
    % mutual information
    MI = (H_T + H_S - H_TS) / log(2);
%     fprintf("\nMutual Information computed using entropies is:"); disp(MI);
    
    
    % Now to calculate single Mutual Information brought by S1, S2, need to
    % calculate some more covariance matrices
    
    % covariance matrices
    sigma_T1S1 = zeros(L1*(p+1),L1*(p+1));
    sigma_T2S2 = zeros(L2*(p+1),L2*(p+1));
    for s = 1:1:p+1
        for t = s:1:p+1
            sigma_T1S1(s*L1-L1+1:s*L1,t*L1-L1+1:t*L1) = Gamma{t-s+1}(1:L1,1:L1);
            sigma_T2S2(s*L2-L2+1:s*L2,t*L2-L2+1:t*L2) = Gamma{t-s+1}(L1+1:end,L1+1:end);
        end
    end

    % concatenating matrices to obtain the full covariance matrix (_tot)
    sigma_T1S1_tot = sigma_T1S1' + triu(sigma_T1S1,1) - tril(sigma_T1S1,-1)';
    sigma_T2S2_tot = sigma_T2S2' + triu(sigma_T2S2,1) - tril(sigma_T2S2,-1)';

%     disp(sigma_T1S1_tot);
%     disp(sigma_T2S2_tot);

    row_T2 = zeros(L2,L1*(p+1));
    row_T1 = zeros(L1,L2*(p+1));
    for l = 1:1:p+1
        row_T2(1:L2,l*L1-L1+1:l*L1) = Gamma{l}(L1+1:end,1:L1);
        row_T1(1:L1,l*L2-L2+1:l*L2) = Gamma{l}(1:L1,L1+1:end);
    end

    addT2 = Gamma{1}(L1+1:end,L1+1:end);
    addT1 = Gamma{1}(1:L1,1:L1);
    for l = 1:1:p+1
        row_T2(1:L2,l*L1-L1+1:l*L1) = Gamma{l}(L1+1:L1+L2,1:L1);
        row_T1(1:L1,l*L2-L2+1:l*L2) = Gamma{l}(1:L2,L2+1:L1+L2);
    end
    
    sigma_TS1 = vertcat( row_T2, sigma_T1S1_tot );
    column = vertcat( addT2, row_T2' );
    sigma_TS1 = horzcat( column, sigma_TS1);

    sigma_TS2 = vertcat( row_T1, sigma_T2S2_tot );
    column = vertcat( addT1, row_T1' );
    sigma_TS2 = horzcat( column, sigma_TS2);

%     disp(sigma_TS2);
%     disp(sigma_TS1);
    
    sigma_S1 = sigma_T1S1_tot(L1+1:end,L1+1:end);
    sigma_S2 = sigma_T2S2_tot(L2+1:end,L2+1:end);

%     disp(sigma_S1);
%     disp(sigma_S2);

    % calculating I(T,S1):
    H_T = 0.5*logdet(sigma_cov_tot_T);
    H_S = 0.5*logdet(sigma_S1);
    H_TS = 0.5*logdet(sigma_TS1);
    
    MI_TX = (H_T + H_S - H_TS) / log(2);
%     fprintf("Mutual Information brought by Source 1 is: " + MI_TX);


    % calculating I(T,S2):
    H_T = 0.5*logdet(sigma_cov_tot_T);
    H_S = 0.5*logdet(sigma_S2);
    H_TS = 0.5*logdet(sigma_TS2);
    
    MI_TY = (H_T + H_S - H_TS) / log(2);
    
%     fprintf("\nMutual Information brought by Source 2 is:"); disp(MI_TY);

end

Covariance = sigma_cov_tot_TS;
if (Red(2)=='y' || Red(3)=='y'), Cov_reord = CovReord_Mult(Covariance,L1,L2,p); end


% now calculating PID quantities
PID = zeros(length(Red),4);
for r = 1:length(Red)
    if (Red(r) == 'y' && r == 1)
        % MMI
        R = min(MI_TX, MI_TY); 
        U_X = MI_TX - R;
        U_Y = MI_TY - R;
        S = MI - U_Y - U_X - R;
        
        print = "Barrett Redundancy is: " + R + ...
            "\nUnique information of X is: " + U_X + ...
            "\nUnique information of Y is: " + U_Y + ...
            "\nSynergy is: " + S + "\n\n";
%         fprintf(print); 

        PID(r,:) = [U_X, U_Y, R, S];
        
    elseif( Red(r) == 'y' && r == 2 )
        % Idep
        PI = num2cell(calc_pi_Idep_mvn(Cov_reord, [p*L1 p*L2 L]));
        [R, U_X, U_Y, S] = PI{:};

        print = "Idep Redundancy is: " + R + ...
            "\nUnique information of X is: " + U_X + ...
            "\nUnique information of Y is: " + U_Y + ...
            "\nSynergy is: " + S + "\n\n";
%         fprintf(print); 

        PID(r,:) = [U_X, U_Y, R, S];

    elseif( Red(r) == 'y' && r == 3 )
        % Iccs 
        PI = num2cell(calc_pi_mvn(lattice2d(), Cov_reord, ...
                                    [p*L1 p*L2 L], @Iccs_mvn_P2).PI);
        [R, U_X, U_Y, S] = PI{:};
        
        print = "Iccs Redundancy is: " + R + ...
            "\nUnique information of X is: " + U_X + ...
            "\nUnique information of Y is: " + U_Y + ...
            "\nSynergy is: " + S + "\n\n";
%         fprintf(print); 

        PID(r,:) = [U_X, U_Y, R, S];

    else
        PID(r,:) = 0;
    end

end


end