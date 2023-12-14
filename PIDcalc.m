function [ U_X,U_Y,R,S,Gamma ] = PIDcalc(p,A,V,target,Red)

% ---------------------------------------------------------------------- %
% Compute PID quantities with the relation between Entropy - MI  %
% ---------------------------------------------------------------------- %

if ~exist('V','var'), V = eye(2); end
assert(iscell(A), "Coefficient matrices must be in a 1xp cell");
assert(0<p, "Enter a valid model order"); 
assert(length(A)==p, "Model order and number of evolution matrices are not the same");

% creating the matrices for the Lyapunov equation

% coefficient matrix
A_L = horzcat( eye(2*(p-1)), zeros(2*(p-1),2) ); 
A_L = vertcat( cell2mat(A) , A_L);
% disp("Coefficient evolution matrix: "); disp(A_L);

% covariance of random noise variables matrix
Lambda_L = horzcat( V, zeros(2, 2*(p-1)) ); 
Lambda_L = vertcat( Lambda_L , zeros(2*(p-1), 2*p) );    
% disp("Noise covariance matrix: "); disp(Lambda_L);

% solving Lyapunov equation for Gamma: 
Gamma = dlyap(A_L,Lambda_L);

% deconstructing and reordering gamma matrices
split = zeros(1,p)+2;
Gamma = Gamma([1 2],:); % only the first line is linear independent
Gamma = mat2cell(Gamma, 2, split);


% now find last gamma matrix (Gamma_p) from YW equation:
Gamma_new = 0;
for l = 1:1:p
    Gamma_new = Gamma_new + A{l}*Gamma{p-l+1};
end
Gamma{end+1} = Gamma_new; 
% disp("Gamma matrices: "); 
% for g = 1:1:p+1
%     print="Gamma"+g+"\n"; fprintf(print);
%     disp(Gamma{g});
% end


% can now construct Sigma matrices

if( target == "single" )

    % T is the target, S the source, TS both of them
    sigma_cov_TS = cell(p+1,p+1);
    sigma_cov_TS(:) = {[0 0; 0 0]};
    for s = 1:1:p+1
        for t = s:1:p+1
            sigma_cov_TS{s,t} = Gamma{t-s+1};
        end
    end
    % going from the upper diagonal to the full matrix 
    sigma_cov_TS = cell2mat(sigma_cov_TS);
    sigma_cov_TS = sigma_cov_TS - tril(sigma_cov_TS,-1);
    sigma_cov_tot_TS= triu(sigma_cov_TS,1) + sigma_cov_TS';

    sigma_cov_tot_S = sigma_cov_tot_TS(1:end-2, 1:end-2);
    
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
    fprintf("\nMutual Information computed using entropies is:"); disp(MI);
    
    
    % Now to calculate single Mutual Information brought by X, Y_past, need to
    % calculate some more covariance matrices (cov(X_n,X_n), cov(Y_n,Y_n))
    
    % covariance matrices
    sigma_X = zeros(p+1,p+1);
    sigma_Y = zeros(p+1,p+1);
    for s = 1:1:p+1
        for t = s:1:p+1
            sigma_X(s,t) = Gamma{t-s+1}(1,1);
            sigma_Y(s,t) = Gamma{t-s+1}(2,2);
        end
    end
    % (can probaby be written more efficiently by tayloring the first row
    % adequately)
    
    % concatenating matrices to obtain the full covariance matrix (_tot)
    sigma_X_tot = sigma_X' + triu(sigma_X,1);
    sigma_Y_tot = sigma_Y' + triu(sigma_Y,1);
    
    row = zeros(1,p+1);
    for l = 1:1:p+1
        row(l) = Gamma{l}(2,1);
    end
    sigma_YX_tot = vertcat( row, sigma_X_tot );
    column = vertcat( Gamma{1}(2,2), row' );
    sigma_YX_tot = horzcat( column, sigma_YX_tot);
    
    sigma_Y_past = sigma_Y_tot(2:end,2:end); 
    
    % calculating I(Y,X_tot):
    H_T = .5*log(det(2*pi*exp(1) * sigma_cov_tot_T));
    H_S = .5*log(det(2*pi*exp(1) * sigma_X_tot));
    H_TS = .5*log(det(2*pi*exp(1) * sigma_YX_tot));
    
    MI_TX = H_T + H_S - H_TS;
    fprintf("Mutual Information brought by X_tot is: " + MI_TX);
    
    
    % calculating I(Y,Y_past):
    H_T = .5*log(det(2*pi*exp(1) * sigma_cov_tot_T));
    H_S = .5*log(det(2*pi*exp(1) * sigma_Y_past));
    H_TS = .5*log(det(2*pi*exp(1) * sigma_Y_tot));
    
    MI_TY = H_T + H_S - H_TS;
    
    fprintf("\nMutual Information brought by Y_past is:"); disp(MI_TY);

elseif( target == "joint")

    % T is the target, S the source, TS both of them
    sigma_cov_TS = cell(p+1,p+1);
    sigma_cov_TS(:) = {[0 0; 0 0]};
    for s = 1:1:p+1
        for t = s:1:p+1
            sigma_cov_TS{s,t} = Gamma{t-s+1};
        end
    end
    % going from the upper diagonal to the full matrix 
    sigma_cov_TS = cell2mat(sigma_cov_TS);
    sigma_cov_TS = sigma_cov_TS - tril(sigma_cov_TS,-1);
    sigma_cov_tot_TS= triu(sigma_cov_TS,1) + sigma_cov_TS';
    
    sigma_cov_tot_S = sigma_cov_tot_TS(1:end-2, 1:end-2);

    sigma_cov_tot_T = Gamma{1};

    % calculating entropies
    H_T = .5*log(det(2*pi*exp(1) * sigma_cov_tot_T));
    H_S = .5*log(det(2*pi*exp(1) * sigma_cov_tot_S));
    H_TS = .5*log(det(2*pi*exp(1) * sigma_cov_tot_TS));
    
    % mutual information
    MI = H_T + H_S - H_TS;
    fprintf("\nMutual Information computed using entropies is:"); disp(MI);
    
    
    % Now to calculate single Mutual Information brought by X_past, Y_past, need to
    % calculate some more covariance matrices
    
    % covariance matrices
    sigma_X = zeros(p+1,p+1);
    sigma_Y = zeros(p+1,p+1);
    for s = 1:1:p+1
        for t = s:1:p+1
            sigma_X(s,t) = Gamma{t-s+1}(1,1);
            sigma_Y(s,t) = Gamma{t-s+1}(2,2);
        end
    end
    % concatenating matrices to obtain the full covariance matrix (_tot)
    sigma_X_tot = sigma_X' + triu(sigma_X,1);
    sigma_Y_tot = sigma_Y' + triu(sigma_Y,1);

%     disp(sigma_X_tot);
%     disp(sigma_Y_tot);

    row_YX = zeros(1,p+1);
    row_XY = zeros(1,p+1);
    for l = 1:1:p+1
        row_YX(l) = Gamma{l}(2,1);
        row_XY(l) = Gamma{l}(1,2);
    end
    sigma_YX_tot = vertcat( row_YX, sigma_X_tot );
    column = vertcat( Gamma{1}(2,2), row_YX' );
    sigma_YX_tot = horzcat( column, sigma_YX_tot);

    sigma_XY_tot = vertcat( row_XY, sigma_Y_tot );
    column = vertcat( Gamma{1}(1,1), row_XY' );
    sigma_XY_tot = horzcat( column, sigma_XY_tot);

%     disp(sigma_YX_tot);
%     disp(sigma_XY_tot);
     
    sigma_Y_past = sigma_Y_tot(2:end,2:end);
    sigma_X_past = sigma_X_tot(2:end,2:end);

%     disp(sigma_Y_past);
%     disp(sigma_X_past);

    % calculating I(XY,X_past):
    H_T = .5*log(det(2*pi*exp(1) * sigma_cov_tot_T));
    H_S = .5*log(det(2*pi*exp(1) * sigma_X_past));
    H_TS = .5*log(det(2*pi*exp(1) * sigma_YX_tot));
    
    MI_TX = H_T + H_S - H_TS;
    fprintf("Mutual Information brought by X_tot is: " + MI_TX);
    
    
    % calculating I(XY,Y_past):
    H_T = .5*log(det(2*pi*exp(1) * sigma_cov_tot_T));
    H_S = .5*log(det(2*pi*exp(1) * sigma_Y_past));
    H_TS = .5*log(det(2*pi*exp(1) * sigma_XY_tot));
    
    MI_TY = H_T + H_S - H_TS;
    
    fprintf("\nMutual Information brought by Y_past is:"); disp(MI_TY);

end

Covariance = sigma_cov_tot_TS;

% now calculating PID quantities

if (Red == "MMI")
    R = min(MI_TX, MI_TY); 
    U_X = MI_TX - R;
    U_Y = MI_TY - R;
    S = MI - U_Y - U_X - R;
elseif(Red == "Idep")
    [R, U_X, U_Y, S] = IdepRed(Covariance, [p p 2]);
elseif(Red == "Iccs")
    sources = {[1:p], [1:p]};
    R = Iccs_mvn_P2(sources, Covariance, [p p 2]);
end

print = "Barrett Redundancy is: " + R + ...
        "\nUnique information of X is: " + U_X + ...
        "\nUnique information of Y is: " + U_Y + ...
        "\nSynergy is: " + S + "\n";
fprintf(print); 

end