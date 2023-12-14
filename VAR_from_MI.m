%% Get random VAR system with specific Mutual Information
% MI - total mutual information between sources and targets
% p - model order of VAR(p)
% Sources - Number of sources


function [Coef,V,flag,g] = VAR_from_MI(MI,p,S)
    
    opts = optimset('display', 'none');
    Coef = cell(1,p);
    % initial conditions for the optimisation
    x0 = [0 -1 2 -2 5 -5];
    flag=0;

    % sample the matrices
    % NB. for Wishart degrees of freedom must be > size - 1
    
    % random VAR matrix
    A = var_rand(S,p,rand(1));
    while any(abs(A(:))<1e-18)
        A = var_rand(S,p,rand(1));
    end
    
    % target conditional covariance
    V = wishrnd(eye(S),S+1);
    %         A_old = A; V_old=V;
    %         [A,V] = var_normalise(A,V);
    
    % find g such that mutual information is MI_tot
    l = 1;
    %         t0 = tic();
    [g, ~, code, out] = fzero(@(x) fun(x,A,V,MI), x0(l), opts);
      
    while( code ~= 1 && l+1 <= length(x0) )
        l = l+1;
        [g, ~, code, out] = fzero(@(x) fun(x,A,V,MI), x0(l), opts);
    end
    if( code ~= 1 ), fprintf("optimisation failed "+sprintf('%01d',code)+"\n"); 
         flag=1; return; 
    end 

%         fprintf('%f\t%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t', MI, p, l, g, gg, out.iterations, ...
%            out.intervaliterations, out.funcCount, cc); toc(t0);
    
    % now do one more calculation
    gg = Sigmoid(g,1);
    [B,~] = specnorm(A,gg);
    for s = 1:p
        Coef{s} = B(:,:,s);
    end
    %         Coef = num2cell(B,[1,2]);
    %         Coef = squeeze(Coef);
    if ~isfinite(fun(g,B,V,MI))
        disp("Optimisation function returned NaN value!");
        flag=1;
        return;
    end

end


function y = fun(x,A,V,I)

    x = Sigmoid(x,1);
    [BB,~] = specnorm(A,x);

    % sometimes for large p and specific (small) values of the spectral radius 
    % the Lyapunov equation does not have a unique solution
    try
        G = var_to_autocov(BB,V,0);
        MI_value = (0.5*log(det(G(:,:,1)))-0.5*log(det(V)))/log(2);
        if ~isreal(MI_value)
            % when the spectral radius is close to 1 the MI can become complex
            error("error in optimising the g!");
        end
        y = MI_value-I;

    catch 
        y=NaN;
        disp("NAN value encountered");

    end
%     disp(y);

end

function y = Sigmoid(x,M)
    
    y = M/(1+exp(-x));

end