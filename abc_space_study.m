

I_list = 1;%:0.5:3;
B = -0.99:0.1:0.99;
A_true = []; C_true = []; B_true = [];

for i = 1:length(I_list)

    A = -sqrt(1-exp(-2*I)):0.1:sqrt(1-exp(-2*I));
    C = -sqrt(1-exp(-2*I)):0.1:sqrt(1-exp(-2*I));
    
    for n = 1:length(A)
        for m = 1:length(C)
            for l = 1:length(B)

                a = A(n); b = B(l); c = C(m);

                % calculate I and check the difference with the chosen onw
                if(1-(a^2+b^2+c^2)+2*a*b*c<1e-12)
%                     count = count + 1; 
                    continue
                end

                I = 0.5*log((1-b^2)/(1-(a^2+b^2+c^2)+2*a*b*c));
                if abs(I-I_list(i))>1e5
                    continue
                end
                
                A_true(end+1) = a;
                B_true(end+1) = b;
                C_true(end+1) = c;

            end
        end
    end
end

disp("done");
%%
figure();
scatter3(real(A_true),real(C_true),real(B_true),'.');
xlabel("a");
ylabel("c");
zlabel("b");

%% --------------------------------- second try



I_list = 1;%:0.5:3;
B = -0.99:0.1:0.99;
A_true = []; C_true = []; B_true = [];

for i = 1:length(I_list)

    A = -sqrt(1-exp(-2*I)):0.1:sqrt(1-exp(-2*I));
    C = -sqrt(1-exp(-2*I)):0.1:sqrt(1-exp(-2*I));
    
    for n = 1:length(B)
        for m = 1:length(C)
            for l = 1:length(A)

                b = B(n); a = A(l); c = C(m);

                % calculate I and check the difference with the chosen onw
                if(1-(a^2+b^2+c^2)+2*a*b*c<1e-12)
%                     count = count + 1; 
                    continue
                end

                I = 0.5*log((1-b^2)/(1-(a^2+b^2+c^2)+2*a*b*c));
                if abs(I-I_list(i))>1e5
                    continue
                end
                
                A_true(end+1) = a;
                B_true(end+1) = b;
                C_true(end+1) = c;

            end
        end
    end
end

disp("done");
%%
figure();
scatter3(real(A_true),real(C_true),real(B_true),'.');
xlabel("a");
ylabel("c");
zlabel("b");