rng('default') % for reproducibility of random variables

disp(" "); disp("---- Start of compilation ----"); disp(" ");

% loading data (if necessary)
if exist('data',     'var') == 0
    disp("Loading data...");
    
    % !!! IMPORTANT: Monkey data must be in a folder called "Monkey_data"
    data = KTMDLoader("Kin2", "W", "Monkeys_data");
    
    disp("... done.");
end

% setting parameters of the data
tot_epoques = length(data);
tot_channels = size(data{1},2);
tot_period = size(data{1},1);

% setting parameter for the analyses
epoques = 1;

% choosing 2 different random channels and some random epoques
c = randi([1 tot_channels],1);
e = randi([1 tot_epoques],1,epoques);

X = zeros(1, tot_period, epoques);

for n = 1:1:size(X,3)
        X(:,:,n) = data{e(n)}(:,c(1))';
end

[p,f] = PowerSpec(X);

disp(" "); disp("---- End of compilation ----"); disp(" ");



function [p,f] = PowerSpec(X)

    freqs = 1:1:100; Fs = 1e3;
    
    for n = 1:1:size(X,3)
        [p, f] = pwelch(X(:,:,n), [], [], freqs, Fs);
        figure(); semilogy(f,p);
        xlabel("Frequency"); ylabel("Power spectrum");
    end

end