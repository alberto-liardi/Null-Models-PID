rng('default') % for reproducibility of random variables
disp("begin run");
gate = "XOR";


if gate == "XOR"
    boolfun = @XORgate;
elseif gate == "OR"
    boolfun = @ORgate;
end 

% ps = [0.8147, 0.1678, 0.0022, 0.0152];

MIx=99; MIy=99; i=0;
while MIx>0 && MIy>0
    ps1 = rand(); ps2=(1-ps1)*rand(); ps3=(1-ps1-ps2)*rand(); ps4=1-ps1-ps2-ps3;
    ps = [ps1, ps2, ps3, ps4];
    pe = 0.000001;
    
    [pz0,pz1] = boolfun(ps);
    
    pt0 = pz0*(1-pe) + pz1*pe;
    pt1 = pz1*(1-pe) + pz0*pe; % 1-pt0;
    
    HT = -pt0*log(pt0)-pt1*log(pt1);
    HTcondZ = -pe*log(pe)-(1-pe)*log(1-pe);
    % HTcondZ = 0;
    MI = HT - HTcondZ;
    
    [MIx, MIy] = DiscreteMarginalMI(ps,pe,HT,gate);
    
    i = i+1;
    if mod(i,10000)==0, disp(i); end
end

disp("found some!")

function [pz0,pz1] = XORgate(ps)
    pz0 = ps(1)+ps(4);
    pz1 = ps(2)+ps(3);
end

function [pz0,pz1] = ORgate(ps)
    pz0 = ps(1);
    pz1 = 1-ps(1);
end