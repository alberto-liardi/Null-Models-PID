
function [UnX, UnY, Red, Syn] = PID_Idep_Gaussian(Sigma_full,S,T)

    PI = num2cell(calc_pi_Idep_mvn(Sigma_full, [S/2 S/2 T])/log2(exp(1)));
    [Red, UnX, UnY, Syn] = PI{:};

end