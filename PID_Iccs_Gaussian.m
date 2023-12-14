
function [UnX, UnY, Red, Syn] = PID_Iccs_Gaussian(Sigma_full,S,T)

PI = num2cell(calc_pi_mvn(lattice2d(), Sigma_full, [S/2 S/2 T], @Iccs_mvn_P2).PI/log2(exp(1)));
[Red, UnX, UnY, Syn] = PI{:};

end