function [U_X, U_Y, R, S, MI] = NormPID(Un_X, Un_Y, Red, Syn, MI)

for k = 1:2

assert(length(MI{k})==length(Un_X{k}) && length(MI{k})==length(Un_Y{k}) ...
       && length(MI{k})==length(Red{k}) && length(MI{k})==length(Syn{k}))

U_X{k} = Un_X{k} ./ MI{k};
U_Y{k} = Un_Y{k} ./ MI{k};
R{k} = Red{k} ./ MI{k};
S{k} = Syn{k} ./ MI{k};
MI{k}(:) = 1;

end