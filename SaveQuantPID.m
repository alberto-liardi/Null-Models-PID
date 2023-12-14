function SaveQuantPID(UnX, UnY, Red, Syn, MI, path, k)

if ~exist('k','var') || (k~=1 && k~=2), error("insert valid model order"); end

app = ["_1","_p"];

% save different p in different files 
% (useful for independent computation of the quantiles)

assert(length(MI)==length(UnX) && length(MI)==length(UnY) ...
   && length(MI)==length(Red) && length(MI)==length(Syn), ...
   "Dimensions of PID atoms do not match!");

save(path+"Quantiles"+app(k)+".mat", "UnX", "UnY", "Red", "Syn", "MI")

end