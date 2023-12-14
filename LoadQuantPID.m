function [UnX, UnY, Red, Syn, MI] = LoadQuantPID(file)

if ~isfile(file), error(file+" directory has not been found!"); end

data = load(file);
UnX = data.UnX;     UnY = data.UnY;
Red = data.Red;     Syn = data.Syn;
MI = data.MI;
