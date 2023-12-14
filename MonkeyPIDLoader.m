function [Un_X, Un_Y, Red, Syn, MI] = MonkeyPIDLoader(file_name)

data = readmatrix(file_name);

Red = {data(:,1)', data(:,2)'};  
Syn = {data(:,3)', data(:,4)'};
Un_X = {data(:,5)', data(:,6)'};
Un_Y = {data(:,7)', data(:,8)'};

MI = { Red{1} + Syn{1} + Un_X{1} + Un_Y{1}, ...
       Red{2} + Syn{2} + Un_X{2} + Un_Y{2} };

end