function [Y] = LowFilter(X,ds,fs,fc)

% default parameters
if ~exist('fs', 'var'), fs = 1000; end % sampling frequency 
if ~exist('fc', 'var'), fc = 100; end % frequency cut off

% getting coefficient of the filter
[b, a] = butter(2, fc/(fs/2));

% filtering along rows
Y = filter(b,a,X,[],2); 

% downsampling
Y = Y(1:ds:end);