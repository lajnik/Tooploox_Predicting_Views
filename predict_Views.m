clear all;
close all;

% Read data
fileID = fopen('data.csv');
C1 = textscan(fileID,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',');
fclose(fileID);
C2 = C1(1, 2:169);
C = cell2mat(C2);
[noRows,noColumns] = size(C);
v168 = C(1:noRows,noColumns);

% Basic statistics
% Arithmetic mean
meanVn = mean(v168);

% Geometric mean
geoMeanVn = geomean(v168);

% Mode
modeVn = mode(v168);

% Median
medianVn = median(v168);

% Quartiles
quartileVn1 = quantile(v168,0.25);
quartileVn3 = quantile(v168,0.75);

% Standard deviation
stdVn =  std(v168);

% Variance
varVn = stdVn^2;

% Distribution of v(168)
figure(1);
normplot(v168);
figure(2);
histfit(v168);
grid on;

% Distribution of the log transformed v(168)
figure(3);
v168Log = log(v168);
normplot(v168Log);
figure(4);
histfit(v168Log);
grid on;
% Conclusion: distribution of data is log-normal distribution

% Removing outliers

