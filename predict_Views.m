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

meanV168log = mean(v168Log);
stdV168log =  std(v168Log);

Cwo = C;
for i=1:noColumns
    if v168Log(i) > meanV168log+3*stdV168log | v168Log(i) < meanV168log-3*stdV168log
        Cwo = Cwo(setdiff(1:size(Cwo,1),[i]),:);
    end
end

% Correlation coefficients
[noRowsWO,noColumnsWO] = size(Cwo);
vN = Cwo(1:noRowsWO,noColumnsWO);
vNlog = log(vN);
n=24;

corrCoefficients = zeros(1,n);

for i=1:n
    vTemp = Cwo(1:noRowsWO,i);
    vTempLog = log(vTemp);
    corrCoefficients(i) = corr(vTemp,vNlog);
end

% Spliting dataset randomly

CwoLog = log(Cwo);
trainingSet = CwoLog;
l = round(0.1*noRowsWO);
r = randperm(noRowsWO,l);

for i=1:size(r,2)
    for n=1:168
    testSet(i,n) = CwoLog( r(i),n);
    end
    trainingSet = trainingSet(setdiff(1:size(trainingSet,1),[i]),:);
end

% Regression model
y = CwoLog(:,168);
x = CwoLog(:,16);

H = [ones(length(y),1),x]

Astar = inv(H'*H)*H'*y
Ytilde = H*Astar

figure(5);
scatter(x,y);
hold on;
plot(x,Ytilde);

% scatter(x,y);
% hold on;
% [a, b] = polyfit(x,y,1) %finds polynomial that fits best in least square sense
% refline(a(1,1), a(1,2));

% Multiple regression model

%wykasowa� - inf 
noInputs = 5;
x = CwoLog(:,2:noInputs);

H = [ones(length(y),1),x]
Astar = inv(H'*H)*H'*y;

%Relative squared error


