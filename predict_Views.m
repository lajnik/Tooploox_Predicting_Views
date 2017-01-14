clear all;
close all;

%% Read data
fileID = fopen('data.csv');
C1 = textscan(fileID,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',');
fclose(fileID);
C2 = C1(1, 2:169);
C = cell2mat(C2);
[noRows,noColumns] = size(C);
v168 = C(1:noRows,noColumns);

%% Basic statistics
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

%% Distribution of the log transformed v(168)
figure(3);
v168Log = log(v168);
normplot(v168Log);
figure(4);
histfit(v168Log);
grid on;
% Conclusion: distribution of data is log-normal distribution

%% Removing outliers

meanV168log = mean(v168Log);
stdV168log =  std(v168Log);

Cwo = C;
for i=1:noColumns
    if v168Log(i) > meanV168log+3*stdV168log | v168Log(i) < meanV168log-3*stdV168log
        Cwo = Cwo(setdiff(1:size(Cwo,1),[i]),:);
    end
end

%% Correlation coefficients
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

%% Spliting dataset randomly

[r c] = find(Cwo==0);
for i=1:size(r,1)
    Cwo(r(i),1) = 0.1;
end

CwoLog = log(Cwo);

trainingSet = CwoLog;
l = round(0.1*noRowsWO);
r = randperm(noRowsWO,l);

for i=1:size(r,2)
    for n=1:168
    testSet(i,n) = CwoLog(r(i),n);
    end
    trainingSet = trainingSet(setdiff(1:size(trainingSet,1),i),:);
end

%% Regression model
y = trainingSet(:,168); %traing set

for i = 1:167
    x = trainingSet(:,i);
    H = [ones(length(y),1),x];
    AstarVal = inv(H'*H)*H'*y;
    Astar(1,i)=AstarVal(1,1);
    Astar(2,i)=AstarVal(2,1);
end;

Ytilde = H*Astar(:,167);
figure(5);
scatter(x,y);
hold on;
plot(x,Ytilde);

% scatter(x,y);
% hold on;
% [a, b] = polyfit(x,y,1) %finds polynomial that fits best in least square sense
% refline(a(1,1), a(1,2));
 
%% Multiple regression model

AstarMulti = zeros(168);

for i=1:167
    noInputs = i;
    x = trainingSet(:,1:noInputs);

    H = [ones(length(y),1),x];
    AstarMultiVal = inv(H'*H)*H'*y;
    for n =1:noInputs+1
        AstarMulti(n,i)=AstarMultiVal(n,1);
    end
   
end

%% Relative squared error

T = size(testSet,1);

%prediction of views in 16th hour
x = testSet(:,16); 
H = [ones(length(x),1),x];

Ytilde = H*Astar;
Yreal = testSet(:,168); 

mSRE = 0;
for i = 1:T
    mSRE = mSRE + (Ytilde(i)/Yreal(i) - 1)^2;
end
mSRE = mSRE/T;

%% Plotting mSREs
% Single input

noHours = 24
for n=1:noHours
    x = testSet(:,n); 
    H = [ones(length(x),1),x];

    AstarTemp = Astar(1:2, n);
    Ytilde = H*AstarTemp;
    Yreal = testSet(:,n); 

    mRSE = 0;
    for i = 1:T
        mRSE = mRSE + (Ytilde(i)/Yreal(i) - 1)^2;
    end
    mRSE = mRSE/T;
    mRSESingle(n) = mRSE;
end

figure(9);
hours = 1:1:noHours;
plot( hours, mRSESingle, 'b*-');
grid on;
hold on;

% Multiple input

for n=1:noHours
    noInputs = n;
    x = testSet(:,1:noInputs);

    H = [ones(length(x),1),x];
    AstarTemp = AstarMulti(1:n+1,n);
    
    Ytilde = H*AstarTemp;
    Yreal = testSet(:,n); 

    mRSE = 0
    for i = 1:T
        mRSE = mRSE + (Ytilde(i)/Yreal(i) - 1)^2;
    end
    mRSE = mRSE/T
    mRSEMulti(n) = mRSE;
end

plot(hours, mRSEMulti, 'r*-');
 
 
 
