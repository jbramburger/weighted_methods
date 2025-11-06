% -------------------------------------------------------------------------
% Weighted Diffusion Forecasts
%
% This script uses diffusion maps to produce non-parametric probabilitistic
% forecasts for data. The method is augmented to compare the performance
% with and without weights. It is applied to data from the Nino 3.4 index 
% which measures monthly sea-surface temperature anomalies in a region of
% the equatorial Pacific commonly used to characterise El Nino/La Nina 
% variability.
%
% Nino 3.4 index data is stored and loaded from nino34.mat
%
% This script requires a number of scripts to implement diffusion maps and
% produce the non-parametric forecast. These are all stored in the
% accompanying folder "NonparametricModel" and were written by Tyrus Berry
% of George Mason University: https://math.gmu.edu/~berry/
% The exception are the scripts that have a "w" suffix on them as they are
% augemented to include the weighting into the method.
%
% This script accompanies Section 3.5 of Weighted Birkhoff Averages 
% Accelerate Data-Driven Methods. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Add path and load data

% Add a path to the modules created by Tyrus Berry
addpath('NonparametricModel');

% Load Nino 3.4 index data
load nino34

%% Process data

x = X(:,2:end);
allData = reshape(x',12*144,1);
allData = allData - mean(allData);

% Train on historical data from January 1920 to December 1999
Tstart = 51*12; % January 1920
Tend = 130*12; % December 1999
trainingData = allData(Tstart:Tend);

% Verify on data from January 2000 to December 2013
Tvers = 131*12; % January 2000
Tvere = 144*12; % December 2013
Tverl = Tvere - Tvers+1;
testingData = allData(Tvers:Tvere);

% Initializations for diffusion map implementation
obs = testingData; % historical observation
truth = testingData; % groud truth for validation
N = size(trainingData,1);
M = size(trainingData,2);

%% Build non-parametric forecast model

% Initializations for diffusion map and forecast model
forecastSteps = 30; 
neigs = 14; % number of eigenfunctions to use
delays = [1 2 3 4 5 6]; % number of delays for diffusion map eigenfunctions
shifts = 1;
R=0.001;

% Lead time for forecast
leadTime=16;

% unweighted forecast model
[ferr,allForecast,allTrue,allForecastVar] = NonparametricForecast(trainingData,obs,truth,R,forecastSteps,delays,neigs);

% weighted forecast model
[ferrw,allForecastw,allTruew,allForecastVarw] = NonparametricForecastw(trainingData,obs,truth,R,forecastSteps,delays,neigs);

%% Plot RMSE results

figure(1)
t=0:forecastSteps-1;plotInd=1;
plot(t,ferr(:,plotInd),'Color',[36/255 122/255 254/255],'LineWidth',3);
hold on
plot(t,ones(forecastSteps,1)*std(allData(:,plotInd)),'k--','LineWidth',2);
l=legend('Root Mean Square Error','Standard Deviation','Location','best');
set(l,'FontSize',16);
grid on
xlim([1 24]);
set(gca,'fontsize',20)
xlabel('Forecast Steps (months)','Interpreter','latex','FontSize',24,'FontWeight','Bold')
ylabel('Error','Interpreter','latex','FontSize',24,'FontWeight','Bold')
set(gcf,'position',[10 10 600 300])

outfilename = 'diffusion_unweighted_error';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);

figure(2)
plot(t,ferrw(:,plotInd),'Color',[1 69/255 79/255],'LineWidth',3);
hold on
plot(t,ones(forecastSteps,1)*std(allData(:,plotInd)),'k--','LineWidth',2);
l=legend('Root Mean Square Error','Standard Deviation','Location','best');
set(l,'FontSize',16);
grid on
xlim([1 24]);
set(gca,'fontsize',20)
xlabel('Forecast Steps (months)','Interpreter','latex','FontSize',24,'FontWeight','Bold')
ylabel('Error','Interpreter','latex','FontSize',24,'FontWeight','Bold')
set(gcf,'position',[10 10 600 300])

outfilename = 'diffusion_weighted_error';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);

%% Plot correlation results

% Unweighted correlation function
allForecastCentered = allForecast - repmat(mean(allForecast),size(allForecast,1),1);
allTrueCentered = allTrue - repmat(mean(allTrue),size(allTrue,1),1);
corfun = mean(allTrueCentered.*allForecastCentered)./sqrt(mean(allForecastCentered.^2).*mean(allTrueCentered.^2));

% Weighted correlation function
allForecastCenteredw = allForecastw - repmat(mean(allForecastw),size(allForecastw,1),1);
allTrueCenteredw = allTruew - repmat(mean(allTruew),size(allTruew,1),1);
corfunw = mean(allTrueCenteredw.*allForecastCenteredw)./sqrt(mean(allForecastCenteredw.^2).*mean(allTrueCenteredw.^2));


figure(3)
plot(t,corfun,'Color',[36/255 122/255 254/255],'LineWidth',3);
grid on
xlim([1 24]);
ylim([-0.2 1])
set(gca,'fontsize',20)
xlabel('Forecast Steps (months)','Interpreter','latex','FontSize',24,'FontWeight','Bold')
ylabel('Correlation','Interpreter','latex','FontSize',24,'FontWeight','Bold')
set(gcf,'position',[10 10 600 300])

outfilename = 'diffusion_unweighted_correlation';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);

figure(4)
plot(t,corfunw,'Color',[1 69/255 79/255],'LineWidth',3);
grid on
xlim([1 24]);
ylim([-0.2 1])
set(gca,'fontsize',20)
xlabel('Forecast Steps (months)','Interpreter','latex','FontSize',24,'FontWeight','Bold')
ylabel('Correlation','Interpreter','latex','FontSize',24,'FontWeight','Bold')
set(gcf,'position',[10 10 600 300])

outfilename = 'diffusion_weighted_correlation';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);

%% Plot diffusion forecast results

figure(5)
plot(allTrue(:,leadTime),'k','linewidth',3);
hold on;
plot(allForecast(:,leadTime),'Color',[36/255 122/255 254/255],'linewidth',3);
plot(allForecast(:,leadTime)-allForecastVar(:,leadTime),':','Color',[36/255 122/255 254/255],'linewidth',2);
plot(allForecast(:,leadTime)+allForecastVar(:,leadTime),':','Color',[36/255 122/255 254/255],'linewidth',2);
grid on
set(gca,'xtick',(12*(0:1:13)+1));
set(gca,'xticklabel',2000:1:2013);
set(gca,'FontSize',18);
xlim([1 165-12*3-7]);
l=legend('Truth','16 Month Forecast','Forecast Standard Deviation','Location','SouthWest');
set(l,'FontSize',16);   
set(gca,'fontsize',20)
xlabel('Year','Interpreter','latex','FontSize',24,'FontWeight','Bold')
ylabel('El-Ni\~no 3.4 Index','Interpreter','latex','FontSize',24,'FontWeight','Bold')
set(gcf,'position',[10 10 600 300])

outfilename = 'diffusion_unweighted_forecast';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);

figure(6)
plot(allTrue(:,leadTime),'k','linewidth',3);
hold on;
plot(allForecastw(:,leadTime),'Color',[1 69/255 79/255],'linewidth',3);
plot(allForecastw(:,leadTime)-allForecastVarw(:,leadTime),':','Color',[1 69/255 79/255],'linewidth',2);
plot(allForecastw(:,leadTime)+allForecastVarw(:,leadTime),':','Color',[1 69/255 79/255],'linewidth',2);
grid on
set(gca,'xtick',(12*(0:1:13)+1));
set(gca,'xticklabel',2000:1:2013);
set(gca,'FontSize',18);
xlim([1 165-12*3-7]);
l=legend('Truth','16 Month Forecast','Forecast Standard Deviation','Location','SouthWest');
set(l,'FontSize',16);   
set(gca,'fontsize',20)
xlabel('Year','Interpreter','latex','FontSize',24,'FontWeight','Bold')
ylabel('El-Ni\~no 3.4 Index','Interpreter','latex','FontSize',24,'FontWeight','Bold')
set(gcf,'position',[10 10 600 300])

outfilename = 'diffusion_weighted_forecast';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);