% -------------------------------------------------------------------------
% Weighted Approximation of Spectral Measures
%
% This script computes weighted and unweighted autocorrelations to produce
% Fourier modes for spectral measures associated to the Koopman operator of 
% measure-preserving dynamical systems.
%
% Requires Chebfun to build spectral measure. Download at: 
%     https://www.chebfun.org/download/
% 
% This script accompanies Section 3.4 of Weighted Birkhoff Averages 
% Accelerate Data-Driven Methods. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Load and process data

% load cavity flow data
load('cavity_KE.mat')

% process data
X = KE1; % KEj is different reynolds numbers for j=1,2,3,4
c2 = mean(X);
X = X - c2; % subtract mean (i.e., delta function at theta=0)
c1 = sqrt(mean(X.*X));
X = X/c1; % normalise so signal has L^2 norm 1
X = X(:);

%% Compute autocorrelations

% number of moments
N = 1000; 

% initializations
MU = zeros(1,N+1);
MUw = zeros(1,N+1);

% Maximum number of datapoints used
maxDat = 20000;

for j = 0:N

    MU(j+1) = mean(X(1:maxDat-j).*X(1+j:maxDat)); % ergodic averages - add in weights for comparison

    % Weight function
    w = zeros(1,maxDat-j);
    w(1) = 0;
    w(2:end-1) = 142.250375*exp(-1./((2:maxDat-j-1)/(maxDat-j).*(1 - (2:maxDat-j-1)/(maxDat-j)) ));
    w(end) = 0;

    MUw(j+1) = sum(X(1:maxDat-j).*X(1+j:maxDat).*w')/sum(w);

end

%% Build spectral measure densities

% flip autocorrelations to get negative indices
MU=[conj(fliplr(MU(2:end))),MU]; 
MUw=[conj(fliplr(MUw(2:end))),MUw];

% sharp cosine fiter
phi_cosine = @(x) (1+cos(pi*x))/2;
phi_sharp_cosine = @(x) phi_cosine(x).^4.*(35-84*phi_cosine(x)+70*phi_cosine(x).^2-20*phi_cosine(x).^3);
lenMU = (length(MU)-1)/2;
FILTER=phi_sharp_cosine(abs((-lenMU:lenMU)/lenMU));
FILTER(1)=0;
FILTER(end)=0;

mu1 = chebfun(FILTER(:).*MU(:),[-pi pi],'trig','coeffs');
mu2 = chebfun(FILTER(:).*MUw(:),[-pi pi],'trig','coeffs');

%% Plot results

% To produce figures from the paper use:
%    N = 1000 moments
%    maxDat = 1100 snapshots of the kinetic energy
%    j = 1 is periodic data
%    j = 2 is quasiperiodic data
%    j = 3 is skew periodic data (chaotic and quasiperiodic parts)

figure(3)
semilogy(mu1,'Color',[36/255 122/255 254/255],'linewidth',3)
xlabel('$\theta$','interpreter','latex','fontsize',24)
title('Unweighted Density Approximation','Interpreter','latex','FontSize',18,'FontWeight','Bold')
set(gca,'fontsize',20)
set(gcf,'position',[10 10 600 300])
grid on
xlim([0 pi/2])

figure(4)
semilogy(mu2,'Color',[1 69/255 79/255],'linewidth',3)
xlabel('$\theta$','interpreter','latex','fontsize',24)
title('Weighted Density Approximation','Interpreter','latex','FontSize',18,'FontWeight','Bold')
set(gca,'fontsize',20)
grid on
xlim([0 pi/2])
set(gcf,'position',[10 10 600 300])
grid on
xlim([0 pi/2])

%% Plot errors over amount of data

% load saved error data. Change j = 1,2,3 for different Reynolds numbers
load cavity_errors_j=1.mat

datpts = 1100:10:20000;

figure(1)
semilogy(datpts,Error,'Color',[36/255 122/255 254/255],'LineWidth',3)
hold on
semilogy(datpts,Errorw,'Color',[1 69/255 79/255],'LineWidth',3)
axis tight
xlim([1200 19000])
set(gca,'fontsize',20)
xlabel('$N$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
legend('Unweighted','Weighted','FontSize',16,'Location','NorthEast','interpreter','latex')

outfilename = 'cavity_error_skewperiodic';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);