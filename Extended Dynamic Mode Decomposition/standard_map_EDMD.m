% -------------------------------------------------------------------------
% Weighted Extended Dynamic Mode Decomposition (wtEDMD)  
%
% This code compares weighted extended dynamic mode decomposition (wtDMD) 
% with extended dynamic mode decomposition on data gathered from the 
% standard map. The goal is to show that building weights into EDMD can 
% provide better approximations of the EDMD matrix with fewer data points. 
% 
% This script accompanies Section 3.2 of Weighted Birkhoff Averages 
% Accelerate Data-Driven Methods. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Generate Data

% Number of datapoints
N = 1e6;

% nonlinearity parameter
lambda = 0.25;

dat = zeros(N,2);
dat(1,:) = 2*pi*rand(2,1)';
for n = 2:N+1
    %lambda = rand;
    dat(n,1) = dat(n-1,1) + lambda*sin(dat(n-1,2));
    dat(n,2) = dat(n-1,2) + dat(n-1,1) + lambda*sin(dat(n-1,2)); 
    dat(n,:) = mod(dat(n,:),2*pi);
end

%% Build complex observable matrices

% Select wavenumbers for Fourier basis
maxk = 1;
[kx, ky] = meshgrid(-maxk:maxk,-maxk:maxk);
kx = kx(:);
ky = ky(:);

Phi = [];
Psi = [];
for j = 1:length(kx)

    Phi(j,:) = exp( 1i*( kx(j)*dat(1:N,1) + ky(j)*dat(1:N,2) ) ); 
    Psi(j,:) = exp( 1i*( kx(j)*dat(2:N+1,1) + ky(j)*dat(2:N+1,2) ) ); 

end

%% Weight function and weighted observable matrices
    
w = zeros(1,N);
w(1) = 0;
w(2:end-1) = 142.250375*exp(-1./((2:N-1)/(N).*(1 - (2:N-1)/(N)) ));
w(end) = 0;

Phiw = Phi.*sqrt(w);
Psiw = Psi.*sqrt(w);

%% Construct EDMD matrices

K = Phi*pinv(Psi);
Kw = Phiw*pinv(Psiw);

%% Recreate figures from text

load standard_loop_lambda=quarter.mat
%load standard_loop_lambda=5.mat
%load standard_loop_lambda=random.mat

figure(2)
datapts = 100:1000:1e6;

% Plot the means of the datasets
loglog(datapts,mean(Ker),'Color',[36/255 122/255 254/255],'LineWidth',3)
hold on
loglog(datapts,mean(Kerw),'Color',[1 69/255 79/255],'LineWidth',3)

xlabel('$N$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
legend('Standard EDMD','Weighted EDMD','FontSize',16,'Location','NorthEast','interpreter','latex')
set(gca,'fontsize',20)

% For printing the figure
% outfilename = 'edmd_error_random';
% 
% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'Units','centimeters');
% pos=get(gcf,'Position');
% set(gcf, 'PaperSize', [pos(3) pos(4)]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
% print('-dpdf',outfilename);