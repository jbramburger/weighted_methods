% -------------------------------------------------------------------------
% Weighted Dynamic Mode Decomposition (wtDMD)  
%
% This code compares weighted dynamic mode decomposition (wtDMD) with
% standard dynamic mode decomposition on fluid flow around a cylinder data. 
% The goal is to show that building weights into DMD can provide better 
% approximations of the DMD matrix and its eigenvalues using fewer 
% datapoints. 
% 
% The for this demonstration is stored in Vorticity_data.mat and was
% provided by Matthew Colbrook. We also use pre-saved random orthogonal 
% vectors to project the data onto, which is stored in rand_ortho.mat. 
%
% This script accompanies Section 3.1 of Weighted Birkhoff Averages 
% Accelerate Data-Driven Methods. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Load and process data

% load in data
load rand_ortho.mat
load Vorticity_data.mat

% global variables
r = 11; % target rank. Paper used r = 11 and r = 21
U = U(:,1:r); % random orthogonal vectors for projection
M = 500; % number of data points
ind = (1:M);

% Arrange vorticity data into snapshot matrices
X = VORT(:,ind); 
Y = VORT(:,ind+1);

%% Define the bump function for weighted method

w = zeros(1,M);
w(1) = 0;
w(2:end-1) = 142.250375*exp(-1./((2:M-1)/(M).*(1 - (2:M-1)/(M)) ));
w(end) = 0;

%% Weighted and standard DMD algorithms

% Project data onto orthogonal modes
PX = X'*U;
PY = Y'*U;

% Create weighted snapshot matrices
PXw = sqrt(w').*PX;
PYw = sqrt(w').*PY;

% DMD matrices
K = pinv(PX)*PY; % unweighted
Kw = pinv(PXw)*PYw; % weighted

%% Extract DMD eigenvalues

[V,LAM] = eig(K,'vector');
[Vw,LAMw] = eig(Kw,'vector');

% Sort unweighted eigenvalues
[~,I] = sort(abs(1-LAM),'ascend'); % reorder modes
LAM = LAM(I);

% Sort weighted eigenvalues
[~,Iw] = sort(abs(1-LAMw),'ascend'); % reorder modes
LAMw = LAMw(Iw);

%% Compute errors against long time results

% load in long time results
if r == 11
    load true_eigs_r=11.mat
elseif r == 21 
    load true_eigs_r=21.mat
end

% DMD matrix errors
KEr = norm(K - longK)/norm(longK)
Kerw = norm(Kw - longKw)/norm(longKw)

% Eigenvalue error
EigEr = norm(LAM - longLAM)/norm(longLAM)
EigErw = norm(LAMw - longLAMw)/norm(longLAMw)

%% Plot long-time eigenvalues

figure(1)
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-k')
hold on
plot(real(longLAMw),imag(longLAMw),'.','Color',[1 69/255 79/255],'markersize',30)
axis equal
axis([-1.15,1.15,-1.15,1.15])

% title('EDMD Eigenvalues','interpreter','latex','fontsize',18)
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
set(gca,'fontsize',20)

%% Recreate plots from the paper

% load in results
if r == 11
    load results_r=11.mat
elseif r == 21 
    load results_r=21.mat
end

figure(2)
datapts = 10:10:500;
semilogy(datapts,EigEr,'Color',[36/255 122/255 254/255],'LineWidth',3)
hold on
semilogy(datapts,EigErw,'Color',[1 69/255 79/255],'LineWidth',3)
set(gca,'fontsize',20)
xlim([datapts(1) datapts(end)])
title('Eigenvalue Error','interpreter','latex','fontsize',20)
xlabel('$N$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
legend('Standard DMD','Weighted DMD','FontSize',16,'Location','SouthWest','interpreter','latex')

outfilename = 'DMD_eig_error_r=11';

% getting the size right
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);

figure(3)
semilogy(datapts,KEr/norm(longK),'Color',[36/255 122/255 254/255],'LineWidth',3)
hold on
semilogy(datapts,Kerw/norm(longKw),'Color',[1 69/255 79/255],'LineWidth',3)
set(gca,'fontsize',20)
xlim([datapts(1) datapts(end)])
title('DMD Matrix Error','interpreter','latex','fontsize',20)
xlabel('$N$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
legend('Standard DMD','Weighted DMD','FontSize',16,'Location','SouthWest','interpreter','latex')

outfilename = 'DMD_matrix_error_r=11';

% getting the size right
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);