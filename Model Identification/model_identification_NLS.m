% -------------------------------------------------------------------------
% Weighted Sparse Identification of Nonlinear Dynamics (wtSINDy)
%
% This script takes data generated from simulate_NLS.m and applies model
% identification to approximate the true reduced soliton dynamics. The 
% soliton is known to exhibit harmonic oscillation and so the centre of 
% mass of the soliton is governed by the harmonic oscillator x'' = - x.
% This simple ordinary differential equation provides the ground truth for
% comparison of performance across the various methods.
%
% Paper data can be loaded in from centre_mass_data.mat instead of
% generating it using simulate_NLS.m.
% 
% This script accompanies Section 3.3 of Weighted Birkhoff Averages 
% Accelerate Data-Driven Methods. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Load and process data

% load in data instead of generating it using simulate_NLS.m
load centre_mass_data.mat

% Define data matrices for model identification
x = [];
ddx = [];
N = length(xpk);
x = xpk(1:N);
ddx = ddxpk(1:N);

%% Create candidate right-hand-side of ODE 
% --> uses polynomials up to a fixed degree

% Maximum polynomial degree
deg = 5;

Psi = [];
for d = 0:deg
   Psi(d+1,:) = x.^d; 
end

% Sparsity parameter for (weighted) SINDy later
lam = 1e-2;

%% Find coefficients using pseudoinverse

% Unweighted solution
Xi = ddx*pinv(Psi);

% Weight function
w = zeros(1,N);
w(1) = 0;
w(2:end-1) = 142.250375*exp(-1./((2:N-1)/(N).*(1 - (2:N-1)/(N)) ));
w(end) = 0;

% Weighted solution
Xiw = (ddx.*sqrt(w))*pinv(Psi.*sqrt(w));

%% Apply unweighted SINDy

k = 1;
Xis = Xi;
Xis_new = Xis;
while sum(sum(abs(Xis - Xis_new))) > 0  || k == 1 
    
    Xis = Xis_new;
    
    % find small coefficients
    smallinds = (abs(Xis) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xis_new(smallinds) = 0;  
     
    % Identify the elements with large indices to remain in the library
    biginds = ~smallinds;

    % Find coefficients for reduced library
    Xis_new(biginds) = ddx*pinv(Psi(biginds,:)); 

    k = k + 1;
end

Xis = Xis_new;

%% Weighted SINDy

k = 1;
Xiws = Xiw;
Xiws_new = Xiws;
while sum(sum(abs(Xiws - Xiws_new))) > 0  || k == 1 
    
    Xiws = Xiws_new;
    
    % find small coefficients
    smallinds = (abs(Xiws) < lam); 
    
    % Threshold out small coefficients to eliminate from library
    Xiws_new(smallinds) = 0;  
     
    % Identify the elements with large indices to remain in the library
    biginds = ~smallinds;

    % Find coefficients for reduced library
    Xiws_new(biginds) = (ddx.*sqrt(w))*pinv(Psi(biginds,:).*sqrt(w)); 

    k = k + 1;
end

Xiws = Xiws_new;

%% Compare with the exact solution

% Exact coefficient matrix
Xiexact = zeros(1,length(Xi(1,:)));
Xiexact(2) = -1;

% Print errors
Er = norm(Xi - Xiexact) % EDMD error   
Erw = norm(Xiw - Xiexact) % wtEDMD error   
Ers = norm(Xis - Xiexact) % SINDy error
Erws = norm(Xiws - Xiexact) % wtSINDy error

%% Plot results from the paper

% load in paper results
load identification_errors_1.mat % sparsity parameter 1e-2 results
% load identification_errors_2.mat % sparsity parameter 1e-4 results

% Plot results on semilog plot
figure(1)
N = 50:10:5000;
semilogy(N,Er,'Color',[36/255 122/255 254/255],'LineWidth',3)
hold on
semilogy(N,Erw,'Color',[1 69/255 79/255],'LineWidth',3)
semilogy(N,Ers,'Color',[254/255 174/255 0/255],'LineWidth',3)
semilogy(N,Erws,'Color',[0 120/255 0],'LineWidth',3)
xlim([200 5000])
ylim([1e-6 1e-2]) % for 1e-2 sparsity parameter
%ylim([3e-6 1e-2]) % for 1e-4 sparsity parameter
set(gca,'fontsize',20)
xlabel('$N$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
legend('EDMD','wtEDMD', 'SINDy','wtSINDy','FontSize',16,'Location','NorthEast','interpreter','latex')

% For saving the output
% outfilename = 'wtSINDy';
% set(gcf, 'PaperUnits','centimeters');
% set(gcf, 'Units','centimeters');
% pos=get(gcf,'Position');
% set(gcf, 'PaperSize', [pos(3) pos(4)]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
% print('-dpdf',outfilename);