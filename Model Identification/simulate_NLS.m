% -------------------------------------------------------------------------
% Simulate the NLS equation and find centre of mass  
%
% This code simulates the nonlinear Schrodinger (NLS) partial differential 
% equation and extracts the centre of mass at each time step. This data is 
% then used to seed a model identification procedure to identify a 
% reduced-order model for the soliton dynamics in the full NLS. 
% 
% This script accompanies Section 3.3 of Weighted Birkhoff Averages 
% Accelerate Data-Driven Methods. 
%
% Author: Jason J. Bramburger
% -------------------------------------------------------------------------

% Clean workspace
clear all; close all; clc

%% Generate synthetic data from NLS using spectral methods

% Space
L = 5*pi; 
M = 128; 
Y = linspace(-L/2,L/2,M+1); 
x = Y(1:M);
k = (2*pi/L)*[0:M/2-1 -M/2:-1]';

% Time
dt = 0.01; % number of snapshots (N)
t = 0:dt:101;

% Initial condition
shift = 2;
u = 2*sech(2*(x - shift));
ut = fft(u);

% Simulating the NLS
[t, utsol] = ode45(@(t,y) nls_rhs(x,y,k),t,ut);

% Undo fft
for j = 1:length(t)
   usol(j,:) = ifft(utsol(j,:)); % back to x-space
end

%% Plot NLS solution

figure(1)
[Z, T] = meshgrid(x,t);
surf(Z,T,abs(usol))
shading interp
colorbar
set(gca,'FontSize',16,'Xlim',[-5 5], 'Ylim',[t(1),20])
view(0,90)
xlabel('$x$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
ylabel('$t$','Interpreter','latex','FontSize',24,'FontWeight','Bold')

outfilename = 'NLSE';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);

%% Find centre of mass and approximate second derivative

% Centre of mass computation
for j = 1:length(t)
    
   xpk(j) = trapz(x,x.*abs(usol(j,:)))/trapz(x,abs(usol(j,:))); 

end 

% Second derivative approximation
for j = 2:length(t)-1
    
    ddxpk(j-1) = (xpk(j+1) + xpk(j-1) - 2*xpk(j))/(dt^2);
    
end  

% Trim data endpoints
xpk = xpk(2:end-1);

% Plot peaks
figure(2)
subplot(2,1,1)
plot(t(2:end-1),xpk,'k','LineWidth',2)
axis([0 50 -2.1 2.1])
set(gca,'fontsize',16)
xlabel('$t$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
title('Center of Mass Data','Interpreter','latex','FontSize',18,'FontWeight','Bold')

subplot(2,1,2)
plot(t(2:end-1),ddxpk,'k','LineWidth',2)
axis([0 50 -2.1 2.1])
set(gca,'fontsize',16)
xlabel('$t$','Interpreter','latex','FontSize',24,'FontWeight','Bold')
title('Second Derivative Estimation of Center of Mass Data','Interpreter','latex','FontSize',18,'FontWeight','Bold')

outfilename = 'NLSE_centre_mass';

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'Units','centimeters');
pos=get(gcf,'Position');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print('-dpdf',outfilename);

%% NLS Right-Hand-Side

function rhs = nls_rhs(x,ut,k)

    Omega = 1.0;
   
    u = ifft(ut);
  
    rhs = -(1i/2)*(k.^2).*ut + 1i*fft( (abs(u).^2).*u ) - (1i/2)*Omega^2*fft(u.*x'.^2);
end