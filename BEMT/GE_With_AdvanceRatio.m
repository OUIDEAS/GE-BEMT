%==========================================================================
% Jeremy Browne
% 12/18/19
% Looking at Curtis1984 paper on forward flight of a Helicopter blade IGE.
% 
% I also attempted to modify Gina's coeffitcent of thrust ratio equation to
% include the advance ratio for forward flight
%==========================================================================

close all; clear; clc;

R = 0.195; %0.139;              % Prop blade radius - meters
D = 2*R;                % prop blade diameter - meters
z = 0.1:0.1:(2*D);        % range of heights for height ratio - meters
HR = z/D;               % Height ratio using the prop diameter
HR2 = z/R;              % height ratio using the prop radius
m = 1.4;                % mass of the quadcopter - kg - used for thrust
g = 9.81;               % acceleration of gravity - m/s^2

T = (m*g);             % Propeller thrust - Newtons - 1/4 the quad weight
A = pi*R^2;              % Area of the propeller - m^2
rho = 1.225;               % density of air - kg/m^3
vi =  sqrt(T/(2*rho*A)); % induced velocity = velocity to hover

Omega = 287;    % for helicopter
omega = 523;    % for quad copter

lam = vi/(omega*R);       % induced velocity at the rotor - vi/(RPM*R) - m/s   
u = 0:0.01:1;            % range of advance ratios
v_free = omega*R*u;        % free stream velocity - m/s
                         % 523 is 5000 RPM converted to Rad/s
                       
FWrange = v_free/vi;                         
% Coefficents for curve fit
A2 = [0.3101, -3.685, 1.001, -0.002];
c = [1, -1.6];

% Thrust output for forward flight Cheeseman1955
FWratio =  [ 0    0.4  0.8  1.2  1.6  2.0  ];
FWTratio = [ 1.31 1.26 1.16 1.08 1.04 1.005];  

% Thrust output for forward flight Cheeseman1955
% FWratio =  [ 0   0.2 0.4  0.8  1.2  1.6  2.0  ]; % forward velocity
% FWTratio = [ 1.31 1.26 1.16 1.08 1.04 1.005]; % thrust genearation
    
n = 10;
a = 1;
bldLoadNorm = 5;
Ct_bldLoad = 0.0589;
% TgTinf = 1 + (0.25*(n)*(bldLoadNorm))/(sqrt(Ct/bldLoad)*(1/(16*(z/R)^2*(1_+(mu)^2))))


for i = 1:length(z)
  
   for j = 1:length(u)
    % Coefficent of thrust in forward flight includeing GE - From Literature
    ctGEFW(i, j) = 1.00 / (1.00 - ( ( R / (4.00*z(i)) )^2.00 / (1.00 + ( v_free(j) / (vi) )^2.00)) );
%     (vi-0.02*j^2)
    % Coefficent of thrust in hover includeing GE  - From Literature
    ctGE(i, j) = 1.00 / (1.00 - ( R / (4.00*z(i)) )^2.00);
    
    % Coefficent of thrust in for forward flight using my curve fit 
    ctGECorr(i,j) = A2(1)*exp(A2(2)*HR(i)) + A2(3)*exp(A2(4)*HR(i));    % + c(1)*exp(c(2)*HR(i));
   
    % /(1+(u(j)/(3*lam))^2) 
    
    % Thrust output for forward flight Cheeseman1955
    TgTinf(i,j) = 1 + (0.25*n*a*sqrt(0.12)) / ( sqrt(Ct_bldLoad)) * [ 1 / ( 16*(0.5)^2*(1 + (u(j))^2.00)) ] ;
    
   end
    
    
end

% Differences between the first two models
ModelDiff = ctGEFW - ctGE; 

% figure
% plot(u, ctGEFW(30,:))
% grid on; 
% xlabel('u'); ylabel('Ct');

figure
% hold on;
x = HR2;
y = FWrange;
z = ctGEFW;
z2 = ctGE;
surf(x, y, z')
% surf(x, y, z2)
xlabel('HR'); ylabel('V_{inf} / v_i'); zlabel('T_g / T_{inf}');
xlim([0 2.5]); ylim([0 2.5]); zlim([1 1.35]);
% hold off; 
grid on; 
title('Thrust Ratio with GE and Forward flight')
colorbar;

% figure
% x = HR;
% y = u;
% z = ctGE;
% surf(x, y, z')
% xlabel('HR'); ylabel('u'); zlabel('Ct');
% title('Thrust ratio with only GE')
% colorbar;
% 
% figure
% x = HR;
% y = u;
% z = ctGECorr;
% surf(x, y, z')
% xlabel('HR'); ylabel('u'); zlabel('Ct');
% title('Thrust ratio using my Thrust Ratio Curve Fit for GE')
% colorbar;
% 
% % figure
% % x = FWratio;
% % y = FWTratio;
% % plot(x,y)
% % xlabel('Vi/vt'); ylabel('Tg/Tinf');
% % title('Variation of thrust ratio with forward speed')
% 
% figure
% x = HR;
% y = u;
% z = TgTinf;
% surf(x, y, z')
% xlabel('HR'); ylabel('u'); zlabel('Ct');
% title('?????')
% colorbar;

figure
x = FWrange;
y = ctGEFW([1:4],:);
x2 = FWratio;
y2 = FWTratio;

hold on;
plot(x, y)
scatter(x2, y2)
xlabel('V_{inf} / vi'); ylabel('T_g / T_{inf}');
legend('HR = 0.5', 'HR = 1.0', 'HR = 1.5', 'HR = 2.0' )
grid;
ylim([1.0, 1.35]); xlim([0, 2]);
% title('FW flight?')
