%==========================================================================
% Jeremy Browne
% 12/18/19
% Looking at Curtis1984 paper on forward flight of a Helicopter blade IGE.
% 
% I also attempted to modify Gina's coeffitcent of thrust ratio equation to
% include the advance ratio for forward flight
%==========================================================================

close all; clear; clc;

R = 2/2;              % Prop blade radius - meters
D = 2*R;                % prop blade diameter - meters
z = 0.1:1:(2*R);        % range of heights for height ratio - meters
HR = z/R;              % height ratio using the prop radius

M = 2540;
g = 9.81;               % acceleration of gravity - m/s^2

T = (M*g);             % Propeller thrust - Newtons 
A = pi*R^2;              % Area of the propeller - m^2
rho = 1.225;               % density of air - kg/m^3
vi =  sqrt(T/(2*rho*A)); % induced velocity = velocity to hover

Omega = 287;    % for helicopter

lam = vi/(Omega*R);       % induced velocity at the rotor - vi/(RPM*R) - m/s   
u = 0:0.1:2;             % range of advance ratios
v_free = Omega*R*u;      % free stream velocity - m/s
                         % 523 is 5000 RPM converted to Rad/s
                         
% Coefficents for curve fit
A2 = [0.3101, -3.685, 1.001, -0.002];
c = [1, -1.6];

% Thrust output for forward flight Cheeseman1955
% FWratio =  [ 0    0.4  0.8  1.2  1.6  2.0  ];
% FWTratio = [ 1.31 1.26 1.16 1.08 1.04 1.005];  

% Thrust output for forward flight Cheeseman1955
FWratio =  [ 0   0.2 0.4  0.8  1.2  1.6  2.0  ]; % forward velocity
FWTratio = [ 1.31 1.26 1.16 1.08 1.04 1.005]; % thrust genearation
    
n = 10;
a = 1;
bldLoadNorm = 5;
Ct_bldLoad = 0.0589;

FWrange = v_free/vi;                         


for i = 1:length(z)
  
   for j = 1:length(u)
    % Coefficent of thrust in forward flight includeing GE - From Literature
    ctGEFW(i, j) = 1.00 / (1.00 - ( ( R / (4.00*z(i)) )^2.00 / (1.00 + ( u(j) / lam )^2.00)) );
%      ctGEFW(i, j) = 1.00 / (1.00 - ( (1/16)*(R /z(i))^2.00 * 1 / (1.00 + ( v_free(j) / vi )^2.00)) );

    % Thrust output for forward flight Cheeseman1955
    TgTinf(i,j) = 1 + (0.25*n*a*sqrt(0.12)) / ( sqrt(Ct_bldLoad)) * [ 1 / ( 16*(0.5)^2*(1 + (u(j))^2.00)) ] ;
    
   end
    
    
end


figure
% hold on;
x = HR;
y = u;
z = ctGEFW;

surf(x, y, z')
% surf(x, y, z2)
xlabel('HR'); ylabel('u'); zlabel('Ct');
xlim([0 2.5]); ylim([0 1.6]); zlim([1 1.3]);
% hold off; 
grid on; 
title('Thrust Ratio with GE and Forward flight')
colorbar;

figure
x = FWrange;
y = ctGEFW;
plot(x, y)
xlabel('V/vi'); ylabel('TgTinf');
grid;
% ylim([1.0, 1.3]); xlim([0, 2]);
title('FW flight?')
colorbar;

