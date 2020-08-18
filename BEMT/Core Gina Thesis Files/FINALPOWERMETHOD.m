% Gina M. Eberhart
% GE BEMT Thrust Prediction
clear all;
clc;
close all;
%% Initial Conditions and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import Tabulated Airfoil Data
clarky = Airfoilimport('clarkyre100krenncrit9.txt');
%clarky = Airfoilimport('Naca44121ookren.txt');

% Load Twist and Taper Provile for Propeller
load TwistandTaperData10x4.mat;
%load TwistandTaperData.mat;

%Pitch (rads)
th=deg2rad(theta);

% Radial Chord (m)
chord=cinterp;

% Rotor Radius (m)
R=.123;
%R=.139

% Area of Rotor Disk
A=pi*R^2;

% Number of Propellers
NP=4;

% Acceleration Due to Gravity (m/s^2)
g=9.81;

% Payload Weight (N)
WP=0*g;

% Total Weight of Quadcopter (N)
WT= 1.435*g+WP;

% Thrust Required to Hover From Single Propeller (N)
ThrustSP=WT/NP;

% Density of Air (kg/m^3)
rho=1.225;

% Number of Elements to Analyze
N=20;

% Length of Element (m)
dr=R/N;

% rstep=(xt-xs)/10
r1=[dr:dr:R];


%Constants for lambdaIO
% 11x7 thrust curve
% c1=0.0974;
% c2=-1.935;
% c3=1.005;
% c4=-0.00049;
% 10x4.7 thrust curve
% c1=0.2685;
% c2=-1.312;
% c3=0.9949;
% c4=0.001212;

c1=0.2685;
c2=-1.312;
c3=0.9949;
c4=0.001212;

% exponential curve fit coeeficients  
  q1 = 0.9925;
  q2 = 0.0026;

% Exponetial half-life curve fit
p1 = 0.9983;
p2 = -2.0994;
p3 = 0.0181;

% My fit
% B = [ 0.9973 0.0808 0.2801 1.2128 1.8216 1.6591 1.1032 0.9139];
 
% Ct coefficents
AA = [0.3101, -3.685, 1.001, -0.002];
CC = [1, -1.6];
% Cq coefficients
BB = [ 0.0364, -2.825, 0.2600, -14.43];


%% Thrust Prediction
% Tip Speed Ratio
lambda=41.35;

%Height Ratio
HR=[0.5:0.1:6];

%Rotor RPM
RPM=[4000:250:10000];

fprintf("Begin Final Power Method... \n");

for l=1:length(HR)
    txt = ['Running...', num2str(l),'/',num2str(length(HR))];
   fprintf(txt);
    
    for k=1:length(RPM)
        fprintf('.');

        omega(k)=(RPM(k)/60)*2*pi; % b?
        for i=1
            

            % Free stream velocity guess
            V(i)=(omega(k)*R)/lambda(i);
            
            thrust=0.0;
            torque=0.0;
            for j=1:length(r1)
   

                %Going Along the Radius of the blade rad=radius
                r=r1(j);
                
                %sigma Nb*c/2piR
                sigma=(2*chord(j))/(pi*r);
                
                %initial guess for a
                a=0.1;
                
                %initial guess for b
                b=0.01;
                
                %setting a value for while loop conditions
                finished=0;
                sum=1;
                ii = 1;
%                 V0   =0; V2   = 0; phi = 0; Vlocal = 0; 
%                 cl   =0; cd   = 0; ctio   = 0; 
%                 ct   =0; cqio = 0; cq     = 0; 
%                 DtDr =0; tem1 = 0; tem2   = 0; 
%                 anew =0; bnew = 0;
                % record = [];
                
                while (finished==0),
                    if j ==4 && ii == 154
                        uuuu = 0;
                    end
                    
                    
                    %solving for V0 and V2 using free stream velocity and induction factor
                    V0=V(i)*(1+a);
                    V2=omega(k)*r*(1-b);
                    %solving for phi using velocity components
                    phi=atan2(V0,V2);
                    
                    %solving for angle of attack
                    alpha=th(j)-phi;
                    
                    %Finding CL and CD based on angle of attack from Tabulated Data
                    if rad2deg(alpha)>=min(clarky.alpha) && rad2deg(alpha)<max(clarky.alpha)
                        cl=interp1(deg2rad(clarky.alpha),clarky.CL,alpha);
                        cd=interp1(deg2rad(clarky.alpha),clarky.CD,alpha);
                    else
                        cl=0;
                        cd=1;
                    end
                    
                    % finding local velocity Sqrt(v0^2 +v2^2)
                    Vlocal=sqrt(V0^2+V2^2);
                    
                    Inflow1 = Vlocal/(omega(k)*r);
                    Inflow2 = Vlocal/(omega(k)*R);
                    
                    %Coefficients of thrust and torque
%                   ctio=c1*exp(c2*HR(l))+c3*exp(c4*HR(l));
%                     if HR(l) > 3
%                         ctio = 1.0;
%                         cqio = 1.0;
%                     else
%                         
%                     end

                    ctio = AA(1)*exp(AA(2)*HR(l)) + AA(3)*exp(AA(4)*HR(l))+ CC(1)*exp(CC(2)*HR(l));
                    ct=(cl*cos(phi)-cd*sin(phi))*ctio;
                    
                    cqio = 1 + BB(1)*(exp(BB(2)*HR(l))) + BB(3)*(exp(BB(4)*HR(l)));                   
                    cq=(cd*cos(phi)+cl*sin(phi))*cqio;
                    
                    % Finding Differential Thrust
                    DtDr=0.5*rho*Vlocal^2*2*chord(j)*ct;
                    
                    % Finding Differential Torque
                    DqDr=0.5*rho*Vlocal^2*2*chord(j)*r*cq;
                    
                    % BEMT thrust equation solving for
                    tem1=DtDr/(4*pi*r*rho*V(i)^2*(1+a));
                    tem2=DqDr/(4*pi*r^3*rho*V(i)*(1+a)*omega(k)); % there is no b?
                    anew=0.5*(a+tem1);
                    bnew=0.5*(b+tem2);
                    if (abs(anew-a)<1.0e-5),
                        if (abs(bnew-b)<1.0e-5),
                            finished=1;
                        end;
                    end;
                    a=anew;
                    b=bnew;
                    sum=sum+1;
                    if (sum>500),
                        finished=1;
                    end;
                    ii = ii+1;

                    
                    
                end;
                phieach(j)=rad2deg(phi);
                cleach(j)=cl;
                alphadegs(j)=rad2deg(alpha);
                a2(j)=a;
                b2(j)=b;
                thrust=thrust+DtDr*dr;
                torque=torque+DqDr*dr;
            end;
            
            thrustis(k,i)=thrust;
            torqueis(k,i)=torque;


            
        end
  
    end
    fprintf("Done \n");
    record(l, :) = [a, b, V0, V2, phi, Vlocal, cl, cd, HR(k), ctio, ct, cqio, cq, DtDr, tem1, tem2, anew, bnew];
   
    % Rotational Speed to produce 1/4 of T Required to Hover (for quadcopter)
    Omegareq(l)=interp1(thrustis,omega,ThrustSP);
    %
    % Torque OGE for Single Prop
    TorqueOGE(l)=interp1(omega,torqueis,Omegareq(l));
    
    %Coefficeint of Torque for Single Prop
    CqOGE(l)=TorqueOGE(l)/(rho*A*Omegareq(l)^2*R^3);
    
    % Equating Coefficient of Power and Torque
    CpOGE(l)=CqOGE(l);
    
    % Power Required OGE for Single Prop
    P=rho*A*(Omegareq(l)*R)^3*CpOGE(l);
    %Motor Efficiency (70%)
    MEff=0.70;
    
    %Total Power Required
    Ptot=NP*P;
    Ptoteff(l)=Ptot+Ptot*MEff;
    
    % Coefficient of Thrust
    Ct(l)=ThrustSP/(rho*(Omegareq(l)*R)^2*A);
    
    % RPM Required
    RPMreq(l)=(Omegareq(l)*60)/(2*pi);
    
    % Input voltage from Battery
    VoltIn=11.1;
    
    %Required Current
    Ireq=Ptoteff(l)/VoltIn;
    
    %Usable Battery Capacity (mAh)
    C=5100;
    
    % Flight Time Calculation (mins)
    Ft(l)=((C/1000)./Ireq)*60;
end

    fprintf("Complete... \n");

 %% Figures
figure
plot(HR,Ft,'o')
xlabel('HR (z/R)')
ylabel('Flight Time (Mins)')
title('Flight Time vs. HR')
grid on

figure
plot(HR,Ptoteff,'o')
xlabel('HR (z/R)')
ylabel('Power IGE (Watts)')
title('Power vs. HR')
grid on

figure
plot(HR,RPMreq,'o')
xlabel('HR (z/R)')
ylabel('RPM(-)')
title('RPM vs. HR')
grid on

