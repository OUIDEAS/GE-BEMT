%Blade Element Momentum/ Vortex Panel Method Combined Theory
%Author:
%   Aaron Michael Harrington
%   Aerospace Engineer
%   University of Maryland
%Revision:
%   2.3
%Last Updated:
%   February 24, 2010
%   3:23 pm
%
%Statement of Purpose:
%   The idea to incorporate the vortex panel method with blade element
%   momentum theory came from a necessity to use any airfoil without the
%   knowledge of its lift and drag properties. The vortex panel method will
%   calculate the pressure distribution around the arbitrary airfoil based
%   upon an assembly of vortecies of appropriate strength between
%   coordinates specified in the airfoil coordinate file. The vortex panel
%   method has been outlined in many books and is an appropriate method for
%   first order approximations ONLY as it does not accurately model flow
%   separation, boundary layers, etc. The panel method used in this program
%   was modified from its original form and was used compliments of L. Sankar
%   codded April 1997. 
%Validation:
%   Validation of this combined method was done using NASA TM 3218,
%   "Full-Scale-Tunnel Investigation of the Static-Thrust Performance of a Coaxial Helicopter Rotor"
%   By Robert D. Harrington
%   The program seems to overestimate both lift and drag and therefore
%   overestimates Thrust AND Power estimates

clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Global Variable Deffinitions%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamm = 1.4;                     %Assumed Air Constant
    R = 286;                        %Known Constant
    Cla = 5.73;                     %Based on NASA TM by Robert D. Harrington
    atmviscosity = 1.7332.*10.^-5;  %Standard Viscosity of Air

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Variable Initialization    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Solidity = 0;
    Ttot = 0;
    Qtot = 0;
    P_itot = 0;
    P_ptot = 0;
    Ptot = 0;
    CT = 0;
    CP = 0;
    Convergence = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Beginning of the Main Program%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Determine the Calculation Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    airfoil =   input ('Enter the Airfoil input file to use ');
    tip     =   input ('Include Prandtl Tip Loss:   (1 = Ignore,    2 = Include) ');
    twist   =   input ('Twist:  (1 = Linear or None, 2 = Ideal) ');
    taper   =   input ('Taper:  (1 = Linear or Constant, 2 = Ideal ');

%%%%%%% Determine the Blade Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nelements       = input ('Number of Elements: ');
    bladeNb        = input ('Number of Blades: ');
    bladeradius    = input ('Blade Radius (meters): ');
    if taper == 1
        bladerootchord = input ('Root Chord (meters): ');
    end
    bladetipchord  = input ('Tip Chord (meters): ');
    bladecutout    = input ('Cutout (meters): ');
    bladecollectivepitch = input ('Collective Pitch (degrees): ');
    if twist == 1
        bladetheta75 = input ('Theta .75 Span (degrees): ');
        bladetwist = input('Twist (degrees): ');
        bladetwist = bladetwist.*pi./180;
    elseif twist == 2
        bladethetatip = input ('Theta Tip (degrees): ');
        bladethetatip = bladethetatip.*pi./180;
    end
    bladeomega = input('Rotational Rate (RPM): ');
    %Vclimb = input('Climb Velocity (m./sec): ')
    Vclimb = 0;
    %atm.altitude = input ('Altitude (meters): ')
    atmaltitude = 0;
    
%%%%%%% Beginning of Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    atmtemperature = (15.04-0.00649.*atmaltitude)+273.1;        %Temperature from Standard Atmosphere Equation
    atmpressure = 101.29.*(atmtemperature./288.8).^5.256;       %Pressure from Standard Atmosphere
    atmdensity = atmpressure./(.2869.*atmtemperature);          %Density from Standard Atmosphere
    atmVsound = sqrt(gamm.*atmtemperature.*R);                  %Speed of Sound
    dy = (bladeradius-bladecutout)./Nelements;                  
    dr = (1./Nelements);
    bladecollectivepitch = bladecollectivepitch.*pi./180;       %Conversion from Degrees to Radians
    bladDiskArea = (pi.*bladeradius.^2);                        %Calculation of Disk Area
    bladeomega = (bladeomega./60.*2.*pi);                       %RPM -> RPS
    bladeVtip = bladeradius.*bladeomega;                        %Rotor Tip Speed
    lambdac = Vclimb./bladeVtip;                                %Non-dimensional Climb Velocity
    

%%%%%%% Calculation of Element Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:Nelements
        F = 1;
        old = 0;
        last = 0;
        y_location(i)   = bladecutout+dy.*i;
        r(i)            = y_location(i)./bladeradius;
        
%%%%%%% Calculations Relating to Blade Twist %%%%%%%%%%%%%%%%%%%%%%%%%
            if twist == 1
                twist_dist(i)   = bladetheta75+(r(i)-.75).*bladetwist;  %This Calculation is for Linear Twist
            elseif twist == 2
                twist_dist(i)   = bladethetatip./r(i);                  %This Calculation is for Ideal Twist
            end
        
        local_angle(i)  = bladecollectivepitch+twist_dist(i);
        
%%%%%%% Calculations Relating to Blade Taper %%%%%%%%%%%%%%%%%%%%%%%%%
            if taper == 1
                element_chord(i) = bladerootchord-(bladerootchord-bladetipchord).*r(i); %This Calculation is for Linear Taper
            elseif taper == 2
                element_chord(i) = bladetipchord./r(i);                                 %This Calculation is for Ideal Taper        
            end
%%%%%%% Calculation of Solidity Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%

        local_solidity(i) = (bladeNb.*element_chord(i))./(pi.*bladeradius);                       %Find the Local elemental solidity of the blade
        Solidity = Solidity + local_solidity(i);                                        %Compute the running total to calculate rotor solidity
        
%%%%%%% Calculation of the Inflow Parameters %%%%%%%%%%%%%%%%%%%%%%%%%        
        local_inflow_ratio(i) = sqrt(((local_solidity(i).*Cla)./(16.*F)-lambdac./2).^2+(local_solidity(i).*Cla*local_angle(i).*r(i))./(8.*F))-((local_solidity(i).*Cla)./(16.*F)-lambdac./2);
        V_inflow(i) = local_inflow_ratio(i).*bladeVtip-Vclimb;                          %Inflow Velocity Calculation
        V(i) = sqrt((V_inflow(i)+Vclimb).^2+(y_location(i).*bladeomega).^2);            %Total Velocity Calculation
        inducedangle = local_inflow_ratio(i)./r(i);
        inflow_angle(i) = atan((Vclimb+V_inflow(i))./(bladeomega.*y_location(i)));       %Inflow Angle Calculation
        local_mach_number(i) = r(i).*(bladeVtip./atmVsound);                            %Local Mach Number Calculation
        Re(i) = (atmdensity.*V(i).*element_chord(i))./(atmviscosity);                   %Local Reynold's Number Calculation
        tip;
%%%%%%% Inclusion of Tip Losses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if tip == 2
             if r(i) == 1
                 F = 0.0000000001;
                 local_inflow_ratio(i) = sqrt(((local_solidity(i).*Cla)./(16.*F)-lambdac./2).^2+(local_solidity(i).*Cla*local_angle(i).*r(i))./(8.*F))-((local_solidity(i).*Cla)./(16.*F)-lambdac./2);
                 V_inflow(i) = local_inflow_ratio(i).*bladeVtip-Vclimb;                     %Inflow Velocity
                 V(i) = sqrt((V_inflow(i)+Vclimb).^2+(y_location(i).*bladeomega).^2);       %Total Velocity
                 inflow_angle(i) = atan((Vclimb+V_inflow(i))./(bladeomega.*y_location(i)));
             else
               while(abs(local_inflow_ratio(i)-old)/local_inflow_ratio(i)>0.005)
                    old = local_inflow_ratio(i);
                    f_tip = (bladeNb./2).*((1-r(i))./(r(i).*inflow_angle(i)));
                    f_root = (bladeNb./2).*((r(i))./((1-r(i))*inflow_angle(i)));
                    f = f_root.*f_tip;
                    F = (2./pi).*acos(exp(-1.*(f_tip)));
                    local_inflow_ratio(i) = sqrt(((local_solidity(i).*Cla)./(16.*F)-lambdac./2).^2+(local_solidity(i).*Cla*local_angle(i).*r(i))./(8.*F))-((local_solidity(i).*Cla)./(16.*F)-lambdac./2);
                    V_inflow(i) = local_inflow_ratio(i).*bladeVtip-Vclimb;
                    V(i) = sqrt((V_inflow(i)+Vclimb).^2+(y_location(i).*bladeomega).^2);
                    inflow_angle(i) = atan((Vclimb+V_inflow(i))./(bladeomega.*y_location(i)));
                end
             end
        end
        
        local_AOA(i) = local_angle(i)-inflow_angle(i);
        
%%%%%%% Calculation of Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Vortex Panel Method   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cl(i), cd(i), x, y] = VPM(airfoil,(local_AOA(i)*180/pi));  %The input variables to VPM are the airfoil data file and the angle of attack in degrees
        x = x*element_chord(i);                                     %This will size the airfoil x coordinates based on the elemental chord length at position r(i)
        y = y*element_chord(i);                                     %This will size the airfoil y coordinates based on the elemental chord length at position r(i)

        figure(1);
        hold on;
        lengthofx = length(x);
        z = ones(length(x))*y_location(i);
        plot3(x,y,z);
        view(70,-10);
        axis([-bladeradius bladeradius -bladeradius bladeradius 0 bladeradius]);
        
%         %This is strictly a linear estimate 
%         cl(i) = Cla.*local_AOA(i);
%         cd(i) = 0.0087-0.02167.*local_AOA(i)+.4.*(local_AOA(i)).^2;
%         
%%%%%%% Differentials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dL1(i) = 0.5.*atmdensity.*(V(i)).^2.*element_chord(i).*cl(i).*dy;
        dD1(i) = 0.5.*atmdensity.*(V(i)).^2.*element_chord(i).*cd(i).*dy;
        
        dT(i) = bladeNb.*(dL1(i).*cos(inflow_angle(i))-dD1(i).*sin(inflow_angle(i)));
        Ttot = Ttot + dT(i);
        
        dQ(i) = bladeNb.*(dD1(i).*cos(inflow_angle(i)+dL1(i).*sin(inflow_angle(i))));
        Qtot = Qtot + dQ(i);
        
        dP_i(i) = bladeNb.*(dL1(i).*sin(inflow_angle(i)).*y_location(i).*bladeomega);
        P_itot = P_itot + dP_i(i);
        dP_p(i) = bladeNb.*(dD1(i).*cos(inflow_angle(i)).*y_location(i).*bladeomega);
        P_ptot = P_ptot + dP_p(i);
        dP(i) = dP_i(i) + dP_p(i);
        Ptot = Ptot + dP(i);
        
        dCt(i) = dT(i)./(atmdensity.*pi.*(bladeradius).^2.*(bladeVtip).^2);
        dCp(i) = dP(i)./(atmdensity.*pi.*(bladeradius).^2.*(bladeVtip).^3);
        CT = CT + dCt(i);
        CP = CP + dCp(i);
        
        Disk_Loading = CT.*atmdensity.*(bladeVtip).^2;
        Ct_over_sigma = CT./Solidity;
    end

figure(2);
hold on;
plot(r,dL1,'b-')
plot(r,dD1,'r-')
xlabel('r')
title('Elemental Lift and Drag Along Blade Radius')
legend('dL','dD')

figure(3);
plot(r,local_inflow_ratio)
xlabel('r')
ylabel('{\lambda}')
title('Inflow Along Blade Radius')

Total_Thrust = Ttot*.22481
Total_Torque = Qtot*.73756
Total_Power = Ptot*.001341
CT
CP
Disk_Loading
Ct_over_sigma