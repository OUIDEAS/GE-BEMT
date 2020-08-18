clear all;
clc

%Chord Length
chord=0.10;

% Pitch
pitch=0.0;

% Rotor Diameter
dia=1.6;

% Rotor Radius
R=dia/2.0;

%Rotor RPM
RPM=2100.;

%pitch angle setting at tip
tip=25.0;

% Tip of Blade
xt=R;

% Pitch angle setting at 25% radius
hub=65.0;

% Hub of blade
xs=0.1*R
tonc=0.12*chord;

% Density of Air
rho=1.225;

% Converting RPM to rads/s
n=RPM/60.0;
omega=n*2.0*pi;


coef1=(tip-hub)/(xt-xs);
coef2=hub-coef1*xs;
rstep=(xt-xs)/20
r1=[xs:rstep:xt];
k=0;
eff0=0;

% Free stream velocity guess
V=60;

thrust=0.0;
torque=0.0;
for j=1:length(r1),
    
    %Going Along the Radius of the blade rad=radius
    rad=r1(j);
    
    %calculating theta (but I have theta already)
    theta=coef1*rad+coef2+pitch;
    t2(j)=theta;
    th=theta/180.0*pi;
    
    %sigma Nb*c/2piR
    sigma=2.0*chord/2.0/pi/rad;
    
    %initial guess for a
    a=0.1;
    
    %initial guess for b
    b=0.01;
    
    %setting a value for while loop conditions
    finished=0;
    sum=1;
    while (finished==0),
        
        %solving for V0 and V2 using free stream velocity and induction factor
        V0=V*(1+a);
        V2=omega*rad*(1-b);
        %solving for phi using velocity components
        phi=atan2(V0,V2)
        
        %solving for angle of attack
        alpha=th-phi;
        
        %Finding CL and CD based on angle of attack
        cl=6.2*alpha;
        cd=0.008-0.003*cl+0.01*cl*cl;
        
        % finding local velocity Sqrt(v0^2 +v2^2)
        Vlocal=sqrt(V0*V0+V2*V2);
        ct=(cl*cos(phi)-cd*sin(phi))
        
        % Finding Differential Thrust
        DtDr=0.5*rho*Vlocal*Vlocal*2.0*chord*(cl*cos(phi)-cd*sin(phi));
        
        % Finding Differential Torque
        DqDr=0.5*rho*Vlocal*Vlocal*2.0*chord*rad*(cd*cos(phi)+cl*sin(phi));
        
        % BEMT thrust equation solving for
        tem1=DtDr/(4.0*pi*rad*rho*V*V*(1+a));
        tem2=DqDr/(4.0*pi*rad*rad*rad*rho*V*(1+a)*omega);
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
    end;
    a2(j)=a;
    b2(j)=b;
    thrust=thrust+DtDr*rstep;
    torque=torque+DqDr*rstep;
end;
disp(thrust);
disp(torque);

% Coefficients of thrust and torque
t=thrust/(rho*n*n*dia*dia*dia*dia);
q=torque/(rho*n*n*dia*dia*dia*dia*dia);


% J=V/(n*dia);
% if (t<0),
%     eff=0.
% else
%     eff=t/q*J/(2.0*pi)
% end;