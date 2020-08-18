close all;
clear;
clc;

load Torque_Curve.mat
load normal.mat

Nft = Ft;
Nptot = Ptoteff;
Nrpm = RPMreq;

load myFit.mat

myFt = Ft;
myPtot = Ptoteff;
myRpm = RPMreq;


diffFt = Nft-myFt;
diffPtot = Nptot-myPtot;
diffRpm = Nrpm-myRpm ;

FtperDiff = abs(diffFt./((Nft+myFt)./2)*100);
PtotperDiff = abs(diffPtot./((Nptot+myPtot)/2)*100);
RPMperDiff = abs(diffRpm./((Nrpm+myRpm)/2)*100);


figure
hold on;
plot(HR, Nft,"*", "MarkerSize", 10, "Color", "k")
plot(HR, myFt,"-","LineWidth",2,"Color",'b')
hold off;
grid on;
title("Flight Time vs HR")
xlabel("HR")
ylabel("Flight Time (mins)")
legend('Orrignal Model', 'With Torque')


figure
hold on;
plot(HR, Nptot,"*", "MarkerSize", 10, "Color", "k")
plot(HR, myPtot,"-","LineWidth",2,"Color",'b')
hold off;
grid on;
title("Power vs HR")
xlabel("HR")
ylabel("Power (Watts)")
legend("Orrignal Model", "With Torque")


figure
hold on;
plot(HR, Nrpm,"*", "MarkerSize", 10, "Color", "k")
plot(HR, myRpm,"-","LineWidth",2,"Color",'b')
hold off;
grid on;
title("RPMvs HR")
xlabel("HR")
ylabel("RPM(-)")
legend("Orrignal Model", "With Torque")

% -------

% figure
% plot(HR, FtperDiff)
% grid on;
% title("Flight Time Percent Difference")
% xlabel("HR")
% ylabel("%")
% 
% figure
% plot(HR, PtotperDiff)
% grid on;
% title("Power Percent Difference")
% xlabel("HR")
% ylabel("%")
% 
% figure
% plot(HR, RPMperDiff)
% grid on;
% title("RPM Percent Difference")
% xlabel("HR")
% ylabel("%")

% c1=0.2685;
% c2=-1.312;
% c3=0.9949;
% c4=0.001212;
% 
% %Height Ratio
HR=[0.1:0.1:3];
%   p1 = 0.040746;
%   p2 = -0.44068;
%   p3 = 1.9212;
%   p4 = -4.3359;
%   p5 = 5.3981;
%   p6 = -3.6263;
%   p7 = 1.1677;
%   p8 = 0.86936;
% 
% 
% 

% power curve fit coeeficients  
  q1 = 0.9925;
  q2 = 0.0026;
  
% % Hyperbolic curve fit
% b1 = 0.9970;
% b2 = 0.0005;
% b3 = -0.1135;  

% Exponetial half-life curve fit
b1 = 0.9983;
b2 = -2.0994;
b3 = 0.0181;  


for l=1:length(HR)
%ctio2(l)=c1*exp(c2*HR(l))+c3*exp(c4*HR(l));
%cqio2(l) = p1*HR(l)^7 + p2*HR(l)^6 + p3*HR(l)^5 + p4*HR(l)^4 + p5*HR(l)^3 + p6*HR(l)^2 + p7*HR(l) + p8; 
%cqio2(l) = q1*exp(q2*HR(l));
cqio2(l) = b1 + b2./(2.^(HR(l)/b3));

end


% figure
% hold on;
% plot(x, TorqueC)
% plot(HR, cqio2)
% grid on
% title("Exponential Half-life Curve fit to RC Benchmark")
% xlabel("HR")
% ylabel("Cq")
% legend("Rc Benchmark", "Curve Fit")

 
% x = x;
% y = TorqueC;
% hyprb = @(b,x) b(1) + b(2)./(x + b(3));                 % Generalised Hyperbola
% NRCF = @(b) norm(y - hyprb(b,x));                       % Residual Norm Cost Function
% B0 = [1; 1; 1];
% B = fminsearch(NRCF, B0);% Estimate Parameters
% 
% figure(1)
% plot(x, y, 'pg')
% hold on
% plot(x, hyprb(B,x), '-r')
% scatter(x, hyprb(B,x))
% hold off
% grid
% text(0.7, 0.52, sprintf('y = %.4f %+.4f/(x %+.4f)', B))

% x = x;
% y = TorqueC;
% nonLin = @(b,x) b(1) + b(2)./(2.^(x/b(3)));              % Generalised exponential half life
% NRCF = @(b) norm(y - nonLin(b,x));                       % Residual Norm Cost Function
% B0 = [1; 1; 1];
% B = fminsearch(NRCF, B0);   % Estimate Parameters
% 
% figure(1)
% plot(x, y, 'pg')
% hold on
% plot(x, nonLin(B,x), '-r')
% hold off
% grid
% text(0.7, 0.52, sprintf('y = %.4f %+.4f/(x %+.4f)', B))





% i = 0; 
% for i = 1:61
%     
%     Avg(i) = mean(ct(i,:));
%     
% end
% 
% figure
% plot(HR(1,:), Avg(1,:))
% title('Plot of Average Ct at each HR for each element')
% xlim([0,3])
% grid on;
% xlabel('HR')
% ylabel('Ct average')
% 
% % i = 0; 
% % for i = 1:61
% %     
% %     Avg3(i) = mean(ctio(i,:));
% %     
% % end
% 
% figure
% plot(HR(1,:), ctio(:,1))
% title('Ctio vs HR Using b and Omega')
% xlim([0,3])
% grid on;
% xlabel('HR')
% ylabel('Ctio')
% 
% figure
% plot(As)
% title('values of a Using b and Omega')
% grid on;
% xlabel('entry')
% ylabel('a')
% 
% figure
% plot(HR(1,:),As)
% title('values of a vs HR Using b and Omega')
% grid on;
% xlabel('HR')
% ylabel('a')
% 
% figure
% plot(Bs)
% title('values of b Using b and Omega')
% grid on;
% xlabel('entry')
% ylabel('b')
% 
% figure
% plot(HR(1,:),Bs)
% title('values of b vs HR Using b and Omega')
% grid on;
% xlabel('HR')
% ylabel('b')
% 
