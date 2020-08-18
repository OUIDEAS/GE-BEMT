

clear; close all; clc;

t = 0;
z = 0; z_dot = 0; z_ddot = 0;
g = 9.81; m =1;theta = 0; phi = 0; dt = 0.05; kt = 1/4*m*g/(5000^2);

z_ge = 0; z_dot_ge = 0; z_ddot_ge = 0; R = 0.127*2; 
b1 = 0.3101; b2 = -3.685; b3 = 1.001; b4 = -0.002; 
% U1_ge = kt*b1*exp(b2)*(W)^2 + kt*b3*exp(b4)*(W)^2;

a=1;
u1 = 4*5100^2;  % z rpm command
u1_ge =4*5100^2;  % z rpm command


% PID position control
P=2; I=0; D=0.8; %1.1; I=0.035; D=2.8;
P_ge=2; I_ge=0; D_ge=0.8;
% Target Altitude
z_target = 0.5;

% Velocity control
P_climb=0.5; I_climb=0; D_climb=0;
P_climb_ge=1.75; I_climb_ge=0.00; D_climb_ge=0;


max_rpm = 4*10000^2;
max_climb = 2; 

for i = 2:100
% =============== Z desired
z_ddot_des(i) = -g/m + 1/m*cos(theta)*cos(phi)*kt*u1; 
% z_dot(i) = z_ddot(i-1)*dt + z_dot(i-1);
% z(i) = z_dot(i)*dt^2 + z(i-1);

% Inverted Dynamics
U1(i) = m/(kt*cos(theta)*cos(phi)) * (z_ddot_des(i) + g/m);
z_ddot(i) = -g/m + 1/m*cos(theta)*cos(phi)*kt*U1(i);

% if z(i-1) < 0 || z(i-1) == 0
%     z(i) = 0;
% end

z_dot(i) = z_ddot(i)*dt + z_dot(i-1);
z(i) = z_dot(i)*dt^2 + z(i-1);

% ============== Z GE desired
z_ddot_des_ge(i) = -g/m + 1/m*cos(theta)*cos(phi)*kt*u1_ge; 

% Inverted Dynamics
U1_ge(i) = m/( kt*cos(theta)*cos(phi)*(b1*exp(b2*z_ge(i-1)/R) + b3*exp(b4*z_ge(i-1)/R)) )*(z_ddot_des_ge(i) + g/m);
z_ddot_ge(i) = -g/m + 1/m*cos(theta)*cos(phi)*(b1*exp(b2*z_ge(i-1)/R) + b3*exp(b4*z_ge(i-1)/R))*kt*U1_ge(i); 

% if z_ge(i-1) < 0 || z_ge(i-1) == 0
%     z_ge(i) = 0;
% end

z_dot_ge(i) = z_ddot_ge(i)*dt + z_dot_ge(i-1);
z_ge(i) = z_dot_ge(i)*dt^2 + z_ge(i-1);

HR(i) = z_ge(i)/R;

%===== Height Controller ===================

error(i) = z_target - z(i);
error_ge(i) = z_target - z_ge(i);

error_sum(i) = sum(error)*dt;
error_sum_ge(i) = sum(error_ge)*dt;

error_rate(i) = (error(i)-error(i-1))/dt;
error_rate_ge(i) = (error_ge(i)-error_ge(i-1))/dt;

PID(i) = P*error(i) + I*error_sum(i) +D*error_rate(i);

PID_ge(i) = P_ge*error_ge(i) + I_ge*error_sum_ge(i) +D_ge*error_rate_ge(i);

if PID(i) < 0
    PID(i) =0;
elseif PID(i) > 1
    PID(i) = 1;
end

if PID_ge(i) < 0
    PID_ge(i) = 0;
elseif PID_ge(i) > 1
    PID_ge(i) = 1;
end

% Used for position control without a velocity controller
u1 = PID(i)*max_rpm;
u1_ge = PID_ge(i)*max_rpm;

% climb_Rate_target(i) = PID(i)*max_climb;
% climb_Rate_error(i) = climb_Rate_target(i) - z_dot(i);
% climb_error_sum(i)= sum(climb_Rate_error(i)); 
% climb_error_rate(i) = climb_Rate_error(i) - climb_Rate_error(i-1);
% 
% 
% climb_Rate_target_ge(i) = PID_ge(i)*max_climb;
% climb_Rate_error_ge(i) = climb_Rate_target_ge(i) - z_dot_ge(i);
% 
% 
% climb_PID(i) = P_climb*climb_Rate_error(i) + I_climb*climb_error_sum(i) +D_climb*climb_error_rate(i);
% 
% climb_PID_ge(i) = P_climb_ge*error_ge(i) + I_climb_ge*error_sum_ge(i) +D_climb_ge*error_rate_ge(i);
% 
% if climb_PID(i) < 0
%     climb_PID(i) = 0;
% elseif climb_PID(i) > 1
%     climb_PID(i) = 1;
% end
% 
% if climb_PID_ge(i) < 0
%     climb_PID_ge(i) = 0;
% elseif climb_PID_ge(i) > 1
%     climb_PID_ge(i) = 1;
% end
% 
% vel_cmd(i) = climb_PID(i)*max_rpm;
% 
% vel_cmd_ge(i) = climb_PID_ge(i)*max_rpm;
% 
% u1 = vel_cmd(i);
% u1_ge = vel_cmd_ge(i);



% if t(i) > 0.5
%     a = 2;
% elseif t(i) > 4
%     a = 3;
% elseif t(i) > 6
%     a =4; 
% elseif t(i) > 8
%     a =5;
% end



t(i) = t(i-1)+dt;
    
end




% x = [z_dot; z_ddot]; 
% A = [ 0 1; mu 0];
% B = [0 k; 0 0];
% U = [ 0; u1];

figure
hold on;
scatter(t, z_dot);
scatter(t, z_dot_ge, '*');

ylabel('Velocity')
xlabel('time')
legend('Normal DI', 'DI with GE')
grid on;
% xlim([0 30]);

figure
hold on;
scatter(t, z_ddot );
scatter(t, z_ddot_ge, '*');
% scatter(t, HR, '^');

ylabel('Acceleration')
xlabel('time')
legend('Normal DI', 'DI with GE')

grid on;

figure
hold on;
scatter(t, U1);
scatter(t, U1_ge, '*')

ylabel('Command (RPM)')
xlabel('time')
legend('Normal DI', 'DI with GE')

grid on;

figure
hold on;
scatter(t, z);
scatter(t, z_ge, '*');

ylabel('Position')
xlabel('time')
legend('Normal DI', 'DI with GE')

grid on;
% Less RPM is needed to get the same acceleration for flight inside GE vs
% flight outside GE
