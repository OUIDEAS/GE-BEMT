function [spline, spline_time] = fitCubicSplines(theta, time, dt)
    
    if ~isvector(theta) || ~isvector(time) % Check if inputs are valid vectors
        
        fprintf('Inputs must be vectors!\n');
        
    elseif length(theta) ~= length(time) % Check if inputs are the same length
        
        fprintf('Inputs must be same length!\n');
    
    else % If everything is fine, continue
        
        theta_dot_calc = zeros(1,length(theta)+1);
        for i = 1:length(theta)-1
            theta_dot_calc(i+1) = (theta(i+1)-theta(i))/(time(i+1)-time(i));
        end
        
        theta_dot = zeros(1,length(theta));
        for i = 1:length(theta_dot_calc)-1
            if(sign(theta_dot_calc(i))~=sign(theta_dot_calc(i+1)))
                theta_dot(i) = 0;
            else
                theta_dot(i) = (theta_dot_calc(i+1) + theta_dot_calc(i))/2;
            end
        end
        
        spline = [];
        spline_time = [];
        for i = 1:length(theta)-1
            
            theta_curr = theta(i);
            theta_next = theta(i+1);
            theta_curr_dot = theta_dot(i);
            theta_next_dot = theta_dot(i+1);
            time_curr = time(i);
            time_next = time(i+1);
                        
            [s_curr, t_curr] = cubicSplines(theta_curr, theta_next, theta_curr_dot, theta_next_dot, time_curr, time_next, dt);
            
            spline = [spline, s_curr];
            spline_time = [spline_time, t_curr]; 
        end
        
    end
    
    figure(1);
    plot(time, theta,'b');
    hold on;
    hl = plot(spline_time,spline,'r');
    hold off;
    grid on;
    set(hl,'LineWidth',1.5);
    title('Fitted Spline and Given Angles vs. Time','FontSize',11,'FontWeight','bold');
    ylabel('Theta');
    xlabel('Time');
    legend('Given Theta', 'Fitted Spline')
    
    figure(2)
    subplot(3,2,1)
    plot(time, theta,'b');
    hold on;
    hl = plot(spline_time,spline,'r');
    hold off;
    xlim([-0.1 0.1])
    ylim([9.9 10.1])
    grid on;
    set(hl,'LineWidth',1.5);
    title('Theta-1 Zoom In','FontSize',11,'FontWeight','bold');
    ylabel('Theta');
    xlabel('Time');
    
    subplot(3,2,2)
    plot(time, theta,'b');
    hold on;
    hl = plot(spline_time,spline,'r');
    hold off;
    xlim([1.9 2.1])
    ylim([34.9 35.1])
    grid on;
    set(hl,'LineWidth',1.5);
    title('Theta-2 Zoom In','FontSize',11,'FontWeight','bold');
    ylabel('Theta');
    xlabel('Time');

    subplot(3,2,[3,4])
    plot(time, theta,'b');
    hold on;
    hl = plot(spline_time,spline,'r');
    hold off;
    xlim([2.9 3.1])
    ylim([24.9 25.1])
    grid on;
    set(hl,'LineWidth',1.5);
    title('Theta-3 Zoom In','FontSize',11,'FontWeight','bold');
    ylabel('Theta');
    xlabel('Time');
    
    subplot(3,2,5)
    plot(time, theta,'b');
    hold on;
    hl = plot(spline_time,spline,'r');
    hold off;
    xlim([5.9 6.1])
    ylim([29.9 30.1])
    grid on;
    set(hl,'LineWidth',1.5);
    title('Theta-4 Zoom In','FontSize',11,'FontWeight','bold');
    ylabel('Theta');
    xlabel('Time');
        
    subplot(3,2,6)
    plot(time, theta,'b');
    hold on;
    hl = plot(spline_time,spline,'r');
    hold off;
    xlim([7.9 8.1])
    ylim([59.9 60.1])
    grid on;
    set(hl,'LineWidth',1.5);
    title('Theta-5 Zoom In','FontSize',11,'FontWeight','bold');
    ylabel('Theta');
    xlabel('Time');
end

function [spline, time] = cubicSplines(theta_curr, theta_next, theta_curr_dot, theta_next_dot, time_0, time_final, dt)
    dur = time_final - time_0;
    a0 = theta_curr;
    a1 = theta_curr_dot;
    a2 = ((3/(dur^2))*(theta_next-theta_curr))-((2/dur)*theta_curr_dot)-((1/dur)*theta_next_dot);
    a3 = (-(2/(dur^3))*(theta_next-theta_curr))+((1/dur^2)*(theta_next_dot+theta_curr_dot));
    
    spline = zeros(1,length(time_0:dt:time_final));
    time = time_0:dt:time_final;
    t_interval = 0:dt:dur;
    for i = 1:length(t_interval)
        t = t_interval(i);
        spline(i) = a3*(t^3)+a2*(t^2)+a1*t+a0;
    end
end
