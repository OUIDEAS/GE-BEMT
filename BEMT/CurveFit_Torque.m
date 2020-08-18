clear
clc
close all

x_new = 0.1:0.2:3.0;




load Torque_Curve.mat
load ImageGuess.mat

% X = [X, 5.5, 8.0]
% Q = [Q, 1, 1]


for i=1:length(Q)
    if X(i) > 1.1
        Q(i) = 1.000;
    end
%     if X(i) >=2 && X(i) > 1.5
%         Q(i) = 0.9998;
%     end
    
end

XX = X(1:12);
QQ = Q(1:12);
Vq = interp1(X, Q, X, "pchip");

x = X;
y = Q;
nonLin = @(b,x) b(1) + b(2)*(exp(-x./b(3))) + b(4)*b(6)./(b(5)*(((x).^0.909))) - b(7)./(x).^(b(8)); % ;              % Generalised exponential half life
% nonLin = @(b,x) b(1) + b(2)./(exp(x./b(3)));
NRCF = @(b) norm(y - nonLin(b,x));                       % Residual Norm Cost Function
B0 = [1; 1; 1; 1; 1; 1; 1;1];
B = fminsearch(NRCF, B0);   % Estimate Parameters

figure(1)
plot(x, y, 'pg')
hold on
plot(x, nonLin(B,x), '-r')
scatter(x, nonLin(B,x))
hold off
grid
ylim([0.94, 1.02])
%xlim([0 1.5])
text(0.7, 0.52, sprintf('y = %.4f %+.4f/(x %+.4f)', B))
b = B


HR=[0.1:0.1:3];

for l=1:length(HR)
cqio(l) = b(1) + b(2)*(exp(-HR(l)./b(3))) + b(4)*b(6)./(b(5)*(((HR(l)).^0.909))) - b(7)./(HR(l)).^(b(8));

end

figure
hold on
plot(X, Q, "*", "MarkerSize", 10, "Color", "k")
plot(HR, cqio,"-","LineWidth",2,"Color",'r')
hold off
grid on;
title("Curve fit")
xlabel("HR")
ylabel("Cq")
legend("Rc Benchmark", "Curve Fit");
% 






