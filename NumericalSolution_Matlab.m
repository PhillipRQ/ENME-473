clc; 
clear; 
close all;

[R1,R2,R3,R4,theta1,theta2] = deal(5,3,7,6,0,60*pi/180);

%% Newton-Raphson method
theta3 = input("Input theta3: ");
theta4 = input("Input theta4: ");
% x = [5*pi/180;200*pi/180]; %if initial guess theta3=5, theta4=200
x = [theta3*pi/180; theta4*pi/180]; %if initial guess theta3=10, theta4=100

fprintf('Initial Guess is: theta3 = %.4f degrees, theta4 = %.4f degrees',x(1),x(2));

theta3 = x(1); theta4 = x(2);
i = 1;
j=1;
while(true)
    f = [R2*cos(theta2)+R3*cos(theta3)-R4*cos(theta4)-5;
        R2*sin(theta2)+R3*sin(theta3)-R4*sin(theta4)];
    J = [-R3*sin(theta3) R4*sin(theta4); R3*cos(theta3) -R4*cos(theta4)];
    x = x - inv(J)*f;
    
    fprintf('\nIteration number: %d \n',i);
    fprintf('Newton-Raphson Guess: \n theta3=%.4f degrees, \n theta4=%.4f degrees ',x(1)*180/pi,x(2)*180/pi);

    if (abs(theta3*180/pi-x(1)*180/pi)<0.01 && abs(theta4*180/pi-x(2)*180/pi)<0.01)
        break;
    end;
    theta3 = x(1); theta4 = x(2);

    i=i+1;
end;

fprintf('\nNewton Raphson solutions:\n theta3 = %.4f degrees, theta4 = %.4f degrees', theta3*180/pi, theta4*180/pi)



