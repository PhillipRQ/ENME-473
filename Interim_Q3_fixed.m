% ENME 473 Project Deliverable 1 - Questions 3a & 3c
clc; clear;

%% Known Values
R1 = sqrt(237.2^2 + 70.6^2);
R2 = 272.4;
R23 = 236.9;
R14 = 279.0;
R3 = 489.9;
R4 = 539.8;
R36 = 342.1;
R46 = 179.3;
R5 = 595.5;
R6 = 378.9;
R7 = 198.2;
R8 = 254.4;
RA = 336.9;         % distance from joint 7/8 to Pin A along link 8
theta1 = deg2rad(180 - atand(70.6/237.2));
alpha = deg2rad(3.2);
beta = deg2rad(3.5);

% create the range of input angles
in = 0:120;
theta2_vals = deg2rad(in);

%% Symbolic equations and Jacobian
syms t23 t14 t46 t36 t8 t7 t2
f1 = R2*cos(t2) + R23*cos(t23) - R14*cos(t14) - R1*cos(theta1);
f2 = R2*sin(t2) + R23*sin(t23) - R14*sin(t14) - R1*sin(theta1);
f3 = R2*cos(t2) + R4*cos(t23) + R46*cos(t46) - R36*cos(t36) - R3*cos(t14+alpha) - R1*cos(theta1);
f4 = R2*sin(t2) + R4*sin(t23) + R46*sin(t46) - R36*sin(t36) - R3*sin(t14+alpha) - R1*sin(theta1);
f5 = R2*cos(t2) + R4*cos(t23) + R6*cos(t46) + R8*cos(t8) - R7*cos(t7) - R5*cos(t36+beta) - R3*cos(t14+alpha) - R1*cos(theta1);
f6 = R2*sin(t2) + R4*sin(t23) + R6*sin(t46) + R8*sin(t8) - R7*sin(t7) - R5*sin(t36+beta) - R3*sin(t14+alpha) - R1*sin(theta1);
f_sym = [f1; f2; f3; f4; f5; f6];
vars = [t23 t14 t46 t36 t8 t7];
J_sym = jacobian(f_sym, vars);
% Convert symbolic equations to numerical functions
f_fun = matlabFunction(f_sym, 'Vars', {t23, t14, t46, t36, t8, t7, t2});
J_fun = matlabFunction(J_sym, 'Vars', {t23, t14, t46, t36, t8, t7, t2});

%% Posture Analysis (Newton-Raphson)

% store initial guesses
x = deg2rad([140 60 150 150 210 150])';

% storage
angles = zeros(6, length(in));
Ax = zeros(1, length(in));
Ay = zeros(1, length(in));

% loop through input values
for k = 1:length(theta2_vals)
    theta2 = theta2_vals(k);

    while true
        F = f_fun(x(1), x(2), x(3), x(4), x(5), x(6), theta2);
        J = J_fun(x(1), x(2), x(3), x(4), x(5), x(6), theta2);
        dx = J\F;
        x_new = x - dx;
        x_new = atan2(sin(x_new), cos(x_new));

        if norm(F,inf) < 1e-6 && norm(x_new - x,inf) < 1e-6
            x = x_new;
            break;
        end
        x = x_new;
    end

    angles(:,k) = x;

    % compute Pin A position
    Ax(k) = R2*cos(theta2) + R4*cos(x(1)) + R6*cos(x(3)) + (R8 + RA)*cos(x(5));
    Ay(k) = R2*sin(theta2) + R4*sin(x(1)) + R6*sin(x(3)) + (R8 + RA)*sin(x(5));
end

%% Post-process Angles

% convert to degrees with unwrapping
angles = rad2deg(unwrap(angles, [], 2));
theta3_vals = angles(2,:) + rad2deg(alpha);
theta4_vals = angles(1,:);
theta5_vals = angles(4,:) + rad2deg(beta);
theta6_vals = angles(3,:);
theta7_vals = angles(6,:);
theta8_vals = angles(5,:);


%% Plot 1: Link Angles vs Input Angle (Q3a)
figure
plot(in, theta3_vals, 'LineWidth', 2)
hold on
plot(in, theta4_vals, 'LineWidth', 2)
plot(in, theta5_vals, 'LineWidth', 2)
plot(in, theta6_vals, 'LineWidth', 2)
plot(in, theta7_vals, 'LineWidth', 2)
plot(in, theta8_vals, 'LineWidth', 2)
xlabel('\theta_2 (deg)')
ylabel('Link Angle (deg)')
title('Link Angles vs Input Angle')
legend('\theta_3','\theta_4','\theta_5','\theta_6','\theta_7','\theta_8','Location','best')
grid on
box on
%exportgraphics(gcf, 'ENME_473_Project_3a.png', 'Resolution', 600);

%% Plot 2: Path Traced by Pin A (Q3c)
figure
plot(Ax, Ay, 'b-', 'LineWidth', 2)
hold on
% mark start and end points
plot(Ax(1),   Ay(1),   'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', '\theta_2 = 0°')
plot(Ax(end), Ay(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', '\theta_2 = 120°')
% mark every 30 degrees
for idx = [31 61 91]  % theta2 = 30, 60, 90
    plot(Ax(idx), Ay(idx), 'kd', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off')
    text(Ax(idx)+10, Ay(idx)+10, sprintf('\\theta_2 = %d°', in(idx)), 'FontSize', 9)
end
xlabel('X (mm)')
ylabel('Y (mm)')
title('Path Traced by Pin A')
legend('Path of Pin A', '\theta_2 = 0°', '\theta_2 = 120°', 'Location', 'best')
grid on
box on
axis equal
%exportgraphics(gcf, 'ENME_473_Project_3c.png', 'Resolution', 600);

%% Print tabular results (Q3c)
fprintf('\n=== Pin A Position (mm) ===\n');
fprintf('%10s %12s %12s\n', 'Theta 2', 'X_A', 'Y_A');
for k = 1:length(in)
    fprintf('%10d %12.2f %12.2f\n', in(k), Ax(k), Ay(k));
end

%% Display Angles Table
angles_table = table(in', ...
    theta3_vals', theta4_vals', theta5_vals', ...
    theta6_vals', theta7_vals', theta8_vals', ...
    'VariableNames', {'Input Angle','theta3','theta4','theta5','theta6','theta7','theta8'});
disp('Link Angles (deg)')
disp(angles_table)
