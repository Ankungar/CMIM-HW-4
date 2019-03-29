%Computational Methods in Mechanics
%Homework set 4, exercise 1
%Daphne van Dijken

clear all
close all
clc

syms t

%Initial values
a = 0.1; %meter
b = 0.2; %meter
omega = 1; %rad/s
tol = 1e-4; %Tolerance
x0 = [0.2527; 0.2803]; %Initial condition, found when running at t = 0.


%Function f with time dependence
f = @(x,t) [a*cos(pi/6 + omega*t) + b*cos(x(1)) - x(2);
             a*sin(pi/6 + omega*t) - b*sin(x(1));];

%Jacobian J
J = @(x) [-b*sin(x(1)), -1;
          -b*cos(x(1)),0];

dfdt = matlabFunction(diff(f(x0,t)));

%Run NR method with time dependence
x_time = [];
x_dot = [];
x1 = x0;
for t_val = 0:0.01:50
    [x1, n1] = NR_method_t(f,J,x1,t_val,tol);
    x_time = [x_time, x1];
    x_dot = [x_dot, -1*J(x1)\dfdt(t_val)];
end

t_val = [0:0.01:50];
theta = x_time(1,:);
theta_dot = x_dot(1,:);
d = x_time(2,:);
d_dot = x_dot(2,:);

%% Plot results
figure
hold on
plot(t_val,theta)
plot(t_val,theta_dot)
xlabel('time (s)')
ylabel('theta (rad)')
legend('theta','theta dot')
hold off

figure
hold on
plot(t_val,d)
plot(t_val,d_dot)
xlabel('time (s)')
ylabel('d (m)')
legend('d','d dot')
hold off

