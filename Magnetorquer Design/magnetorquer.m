%% Function that converts a PWM input to a Torque Output
function [tau,Pconsumed] = magnetorquer(pwm,axis)
% Simulating the Magnetorquer
if axis == 'x' || axis == 'y'
    % Defining the Magnetorquer Torque Generated
    tau = 2.65317e-07;
    % Defining the Maximum Power Consumption
    Pmax = 0.312649; %[W]
    % Adjusting the torque by the PWM signal
    tau = pwm*tau;
    % Adjusting the Power Consumption
    Pconsumed = pwm*Pmax;
elseif axis == 'z'
    % Defining the Magnetorquer Torque Generated
    tau = 2.57520e-07;
    % Defining the Maximum Power Consumption
    Pmax = 0.123075; %[W]
    % Adjusting the torque by the PWM signal
    tau = pwm*tau;
    % Adjusting the Power Consumption
    Pconsumed = pwm*Pmax;
else
    error('Define an axis for the magnetorquer (x, y or z)');
end
end