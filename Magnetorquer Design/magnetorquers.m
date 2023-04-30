% Function for sizing the Magnetorquers
function tau_req = magnetorquers(sat)
%% Defining Satellite Design Parameters
% Pulling out the Inertias
Ix_max = sat.geom.Ix_max;
Iy_max = sat.geom.Iy_max;
Iz_max = sat.geom.Iz_max;
% Finding the Maximum Moment of Inertia as this will yield the largest 
% torque required
I_max = max([Ix_max,Iy_max,Iz_max]);

%% System Requirements
% Rotation Rate
omega_req = sat.reqs.ADCSreqrotrate; %[rad/s]
% Timeframe for Rotation
t_rot = sat.reqs.ADCStimeforrot; %[sec]

%% Design Calculations
% Calculating the required Rate of Decay of the Rotation Rate (ie Angular
% Deceleration)
alpha_req = omega_req/t_rot; %[rad/s^2]
% Defining the Torque required from the Magnetorquer
tau_req = alpha_req*I_max;
end