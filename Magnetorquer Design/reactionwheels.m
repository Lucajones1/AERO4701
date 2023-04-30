%% Script for sizing the Reaction Wheels
clear;
clc;

close all;

GraphSetup(12);

%% Defining Satellite Parameters
% Loading Moments of Inertia data from the Icalc.m script
load('data.mat');
% Finding the Maximum Moment of Inertia as this will yield the largest 
% torque required
I_max = max([Ix_max,Iy_max,Iz_max]);
I_max = 0.00192;

%% Initial Design Parameters
% Reaction Wheel initial dimension estimates
% Disk Radius
d_disk = 30; %[mm]
d_disk = d_disk/1000; %[m]
% Disk Thickness
t_disk = 10; %[mm]
t_disk = t_disk/1000; %[m]

% Material Densities
% 304 Stainless Steel
rho_SS = 7930; %[kg/m^3]

% Loaded Motor Speed
omega_m = 4500; %[RPM]
omega_m = omega_m*2*pi/60; %[rad/s]

%% System Requirements
% Rotation Rate 
omega_req = 90; %[deg/s]
% Timeframe for Rotation
t_rot = 2; %[days]
% Converting units to SI
omega_req = omega_req*pi/180; %[rad/s]
t_rot = 2*24*60*60; %[sec]

%% Design Calculations
% Calculating the Reaction Wheel System Masses
m_disk = pi/4*d_disk^2*t_disk*rho_SS; %[kg]

% Calculating the Reaction Wheels' Mass Moments of Inertia
I_disk = 1/2*m_disk*(d_disk/2)^2;  %[kgm^2]

% Calculating the Maximum Angular Momentum vector that can be generated
% with the current design
L_rw = I_disk*omega_m; %[kgm^2/s]

% Calculating the resulting Rotation Rate of the Satellite
omega_sat = L_rw/I_max; %[rad/s]

% Finding the time taken to recover from the Tumble Rate Mission 
% Requirement
t_detumble = omega_req/omega_sat; %[sec]
fprintf('Time taken to detumble the Cubesat: %.2f seconds\n',t_detumble);
