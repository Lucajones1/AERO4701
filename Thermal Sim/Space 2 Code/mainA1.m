% AERO3760 Assignment 1 Script
clear;
clc;

%% Orbital Elements of the Cubesat
% Satellite Data
sat_number = "BENTsat";
launch_year = 23;
launch_no = 1;
epoch_year = 23;
epoch_day = 79;

% First 5 Orbital Elements
inc = 45;
rasc = 0;
ecc = 0;
argp = 0;
manom = 0;

% Altitude
alt = 600; %[km]
% Defining the elliptical properties of Earth
semimajor = 6378.137; %[km]
semiminor = 6356.752; %[km]

% Final Orbital Element
sm_axis = semimajor+alt; %[km]

%% Calculating the State Vector from the Orbital Elements & performing Transformations
% Calculating the specific angular momentum
mu = 398600;
angm = sqrt(sm_axis*mu*(1-ecc^2));

% Defining the spatial step between each True Anomaly degree
dtheta = 1;

% Defining the number of periods for the plots
no_periods = 1*2.0774;

% Finding the timestep between each True Anomaly degree
period = (2*pi/sqrt(mu))*sm_axis^(3/2);
t_o = manom*(pi/180)*(period/(2*pi));
t_f = t_o + no_periods*period;

% Creating a time vector and defining its timestep
time = linspace(t_o,t_f,10000);
dt = time(2)-time(1);

% Calculating the Time since Vernal Equinox
% Sidereal Period of the Earth
siderealp = 23.9345; %[hr]
% Finding the days since Vernal Equinox
dayfrac = epoch_day - (31+28+20+(9*60^2+37*60)/(siderealp*60^2)); %[days]
% Converting days into seconds
t_gmt = dayfrac*siderealp*60^2; %[s]
% Adding 12 hours onto the time since Vernal Equinox to account for the 
% World Map starting at UTC-12 rather than UTC+0
t_gmt = t_gmt + (12*60^2);
% Adding the calibrated GMT time to the time vector
t_sgmt = time + t_gmt;

% Looping through values for the time vector and calculating the state 
% vector at the given times
r_ECI = [];
v_ECI = [];
r_ECEF = [];
for n = 1:length(time)
    %Calculating the variables required for the position vector
    t = time(n);
    M = ((2*pi)/period)*t;
    % Finding the Eccentric Anomaly for the current timestep
    E = eccentric_anomaly(M,ecc);
    % Finding the True Anomaly for the current timestep
    theta = 2*atan(tan(E/2)*sqrt((1+ecc)/(1-ecc)));
    % Calculating the Position Vector in the Perifocal Reference Frame
    rxbar(:,n) = angm^2/mu*(1/(1+ecc*cos(theta)))*[cos(theta);sin(theta);0];
    % Calculating the Velocity Vector in the Perifocal Reference Frame
    vxbar(:,n) = mu/angm*[-sin(theta);ecc+cos(theta);0];
    % DCM for Perifocal to ECI Reference Frame
    Q = [cosd(argp),sind(argp),0;-sind(argp),cosd(argp),0;0,0,1]*[1,0,0;0,cosd(inc),sind(inc);0,-sind(inc),cosd(inc)]*[cosd(rasc),sind(rasc),0;-sind(rasc),cosd(rasc),0;0,0,1];
    % Tranforming the State Vector from Perifocal to ECI
    r_ECI = [r_ECI,Q*rxbar(:,n)];
    v_ECI = [v_ECI,Q*vxbar(:,n)];
    % Transforming the ECI position vector to ECEF
    r_ECEF = [r_ECEF,ECI2ECEF(r_ECI(:,n),t_sgmt(n))];
end




% Accounting for the J2 effect
% Calculating Earth's perturbation
J2 = 0.00108263;
s_E = 6378.137; % Earth's Semimajor Axis [km]
pertb = -3/2*sqrt(mu)*J2*s_E^2/((1-ecc^2)^2*sm_axis^(7/2));
% Calculating the relevant Euler angle rotation rates
rasc_dot = pertb*cosd(inc);
argp_dot = pertb*((5/2)*sind(inc)^2 - 2);

% Looping through the State Vector and adjusting it to account for the J2
% effect
r_ECI_J2 = zeros(size(r_ECI));
v_ECI_J2 = r_ECI_J2;
r_ECEF_J2 = r_ECI_J2;
r_LLH = r_ECEF_J2;
for n = 1:length(time)
    % Calculating the new Euler Angles
    rasc_J2 = rasc + rasc_dot*(time(n)-t_o)*180/pi;
    argp_J2 = argp + argp_dot*(time(n)-t_o)*180/pi;
    % DCM for Perifocal to ECI Reference Frame and accounting for the J2 effect
    Q = [cosd(argp_J2),sind(argp_J2),0;-sind(argp_J2),cosd(argp_J2),0;0,0,1]*[1,0,0;0,cosd(inc),sind(inc);0,-sind(inc),cosd(inc)]*[cosd(rasc_J2),sind(rasc_J2),0;-sind(rasc_J2),cosd(rasc_J2),0;0,0,1];
    % Tranforming the State Vector from Perifocal to ECI
    r_ECI_J2(:,n) = Q*rxbar(:,n);
    v_ECI_J2(:,n) = Q*vxbar(:,n);
    % Transforming ECI to LLH for use on the Groundtrack
    r_ECEF_J2(:,n) = ECI2ECEF(r_ECI_J2(:,n),t_sgmt(n));
    r_LLH(:,n) = ECEF2LLH(r_ECEF_J2(:,n));
end

%% Generating the required Plots
% 3D Orbit plot and Groundtrack
earth_3Dgroundtrack(r_ECI,r_ECEF,r_LLH,epoch_day,dt,t_sgmt);
% Polar Ground Plot
polargroundplot(r_ECEF,no_periods);

% Saving all Data
save('data.mat');