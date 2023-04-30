% Function to get the Satellite's Position Vector from the Orbital
% Elements. This function also transforms the Normal Vectors of each Face
% of the Cubesat to the same reference frames as the Position Vector.
function [time,r_ECI_J2,v_ECI_J2,r_ECEF_J2,r_sun,mlong_s,epsilon,theta] = getrsat(sm_axis,inc,ecc,rasc,argp,manom,date_time,no_periods)
%% Orbital Elements of the Cubesat
% Determing whether rasc, manom or/and argp are provided
if nargin <= 5
    manom = 0; %[deg]
end
if nargin <= 4
    argp = 0; %[deg]
end
if nargin <= 3
    rasc = 0; %[deg]
end

% Pulling out the Date and Time
% Time
hours = date_time(1,1);
mins = date_time(1,2);
sec = date_time(1,3);
% Universal Time (UT)
UT = (hours-12)/24+mins/1440+sec/86400;
% Date
day = date_time(2,1);
month = date_time(2,2);
year = date_time(2,3);
% Calculating the Julian Date
JD = 367*year-((7*year+(month+9)/12)/4)+275*month/9+day+UT+1721013.5;
% Number of Days since J2000
nd = JD-2451545; %[days]

% Calculating days since Vernal Equinox
epoch_day = 0; %[days]
i = month;
leap_year = ~any(mod(year-2000,4)/4);
while i ~= 1
    if ~isempty(find([1,3,5,7,8,10,12]==i,1))
        epoch_day = epoch_day+31; %[days]
    elseif ~isempty(find([4,6,9,11]==i,1))
        epoch_day = epoch_day+30; %[days]
    else
        epoch_day = epoch_day+28+leap_year; %[days]
    end
    i = i-1;
end
% Adding final Number of Days for the last month
epoch_day = epoch_day+day; %[days]

%% Calculating the State Vector from the Orbital Elements & performing Transformations
% Calculating the Specific Angular Momentum
mu = 398600;
angm = sqrt(sm_axis*mu*(1-ecc^2));

% Determining the Number of Periods
if nargin <= 6
    no_periods = 1*2.0774;
end

% Finding the timestep between each True Anomaly degree
period = (2*pi/sqrt(mu))*sm_axis^(3/2);
t_o = manom*(pi/180)*(period/(2*pi));
t_f = t_o + no_periods*period;

% Creating a time vector and defining its timestep
time = linspace(t_o,t_f,10000);
dt = time(2)-time(1);

% Calculating the Time since Vernal Equinox
% Sidereal Period of the Earth (Hours)
siderealp_h = 23.9345; %[hr]
% Sidereal Period of the Earth (Days)
siderealp_d = 365.25636; %[days]
% Finding the days since Vernal Equinox
dayfrac = epoch_day - (31+28+20+(9*60^2+37*60)/(siderealp_h*60^2)); %[days]
% Converting days into seconds
t_gmt = dayfrac*siderealp_h*60^2; %[s]
% Adding 12 hours onto the time since Vernal Equinox to account for the 
% World Map starting at UTC-12 rather than UTC+0
t_gmt = t_gmt + (12*60^2);
% Adding the calibrated GMT time to the time vector
t_sgmt = time + t_gmt;

% Determing the Sun's Position via the Paper's Method
% Distance from Earth to Sun (1 AU)
dist_e2s = 1.4959787e8; %[km]
% Creating Number of Days Vector
nd = nd*ones(size(time))+time./(siderealp_h*60^2);
% Mean Anomaly of the Sun
manom_s = 357.529+0.98560023.*nd; %[deg]
% Mean Longitude of the Sun
mlong_s = 280.459 + 0.98564736.*nd; %[deg]
% Ecliptic Longitude
lambda = mlong_s+1.915.*sind(manom_s)+0.02.*sind(2.*manom_s); %[deg]
% Obliquity
epsilon = 23.439-3.56e-7.*nd;
% Unit Vector from the Earth to the Sun
uhat = [cosd(lambda);...
    cosd(epsilon).*sind(lambda);...
    sind(epsilon).*sind(lambda)];
% Mangitude of the r_sun vector
magr_sun = ((1.00014-0.01671.*cosd(manom_s))-(0.000140.*cosd(2.*manom_s)))*dist_e2s;
% Calculating the r_sun vector
r_sun = magr_sun.*uhat;

% Looping through values for the time vector and calculating the state 
% vector at the given times
r_ECI = [];
v_ECI = [];
r_ECEF = [];
theta = [];
for n = 1:length(time)
    % Calculating the variables required for the position vector
    t = time(n);
    M = ((2*pi)/period)*t;
    % Finding the Eccentric Anomaly for the current timestep
    E = eccentric_anomaly(M,ecc);
    % Finding the True Anomaly for the current timestep
    theta = [theta,2*atan(tan(E/2)*sqrt((1+ecc)/(1-ecc)))]; %[rad]
    % Calculating the Position Vector in the Perifocal Reference Frame
    rxbar(:,n) = angm^2/mu*(1/(1+ecc*cos(theta(n))))*...
        [cos(theta(n));sin(theta(n));0];
    % Calculating the Velocity Vector in the Perifocal Reference Frame
    vxbar(:,n) = mu/angm*[-sin(theta(n));ecc+cos(theta(n));0];
    % DCM for Perifocal to ECI Reference Frame
    Q = [cosd(argp),sind(argp),0;-sind(argp),cosd(argp),0;0,0,1]*...
        [1,0,0;0,cosd(inc),sind(inc);0,-sind(inc),cosd(inc)]*...
        [cosd(rasc),sind(rasc),0;-sind(rasc),cosd(rasc),0;0,0,1];
    % Transforming the State Vector from Perifocal to ECI
    r_ECI = [r_ECI,Q*rxbar(:,n)];
    v_ECI = [v_ECI,Q*vxbar(:,n)];
    % Transforming the ECI position vector to ECEF
    r_ECEF = [r_ECEF,ECI2ECEF(r_ECI(:,n),t_sgmt(n))];
end

%% Accounting for the J2 effect
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
    Q = [cosd(argp_J2),sind(argp_J2),0;-sind(argp_J2),cosd(argp_J2),0;0,0,1]*...
        [1,0,0;0,cosd(inc),sind(inc);0,-sind(inc),cosd(inc)]*...
        [cosd(rasc_J2),sind(rasc_J2),0;-sind(rasc_J2),cosd(rasc_J2),0;0,0,1];
    % Tranforming the State Vector from Perifocal to ECI
    r_ECI_J2(:,n) = Q*rxbar(:,n);
    v_ECI_J2(:,n) = Q*vxbar(:,n);
    % Transforming ECI to ECEF
    r_ECEF_J2(:,n) = ECI2ECEF(r_ECI_J2(:,n),t_sgmt(n));
    r_LLH(:,n) = ECEF2LLH(r_ECEF_J2(:,n));
end

%% Generating the required Plots
% 3D Orbit plot and Groundtrack
earth_3Dgroundtrack(r_ECI,r_ECEF,r_LLH,epoch_day,dt,t_sgmt);
end