% Function that models the Earth's Magnetic Field
function BCs = magBCs(sat)
%% Satellite Parameters
% Pulling out the Satellite's Mass and Maximum Area
m_sat = sat.geom.mass;
Amax = sat.geom.Amax;

%% Magnetic Field Boundary Conditions - NOAA
% Minimum
Bmin = 25e-6; %[T]
% Maximum
Bmax = 65e-6; %[T]

%% Solar Heat Flux Constants - from Thermal Sim
% Minimum
G_solmin = 1367; %[W/m^2]
% Maximum
G_solmax = 1428; %[W/m^2]

%% Gravitational Acceleration
% Earth Semimajor and Semiminor Axes
% Semiminor
semiminor = 6356.752; %[km]
% Semimajor
semimajor = 6378.137; %[km]

% Newton's Gravitational Constant
G = 6.67430e-11; %[Nm^2/kg^2]
% Earth's Mass
m_e = 5.972e24; %[kg]
% Minimum
gmin = G*m_e/semimajor^2; %[m/s^2]
% Maximum
gmax = G*m_e/semiminor^2; %[m/s^2]

%% Atmospheric Drag - J77' Model
% Assumptions:
% - Using the Shock Drag formula
% - Assuming Cd to be constant across the entire orbit
% - Assuming Atmospheric Density to be constant across the entire orbit

% Taking the Exospheric Temperature to be 900K (Nominal Temp). This will be
% the minimum BC. For the maximum BC,the exospheric temperature will be
% 1600K (Extreme Solar Emissions)
% http://www.braeunig.us/space/atmmodel.htm
% http://www.braeunig.us/space/pdf/Atmosphere_120-2K.pdf

% Minimum Atmospheric Density at 600km for the assumed Exospheric 
% Temperature (500K)
rho_600min = 1.24e-14; %[kg/m^3]
% Maximum Atmospheric Density at 600km for the assumed Exospheric 
% Temperature (1600K)
rho_600max = 1.82e-12; %[kg/m^3]

% Drag Coefficient for CubeSat - Vertical Flat Plate
% https://www.agi.com/getmedia/280be1f7-7d1a-49f3-8d78-3492ce527d8c/A-Critical-Assessment-of-Satellite-Drag-and-Atmospheric-Density-Modeling.pdf?ext=.pdf)
Cd = 3.04;

% Velocity of CubeSat - Norm of velocity vector generated in the thermal
% sim
v_sat = 7.55786e3; %[m/s]

% Calculating the minimum and maximum Drag Force on half of the Cubesat
F_dragmin = 1/2*Cd*m_sat*Amax/2*rho_600min*v_sat^2;
F_dragmax = 1/2*Cd*m_sat*Amax/2*rho_600max*v_sat^2;

% Compiling the BCs into a single output parameter
BCs = [Bmin,Bmax;G_solmin,G_solmax;gmin,gmax;F_dragmin,F_dragmax];
end