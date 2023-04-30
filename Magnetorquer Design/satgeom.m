% Function that calculates the Moments of Inertia for the Satellite
function sat = satgeom
%% Defining Standard Cubesat Parameters
% Standard Cubesat Size
unit_size = 2; %[U]
% Standard Side Length
l_side = 100; %[mm]
l_side = l_side/1000; %[m]
% Standard Height Length
l_height = 113.5; %[mm]
l_height = l_height/1000; %[m]
% Standard Length Uncertainty
l_delta = 0.1; %[mm]
l_delta = l_delta/1000; %[m]
% Standard Mass
m_std = 1; %[kg/U]

%% Working out Satellite Parameters
% Satellite's Side Lengths
x_sat = l_side; %[mm]
y_sat = l_side; %[mm]
z_sat = l_height*unit_size; %[mm]
% Satellite's Mass
m_sat = m_std*unit_size; %[kg]
% Satellite's Areas
Axy = x_sat*y_sat; %[m^2]
Ayz = x_sat*y_sat; %[m^2]
Axz = x_sat*y_sat; %[m^2]
Amax = sqrt(x_sat^2+y_sat^2)*sqrt(x_sat^2+z_sat^2);

%% Satellite's Mass Moment of Inertia
% Assumptions: Uniform Density across the Cubesat (assumes Neutral Plane
% lies at the centre of the Cubesat)
% Calculating Mass Moment of Inertia for the 3 different Axes of Rotation
% that the Rigid Body (Cubesat) can rotate around
% General Formula for a Cuboid: Ix = 1/12*m*(y^2+z^2)
Ix = 1/12*m_sat*(y_sat^2+z_sat^2); %[kgm^2]
Iy = 1/12*m_sat*(z_sat^2+x_sat^2); %[kgm^2]
Iz = 1/12*m_sat*(x_sat^2+y_sat^2); %[kgm^2]
% Propagating the Uncertainty to the Second Moments of Area
Ix_delta = 1/12*m_sat*((2*l_delta/y_sat*y_sat^2)+(2*l_delta/z_sat*z_sat^2));
Iy_delta = 1/12*m_sat*((2*l_delta/z_sat*z_sat^2)+(2*l_delta/x_sat*x_sat^2));
Iz_delta = 1/12*m_sat*((2*l_delta/x_sat*x_sat^2)+(2*l_delta/y_sat*y_sat^2));
% Finding the Maximum value each Second Moment of Area can hold
Ix_max = Ix+Ix_delta;
Iy_max = Iy+Iy_delta;
Iz_max = Iz+Iz_delta;

%% Compiling the Satellite Design Parameters into one output parameter
% Inertias
sat.geom.Ix_max = Ix_max;
sat.geom.Iy_max = Iy_max;
sat.geom.Iz_max = Iz_max;
sat.geom.Ix = Ix;
sat.geom.Iy = Iy;
sat.geom.Iz = Iz;
% Dimensions
sat.geom.x_sat = x_sat;
sat.geom.y_sat = y_sat;
sat.geom.z_sat = z_sat;
sat.geom.Axy = Axy;
sat.geom.Ayz = Ayz;
sat.geom.Axz = Axz;
sat.geom.Amax = Amax;
% Mass
sat.geom.mass = m_sat;
% Size
sat.geom.size = unit_size;
% CubeSat Standards
sat.geom.stds.lw_unit = l_side;
sat.geom.stds.height_unit = l_height;
sat.geom.stds.m_unit = m_std;
end