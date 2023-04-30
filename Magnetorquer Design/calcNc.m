%% Function to calculate the Control Torque from the Magnetorquers
function Mvec = calcNc(r_sat,omega)
% Returning an Error if the sizes of the inputs are wrong
[rowsr,colsr] = size(r_sat);
[rowsw,colsw] = size(omega);
if colsr ~= 1 || colsw ~= 1 || rowsr ~= 3 || rowsw ~= 3
    error('One/Both of the inputs are not 3x1 column vectors');
end

% Defining Orbital Parameters
% Period
T = 5.801235e+03; %[sec]
% Inclination
inc = 45; %[deg]
% Definign Staellite Parameters
% Minimum Mass Moment of Inertia
I_min = 0.00334; %[kgm^2]

% Converting the Satellite's Position vector to Longitude, Latitude and
% Height
r_LLH = ECEF2LLH(r_sat);
% Setting the date as the Vernal Equinox
dyear = decyear(2023,3,20);
% Calculating the Earth's Magnetic Field vector using IGRF-13
bvec = igrfmagm(r_LLH(3,1),r_LLH(1,1),r_LLH(2,1),dyear,13);
% Converting the Magnetic Field vector from nT to T
bvec = bvec/10^9; %[T]
% Defining the Control gain
kw = 2*(2*pi/T)*(1+sind(inc))*I_min;
% Calculating the Magnetic Dipole Moment
mvec = -kw/norm(bvec)*cross(bvec,omega');
% Control Torque
Mvec = cross(mvec,bvec);
end

%% Ancillary Functions
% Function transforming ECEF to LLH
function r_LLH = ECEF2LLH(r_ECEF)
[rows,cols] = size(r_ECEF);
r_LLH = zeros(rows,cols);
r_E = 6378.137;
for n = 1:cols
    R = norm(r_ECEF(:,n));
    alt = R-r_E;
    lambdai = asin(r_ECEF(3,n)/R);
    lambdai = 180/pi*lambdai;
    phi = atan2(r_ECEF(2,n),r_ECEF(1,n));
    phi = 180/pi*phi;
    % Creating the LLH vector of [Latitude;Longitude;Altitude]
    r_LLH(:,n) = [lambdai;phi;alt];
end
end