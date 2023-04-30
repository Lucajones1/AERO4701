function r_LGCV = ECEF2LGCV(rsat_ECEF,rgs_LLH)
r_E = 6378.137; %[km]
R = rgs_LLH(3)+r_E;
phi = rgs_LLH(2)*pi/180;
lambdai = rgs_LLH(1)*pi/180;
% Transforming the Groundstation's position vector to ECEF
rgs_ECEF = [R*cos(lambdai)*cos(phi);R*cos(lambdai)*sin(phi);R*sin(lambdai)];
% Finding the Relative vector
rdiff = rsat_ECEF-rgs_ECEF;
C = [-sin(lambdai)*cos(phi),-sin(phi),-cos(lambdai)*cos(phi);-sin(lambdai)*sin(phi),cos(phi),-cos(lambdai)*sin(phi);cos(lambdai),0,-sin(lambdai)];
% Transforming the Relative vector from ECEF to LGCV
r_LGCV = C'*rdiff;
end