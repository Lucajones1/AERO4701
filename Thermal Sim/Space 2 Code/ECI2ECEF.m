function r_ECEF = ECI2ECEF(r_ECI,t)
omega_ie = 7.292155e-5;
[rows,cols] = size(r_ECI);
r_ECEF = zeros(size(r_ECI));
for n = 1:cols
    C = [cos(omega_ie*t(n)),sin(omega_ie*t(n)),0;-sin(omega_ie*t(n)),cos(omega_ie*t(n)),0;0,0,1];
    r_ECEF(:,n) = C*r_ECI(:,n);
end
end