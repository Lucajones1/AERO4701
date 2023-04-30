%% Calculating the Eccentric Anomaly
function eanom = eccentric_anomaly(manom,ecc)
%Defining an acceptable tolerance for the solution
tol = 1e-9;

%Setting the starting approximation for the Eccentric Anomaly
if manom < pi
    eanom = manom + ecc/2;
else
    eanom = manom - ecc/2;
end

%Repeating the calculation of the Eccentric Anomaly until the tolerance is
%reached
check = 1;
while abs(check) > tol
    check = (eanom - ecc*sin(eanom) - manom)/(1 - ecc*cos(eanom));
    eanom = eanom - check;
end
end