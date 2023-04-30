%% Function transforming ECEF to LLH
% Input must be a 3xn matrix
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