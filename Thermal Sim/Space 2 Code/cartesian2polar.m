function r_polar = cartesian2polar(r_cart)
[rows,cols] = size(r_cart);
r_polar = zeros(rows,cols);
for n = 1:cols
    R = norm(r_cart(:,n));
    psi = atan2(r_cart(2,n),r_cart(1,n))*180/pi;
    theta = atan(-r_cart(3,n)/sqrt(r_cart(1,n)^2+r_cart(2,n)^2))*180/pi;
    % Creating a Polar vector of [range;azimuth;elevation]
    r_polar(:,n) = [R,psi,theta];
end
end