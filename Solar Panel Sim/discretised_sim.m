%% Satellite Thermal Simulation
% This script solves the 4 different heat loads that the Cubesat is
% exposed to and then extracts a maximum and minimum value.
clear;
clc;

close all;

% Setting up the Plotting Environment
GraphSetup(13);

%% Standards
% Defining the elliptical properties of Earth
semimajor = 6378.137; %[km]
semiminor = 6356.752; %[km]
r_earth = semimajor;
% Earth's Albedo
b = 0.3;
% Solar Flux Constant
G_sol = 1367; %[W/m^2]
% Stefan-Boltzmann Constant
sigma = 5.67037e-8; %[W/m^2K^4]

%% Defining the Orbital Elements for the Cubesat
% Orbital Altitude
alt = 600; %[km]
% Semimajor Axis
sm_axis = semimajor+alt;
% Inclination
inc = 45; %[deg]
% Eccentricity
ecc = 0;
% Right Ascension of Ascending Node
rasc = 0; %[deg]
% Argument of Perigee
argp = 0; %[deg]
% Mean Anomaly
manom = 0; %[deg]
% Date and Time of Launch to Orbit
date_timevec = [0,0,0;20,3,2023;...
    0,0,0;22,6,2023]; %[hr:min:sec;day:month:year]
% Standard Gravitational Parameter for Earth
mu = 398600; %[km^3/s^2]
% Orbital Period
T = 2*pi*sqrt(sm_axis^3/mu); %[s]
t_eclipse = (120/360)*T;

% Looping through the two Date and Time vector values
for p = 1:2
date_time = date_timevec((2*p-1):(2*p),:);

% Getting the Time and Satellite's ECI and ECEF Position Vectors
[time,r_ECI,v_ECI,r_ECEF,r_sun,mlong_s,obliquity,theta] = getrsat(sm_axis,inc,ecc,rasc,argp,manom,date_time,1);
% Adding r_sun to the ECI Plot
figure(1);
plot3([0,r_sun(1,1)],[0,r_sun(2,1)],[0,r_sun(3,1)],'LineWidth',3);
axis([-10000,10000,-10000,10000,-10000,10000]);

% Redefining the Satellite's Semimajor Axis as the Orbital Radius (as ecc =
% 0)
r_orbit = sm_axis;
% Adjusting the True Anomaly Vector so that it goes from 0→180→360deg 
% rather than 0→180→-180→0deg
theta = [theta(1:end/2),theta(end/2+1:end)+2*pi].*180/pi; %[deg]

%% Defining Cubesat Parameters
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

% Satellite's Side Lengths
x_sat = l_side; %[mm]
y_sat = l_side; %[mm]
z_sat = l_height*unit_size; %[mm]
% Satellite's Mass
m_sat = m_std*unit_size; %[kg]

% Calculating the Areas of each Face
A_solpan = 0.1*0.08; %[m]
A_nsolpan = 0; %[m]
% Creating an Area Vector containing the Areas of each Face
A_sat = [A_solpan,A_solpan,A_solpan,A_solpan,A_nsolpan,A_nsolpan];

% Cubesat's Solar Panel Absorptivity
alpha = 0.8;

%% Primary Solar Radiation
% Heat load generated from Solar Radiation incident on the Cubesat

% Assumptions:
% - Maximum Area faces the Sun constantly (negates the need to worry about 
% Fview)
% - Sun's Heat Flux is constant and at its Maximum value

% Angle between r_earth and the Cubesat's Position Vector
chi_sat = acos(r_earth./vecnorm(r_ECI));
% Angle between r_earth and the Sun's Position Vector
chi_sun = acos(r_earth./vecnorm(r_sun));
% Angle between the Cubesat's and Sun's Position Vectors
chi = acos(dot(r_ECI,r_sun)./(vecnorm(r_ECI).*vecnorm(r_sun)));

% Looping through the Chi Vectors to define the Zeta Vector
zeta = [];
for i = 1:length(chi)
    % Defining the Zeta Vector (Representing whether the Cubesat is exposed
    % to Sun or Shade)
    if chi_sat(i)+chi_sun(i)<=chi(i)
        zeta = [zeta,0];
    else
        zeta = [zeta,1];
    end
end


% Verifying that the Satellite's time in the Sun is >= time spent in the
% Earth's Shade - Value should be one
fprintf('Does the Satellite Spend more time in the Sun: %d\n',...
    length(zeta(zeta==1))>=length(zeta(zeta==0)));

% Finding the View Factor and Area vectors for each Face
Fview = calcFview(semimajor,alt,inc,rasc,mlong_s,obliquity,T,time);
% Summing up the View Factors to get the Total View Factor across the True
% Anomaly vector
Fviewtot = Fview{1}+Fview{2}+Fview{3}+Fview{4}+Fview{5};

% Looping through the Area Vector
Q_solcomb = 0;
for i = 1:6
    % Finding the Heat Load Vector from the Primary Solar Radiation
    Q_sol{i} = G_sol*A_sat(i)*alpha.*Fview{i}.*zeta; %[W]
    % Changing Negative Load generation terms to 0
    Q_sol{i}(Q_sol{i}<0) = 0;
    % Combined Solar Power across all faces
    Q_solcomb = Q_solcomb+Q_sol{i};
end

% Plot showing the Combined Thermal Energy, incident on the Solar Panels
f1 = figure;
hold on;
plot(theta,Q_solcomb);
xlabel('True Anomaly [deg]');
ylabel('Thermal Power [W]');
legend('Combined Load','Location','north');
xlim([0,360]);

%% Converting the Thermal Power into Electrical Power Generated
% Defining the Solar Panel Parameters
% Efficiency
eta = 0.155;
% Maximum Power
P_max = 0.8689; %[W]
% Degradation Amount over 6-months as a Decimal
deg_amount = 0.025;
% Converting Degradation Amount to a rate per second
deg_rate = 1-deg_amount/(365/2*86400);
% Looping through the faces and calculating Power Generated
P_spcomb = 0;
for i = 1:6
    % Calculating the Electrical Power Generated
    P_solpan{i} = (eta*deg_rate).*Q_sol{i}; %[W]
    % Adjusting the Power Generated to account for the maximum output of a
    % Solar Panel
    P_solpan{i}(P_solpan{i}>P_max) = P_max; %[W]
end
% Multiplying the whole vector by 2 to account for both Solar Panels on
% either face
P_spcomb = P_spcomb.*2; %[W]

% Generating the plots
[f2,f3] = powerplots(P_solpan,theta,deg_amount);

% Saving the figures from the first loop
if p == 1
    figs1 = [f2,f3];
end
end

% Exporting the Plots
fnum = [f1,figs1,f2,f3];
fname = ["Q_solvstheta","Pgenvstheta_best","Pgencombvstheta_best",...
    "Pgenvstheta_worst","Pgencombvstheta_worst"];
GraphExportEPS(fnum,fname);

%% Functions
% Psi Function
function psivec = calcpsi(chi,b)
% Looping through the chi vector of angles and calculating the resulting
% psi value for each angle.
psivec = [];
for i = 1:length(chi)
    if chi(i) >= 0 && chi(i) <= pi/2
        psivec = [psivec,b*cos(chi(i))];
    else
        psivec = [psivec,0];
    end
end
end

% View Factor Function from Mike's Paper
function Fview = calcFview(r_earth,alt,inc,rasc,mlong_s,obliquity,T,time,...
    beta_crit)
% Beta Angle Vector - Angle between the Orbital Plane and the Sun Vector
betavec = asind(cosd(mlong_s).*(sind(rasc)*sind(inc))-...
    sind(mlong_s).*cosd(obliquity).*(cosd(rasc)*sind(inc))+...
    sind(mlong_s).*sind(obliquity).*cosd(inc));
% If the Critical Beta Angle is already defined, then don't recalculate it
if nargin < 9
    % Critical Beta Angle - Angle between the Orbital Plane and the Sun Vector
    % required to plunge part of the Orbit into an Eclipse (Shade)
    beta_crit = asin(r_earth/(r_earth+alt)).*180/pi; %[deg]
end
% Looping through the Beta Angle Vector and calculating the Eclipse Factor
% and then each face's View Factor
Fviewzen = [];
Fviewnad = [];
Fviewposv = [];
Fviewnegv = [];
FviewNaS = [];
for i = 1:length(time)
    % Defining the Piecewise Function
    if abs(betavec(i)) < beta_crit
        eclipse_fn(i) = 1/180*acosd(sqrt(alt^2+2*r_earth*alt)/...
            ((r_earth+alt)*cos(betavec(i)))); %[deg]
    elseif abs(betavec(i)) >= beta_crit
        eclipse_fn(i) = 0; %[deg]
    end
    % Calculating the View Factor of each Face
    % Zenith Side
    if T/4 > time(i) || time(i) > 3*T/4
        Fviewzen = [Fviewzen,cos(2*pi*time(i)/T)*cosd(betavec(i))];
    else
        Fviewzen = [Fviewzen,0];
    end
    % Nadir Side
    if (T/4 < time(i) && time(i) < T/2*(1-eclipse_fn(i))) ||...
            (T/2*(1+eclipse_fn(i)) < time(i) && time(i) < 3*T/4)
        Fviewnad = [Fviewnad,-cos(2*pi*time(i)/T)*cosd(betavec(i))];
    else
        Fviewnad = [Fviewnad,0];
    end
    % Positive Velocity Side
    if time(i) > T/2*(1+eclipse_fn(i))
        Fviewposv = [Fviewposv,-sin(2*pi*time(i)/T)*cosd(betavec(i))];
    else
        Fviewposv = [Fviewposv,0];
    end
    % Negative Velocity Side
    if time(i) < T/2*(1-eclipse_fn(i))
        Fviewnegv = [Fviewnegv,sin(2*pi*time(i)/T)*cosd(betavec(i))];
    else
        Fviewnegv = [Fviewnegv,0];
    end
    % North and South Sides
    if T/2*(1-eclipse_fn(i)) > time(i) && time(i) > T/2*(1+eclipse_fn(i))
        FviewNaS = [FviewNaS,sind(betavec(i))];
    else
        FviewNaS = [FviewNaS,0];
    end
end

% Saving the View Factor Vectors into cells for passing out of the Function
% Order of Faces:
% Zenith Face
% Nadir Face
% Positive Velocity Face
% Negative Velocity Face
% North Face
% South Face
Fview{1} = Fviewzen;
Fview{2} = Fviewnad;
Fview{3} = Fviewposv;
Fview{4} = Fviewnegv;
Fview{5} = FviewNaS;
Fview{6} = FviewNaS;

% Summing up the View Factors to get the Total View Factor vector
Fviewtot = Fview{1}+Fview{2}+Fview{3}+Fview{4}+Fview{5}+Fview{6};
% Dividing through by the Maximum View Factor to ensure the View Factor
% doesn't exceed 1
maxFview = max(Fviewtot);
% Dividing each individual View Factor by the Maximum View Factor as well
Fview{1} = Fview{1}./maxFview;
Fview{2} = Fview{2}./maxFview;
Fview{3} = Fview{3}./maxFview;
Fview{4} = Fview{4}./maxFview;
Fview{5} = Fview{5}./maxFview;
Fview{6} = Fview{6}./maxFview;
end