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
G_sol = 1428; %[W/m^2]
% G_sol = 1367; % Typical
% Earth Flux Constant
G_earth = 237; %[W/m^2]
% Stefan-Boltzmann Constant
sigma = 5.67037e-8; %[W/m^2K^4]
% Outer Space Ambient Temperature
Tinf_DS = 2.7; %[K]

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
date_time = [0,0,0;20,3,2023]; %[hr:min:sec;day:month:year]
% Standard Gravitational Parameter for Earth
mu = 398600; %[km^3/s^2]
% Orbital Period
T = 2*pi*sqrt(sm_axis^3/mu); %[s]
t_eclipse = (120/360)*T;
disp(T);
disp(t_eclipse);
disp(t_eclipse/T); 

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
Azen = z_sat*x_sat; %[m]
Anad = Azen; %[m]
Aposv = z_sat*y_sat; %[m]
Anegv = Aposv; %[m]
ANaS = x_sat*y_sat; %[m]
% Creating an Area Vector containing the Areas of each Face
A_sat = [Azen,Anad,Aposv,Anegv,ANaS,ANaS];
% Calculating the Total Surface Area of the Cubesat
SA_sat = sum(A_sat); %[m^2]

% Cubesat's Average Absorptivity
alpha = 0.7;
% Cubesat's Average Emissivity
epsilon = 0.5;

%% Primary Solar Radiation
% Heat load generated from Solar Radiation incident on the Cubesat

% Assumptions:
% - Maximum Area faces the Sun constantly (negates the need to worry about 
% Fview)
% - Sun's Heat Flux is constant and at its Maximum value

% Calculating the Cubesat's Area facing the Sun (Maximised)
% A_sat = sqrt(x_sat^2+y_sat^2)*sqrt(z_sat^2+x_sat^2);
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
for i = 1:length(A_sat)
    % Finding the Heat Load Vector from the Primary Solar Radiation
    Q_sol{i} = G_sol*A_sat(i)*alpha.*Fview{i}.*zeta; %[W]
    % Changing Negative Load generation terms to 0
    Q_sol{i}(Q_sol{i}<0) = 0;
end

% Plot of View Factors for Verification - No Lines should be above 1 and
% there should be a clear point of all View Factors going to 0
f1 = figure;
plot(theta,Fview{1});
hold on;
plot(theta,Fview{2});
plot(theta,Fview{3});
plot(theta,Fview{4});
plot(theta,Fview{5});
xlabel('True Anomaly [deg]');
ylabel('View Factor');
legend('Zenith Face','Nadir Face','Pos $v$ Face','Neg $v$ Face',...
    'N and S Faces','Location','north');
xlim([0,360]);

% Plots of Thermal Power just from Primary Solar Radiation for each face 
% and combined
[f2,f3] = powerplots(Q_sol,theta);

%% Secondary (Reflected) Solar Radiation
% Heat load generated from Solar Radiation reflected off the Earth, 
% incident on the Cubesat

% Assumptions:
% - Maximum Area faces the Sun constantly (negates the need to worry about 
% Fview)

% Solving for the Psi vector using the chi vector from above
psivec = calcpsi(chi,b);
% Finding the View Factor vectors for each Face with a Critical Beta Angle
% of 0 degrees (ie no eclipse) → it is the psi function that adds the
% eclipse into this example
Fview = calcFview(semimajor,alt,inc,rasc,mlong_s,obliquity,T,time,0);
% Summing up the View Factors to get the Total View Factor across the True
% Anomaly vector
Fviewtot = Fview{1}+Fview{2}+Fview{3}+Fview{4}+Fview{5};

% Looping through the Area Vector
for i = 1:length(A_sat)
    % Finding the Heat Load vector from the Primary Solar Radiation
    Q_solref{i} = G_sol*A_sat(i)*alpha.*Fview{i}.*psivec; %[W]
end

% Plots of Thermal Power just from Reflected Solar Radiation for each face 
% and combined
[f4,f5] = powerplots(Q_solref,theta);

%% Blackbody Radiation from the Earth
% Heat load generated from Radiation emitted from the Earth as its a hot
% body (ie T > 0K). This Radiation uses the same Absorptivity as the
% Primary Solar Radiation and same View Factor as the Secondary Solar 
% Radiation

% Assumptions:
% - Maximum Area faces the Earth constantly
% - Earth's Heat Flux is constant and at its Maximum value

% Finding the View Factor vectors for each Face with a Critical Beta Angle
% of 0 degrees (ie no eclipse)
Fview = calcFview(semimajor,alt,inc,rasc,mlong_s,obliquity,T,time,0);
% Summing up the View Factors to get the Total View Factor across the True
% Anomaly vector
Fviewtot = Fview{1}+Fview{2}+Fview{3}+Fview{4}+Fview{5};

% Looping through the Area Vector
for i = 1:length(A_sat)
    % Finding the Heat Load Vector from the Primary Earth Radiation
    Q_earth{i} = G_earth*A_sat(i)*alpha.*Fview{i}; %[W]
end

% Plots of Thermal Power just from Primary Earth Radiation for each face 
% and combined
[f6,f7] = powerplots(Q_earth,theta);

%% Thermal Load generated by Onboard Electronics
% Heat load generated by the Onboard Electronics in the Cubesat

% Assumptions:
% - Generated Heat Load is distributed uniformly throughout the Cubesat

% Total Power consumed by the Onboard Electronics (=Thermal Heat generated)
P_onboard = 0; %[W]
% Dividing the Total Power consumed by the Total Surface Area to get the 
% Heat Flux
q_onboard = P_onboard/SA_sat; %[W/m^2]
% Transforming the Heat Flux into a vector of constant magnitude
q_onboard = q_onboard.*ones(size(theta)); %[W/m^2]

% Looping through the Area Vector
for i = 1:length(A_sat)
    % Finding the Heat Load on each Face
    Q_onboard{i} = q_onboard.*A_sat(i); %[W]
end

%% Cubesat Total Heat Load Plots
% Plot showing the Heat Loading for each face
f8 = figure;
hold on;
% Summing all individual Heat Load components together
for i = 1:6
    Q_tot{i} = Q_sol{i}+Q_solref{i}+Q_earth{i}+Q_onboard{i};
    % Plotting the Total Heat Load for each face
    plot(theta,Q_tot{i});
end

[minZen,maxZen] = bounds(Q_tot{1});
[minNad,maxNad] = bounds(Q_tot{2});
[minPos,maxPos] = bounds(Q_tot{3});
[minNeg,maxNeg] = bounds(Q_tot{4});
[minNor,maxNor] = bounds(Q_tot{5});
[minSou,maxSou] = bounds(Q_tot{6});

xlabel('True Anomaly [deg]');
ylabel('Thermal Power [W]');
% title('Total Radiation on Each Face');
legend('Zenith Face','Nadir Face','Pos $v$ Face','Neg $v$ Face',...
    'North Face','South Face','Location','north');
xlim([0,360]);

% Plot showing the Combined Heat Load for the Cubesat
f9 = figure;
hold on;
% Summing all Total Heat Loads for each Face together
Q_combtot = Q_tot{1}+Q_tot{2}+Q_tot{3}+Q_tot{4}+Q_tot{5}+Q_tot{6};
% Plotting the Combined Heat Load
plot(theta,Q_combtot);
xlabel('True Anomaly [deg]');
ylabel('Thermal Power [W]');
% title('Total Radiation on the Cubesat');
legend('Combined Load','Location','north');
xlim([0,360]);

%% Temperature Distribution of Cubesat over its Orbit
% Calculating the Temperature of the Cubesat, resulting from each Face's 
% Heat Load, using the Stefan-Boltzmann Law

% Assumptions:
% - Satellite is at Thermal Steady State (allows us to use the Stefan
% Boltzmann Law) 
% ↑ The above assumption CAN NOT be made, invalidating the plot below!

% Setting up a figure for use later
figure;
hold on;

% Looping through each Face
for i = 1:6
    % Calculating the Face's Temperature using the Stefan-Boltzmann Law
    T_combtot{i} = (Q_tot{i}./(epsilon*sigma*A_sat(i))+...
    Tinf_DS^4).^(1/4); %[K]
    % Plotting the Face's Temperature Distribution
    plot(theta,T_combtot{i});
end

% Labelling the Plot
xlabel('True Anomaly [deg]');
ylabel('Temperature [K]');
title('Temperature of Each Face');
legend('Zenith Face','Nadir Face','Pos $v$ Face','Neg $v$ Face',...
    'North Face','South Face','Location','north');
xlim([0,360]);


% Exporting the Plots
fnum = [f1,f2,f3,f4,f5,f6,f7,f8,f9];
fname = ["viewfactor","Qsolvstheta","Qsolcombvstheta","Qsolrefvstheta",...
    "Qsolrefcombvstheta","Qearthvstheta","Qearthcombvstheta",...
    "Qtotvstheta","Qtotcombvstheta"];
GraphExportThermal(fnum,fname);

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