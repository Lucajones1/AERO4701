%% Earth 3D plot and Groundtrack
function earth_3Dgroundtrack(r_ECI,r_ECEF,r_LLH,epoch_day,dt,t_sgmt)
% Initialise the two 3D Earth Plots
figure(1);
clf;
clf reset;
figure(2);
clf;
clf reset;

% Defining the elliptical properties of Earth
semimajor = 6378.137; %[km]
semiminor = 6356.752; %[km]

% Forming the Ellipsoid and setting up Figure 1's plot
[ex,ey,ez] = ellipsoid(0,0,0,semimajor,semimajor,semiminor);
figure(1);
earth = surf(ex,ey,ez);
xlabel('x Displacement [km]');
ylabel('y Displacement [km]');
zlabel('z Displacement [km]');
title('ECI 3D Orbit Plot');
axis equal;
hold on;
set(gca,'Color','k');

% Placing the map of the Earth onto the ellipsoid surface
link = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
world = imread(link);
set(earth, 'FaceColor', 'texturemap', 'CData', world, 'FaceAlpha', 1, 'EdgeColor', 'none');
set(gca,'ZDir','reverse');

% Setting up Figure 2's plot
e1 = gca;
f2 = figure(2);
e2 = copyobj(e1,f2);
title('ECEF 3D Orbit Plot');
figure(1);

% Defining  constants that are needed for the 3D Earth Plot
omega_ie = 7.292115e-5; %[rad/s]

% Initialise the Groundtrack
figure(3);
clf;
clf reset;

% Placing the map of the Earth onto a 2D plot
image([-180,180],[-90,90],world);
title('Groundtrack of Orbit')
grid on;
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
hold on;


% Determining whether the plots are to be animated or not
choice = input('Animate the Plots? (Y/N)','s');

% Plotting stationary ECI & ECEF 3D Orbit and LLH Groundtrack
% ECEF
figure(2);
% Rotating the Earth by the time since 12:00am GMT
rotate(earth,[0,0,1],omega_ie*t_sgmt(1)*180/pi,[0,0,0]);
plot3(r_ECEF(1,:),r_ECEF(2,:),r_ECEF(3,:),'r');
% ECI
figure(1);
% Rotating the Earth by the time since 12:00am GMT
rotate(earth,[0,0,1],omega_ie*t_sgmt(1)*180/pi,[0,0,0]);
or = plot3(r_ECI(1,:),r_ECI(2,:),r_ECI(3,:),'r','LineWidth',2);
% Plotting Groundtrack LLH vector
figure(3);
scatter(r_LLH(2,:),r_LLH(1,:),'r.');
s1 = scatter(r_LLH(2,1),r_LLH(1,1),'ko');
e1 = scatter(r_LLH(2,end),r_LLH(1,end),'ks');
legend([s1,e1],'Start','End');

if lower(choice) == 'n'
    return;
end

% Initialising the animation for the 3D Orbit
figure(1);
x_sat = r_ECI(1,1);
y_sat = r_ECI(2,1);
z_sat = r_ECI(3,1);
s3d = plot3(x_sat,y_sat,z_sat,'go');
legend([or,s3d],'\color{white} Orbit','\color{white} Satellite');
set(s3d,'markersize',8,'linewidth',2);
set(s3d,'xdatasource','x_sat');
set(s3d,'ydatasource','y_sat');
set(s3d,'zdatasource','z_sat');

% Initialising the animation for the Groundtrack
figure(4);
clf;
clf reset;
image([-180,180],[-90,90],world);
title('Groundtrack of Orbit')
grid on;
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
hold on;

% The plot has to be divided into the trail and the current position of the
% satellite
% Trail component of the Groundtrack
lat_trail = r_LLH(1,1);
lon_trail = r_LLH(2,1);
sgt = plot(lon_trail,lat_trail,'r.');
set(sgt,'markersize',3,'linewidth',3);
set(sgt,'xdatasource','lon_trail');
set(sgt,'ydatasource','lat_trail');

% Current position component of Groundtrack
lat_sat = r_LLH(1,1);
lon_sat = r_LLH(2,1);
scp = plot(lon_sat,lat_sat,'go');
set(scp,'markersize',7,'linewidth',2);
set(scp,'xdatasource','lon_sat');
set(scp,'ydatasource','lat_sat');
legend([sgt,scp],'Groundtrack','Satellite');

% Animation code
% Number of frames to skip
f_skipped = 4;

% Looping through the dataset and animating the plot
[rows,cols] = size(r_LLH);
for n = 1:f_skipped:cols
    % Updating the 3D Plot
    x_sat = r_ECI(1,n);
    y_sat = r_ECI(2,n);
    z_sat = r_ECI(3,n);
    refreshdata(s3d,'caller');
    rotate(earth,[0,0,1],omega_ie*dt*f_skipped*180/pi,[0,0,0]);
    
    % Updating the Groundtrack
    % Trail
    lat_trail = r_LLH(1,1:n);
    lon_trail = r_LLH(2,1:n);
    set(sgt,'XData',lon_trail,'YData',lat_trail);
    % Satellite
    lat_sat = r_LLH(1,n);
    lon_sat = r_LLH(2,n);
    refreshdata(scp,'caller');
    drawnow;
end

end