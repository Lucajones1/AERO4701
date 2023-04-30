%% Polar Ground and Observation Plots
function polargroundplot(rsat_ECEF,no_periods)
% Initialise the Polar Ground Plots
figure(6);
clf;
clf reset;
figure(5);
clf;
clf reset;
hold on;
axis equal;
% Creating a vector from 0-360 deg
deg = [0:1:360];
% Defining and plotting the Circular Rings
xpole(1,:) = 30*sind(deg);
xpole(2,:) = 60*sind(deg);
xpole(3,:) = 90*sind(deg);
ypole(1,:) = 30*cosd(deg);
ypole(2,:) = 60*cosd(deg);
ypole(3,:) = 90*cosd(deg);
plot([0,0],[-90,90],'color',[192,192,192]/255);
plot([-90*sind(45),90*sind(45)],[-90*cosd(45),90*cosd(45)],'color',[192,192,192]/255);
plot([-90,90],[0,0],'color',[192,192,192]/255);
plot([90*sind(45),-90*sind(45)],[-90*cosd(45),90*cosd(45)],'color',[192,192,192]/255);
plot(xpole(1,:),ypole(1,:),'color',[192,192,192]/255);
plot(xpole(2,:),ypole(2,:),'color',[192,192,192]/255);
plot(xpole(3,:),ypole(3,:),'k');

% Plot Compass Directions and Angles on Figure
text(0,90,'N');
text(-90,0,'E');
text(0,-90,'S');
text(90,0,'W');
text(0,0,'90^o');
text(30*0.7,30*0.7,'60^o');
text(60*0.7,60*0.7,'30^o');
text(90*0.7,90*0.7,'0^o');
title('Groundplot of Cubesat from Groundstation');
gp1 = gca;
f6 = figure(6);
gp2 = copyobj(gp1,f6);
title('Groundplot of Cubesat from Groundstation with 10^\circ Elevation Loss');
figure(5);

% Defining position of chosen Groundstation: CUAVA Groundstation
long = -33.889; %[deg]
lat = 151.19; %[deg]
alt = 0.05; %[km]
rgs_LLH = [lat,long,alt];

%% Calculating the Observation vector
% Finding the Relative vector from the Groundstation's LLH vector and the
% Satellite's ECEF vector in the LGCV reference frame
rrel_LGCV = ECEF2LGCV(rsat_ECEF,rgs_LLH);
% Transforming LGCV Cartesian vector to a Polar vector
rrel_polar = cartesian2polar(rrel_LGCV);

% Converting Polar coordinates to x & y coordinates for the Polar Ground
% Plot
polarx = -(90 - rrel_polar(3,:)).*sind(rrel_polar(2,:));
polary = (90 - rrel_polar(3,:)).*cosd(rrel_polar(2,:));

% Working out what sections of the Polar Ground Plot are visible by the
% Groundstation
vis_polar = [];
vis = [];
vis_adj = [];
nvis_polar = [];
nvis = [];
n_val = [];
for n = 1:length(rrel_polar)
    if (polarx(n)^2+polary(n)^2)<=90^2
        vis_polar = [vis_polar,[polarx(n);polary(n)]];
        vis = [vis,rrel_polar(:,n)];
        if (polarx(n)^2+polary(n)^2)>=80^2 && polary(n)>=90*sind(235) && polary(n)<=0 && polarx(n)<=0
        else
            vis_adj = [vis_adj,[polarx(n);polary(n)]];
        end
        n_val = [n_val;n];
    else
        nvis_polar = [nvis_polar,[polarx(n);polary(n)]];
        nvis = [nvis,rrel_polar(:,n)];
    end
end 

% Plotting x & y coordinates onto the Polar Ground Plot
plot(vis_polar(1,:),vis_polar(2,:),'r.');

% Plotting the adjusted Groundplot that takes into account a 10 degree 
% elevation loss due to the Topology of CSG.
figure(6);
plot(vis_adj(1,:),vis_adj(2,:),'r.');

% Plotting the Elevation vs Range
figure(7);
clf;
clf reset;
v = scatter(vis(1,:),vis(3,:),'r.');
hold on;
n = scatter(nvis(1,:),nvis(3,:),'k.');
legend([v,n],'Visible by GS','Not Visible by GS');
xlabel('Range [km]');
ylabel('Elevation [deg]');
title('Elevation vs Range from Groundstation');

% Plotting the Azimuth vs Range
figure(8);
clf;
clf reset;
v = scatter(vis(1,:),vis(2,:),'r.');
hold on;
n = scatter(nvis(1,:),nvis(2,:),'k.');
legend([v,n],'Visible by GS','Not Visible by GS');
xlabel('Range [km]');
ylabel('Azimuth [deg]');
title('Azimuth vs Range from Groundstation');

% Determining the maximum elevation and percentage of time the satellite is
% visible for from the Groundstation
% If statement to prevent this loop from running unless the total time for
% the number of periods is greater than 24 hours in length
if no_periods >= 1
    percentvis = length(vis_polar)/length(rrel_polar);
    fprintf('Percentage of time the satellite is in view: %.2f%%\n',percentvis*100);
    adjpercentvis = length(vis_adj)/length(rrel_polar);
    fprintf('Adjusted percentage of time the satellite is in view: %.2f%%\n',adjpercentvis*100);
    maxelev = max(rrel_polar(3,:));
    fprintf('The maximum elevation of the satellite is: %.2f degrees\n',maxelev);
end
end