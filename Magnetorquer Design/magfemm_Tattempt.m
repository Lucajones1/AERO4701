% Script to run FEMM 4.2
clear;
clc;

close all;

GraphSetup(14);

% Adding the FEMM 4.2 Directories
addpath('C:\femm42\mfiles'); % MATLAB Functions
addpath('C:\femm42\bin'); % Executable

%% Design Parameters for the Magnetorquers
% Coil Parameters Number of Turns
nturns = 515; %[Turns]
% Coil Current
I = 0.05; %[A]
% Coil Voltage
V = 5; %[VDC]
% Coil Resistance
R = 70.4237; %[Ω]
% Wire Diameter - 38 AWG
dia_w = 0.1; %[mm]

% Magnetorquer Structure Design Tube Outer Diameter
OD_tube = 20; %[mm]
% Tube Outer Radius
ro_tube = OD_tube/2; %[mm]
% Tube Wall Thickness
t_tube = 2; %[mm]
% Tube Inner Diameter
ID_tube = OD_tube-2*t_tube; %[mm]
ri_tube = ID_tube/2; %[mm]
% Tube Length
l_tube = 62; %[mm]
% Tube Relative Permeability - Assumed as being Air
mur_tube = 1;

% Defining the External Magnetic Field for the simulation
B_e = 50e-6; %[T]
% Permeability of Free Space
mu0 = 1.2566e-6; %[N/A^2]
% Converting the External Magnetic Field to a Coercivity
Hc = 2*B_e/mu0;


%% Building the Geometry of the Magnetorquer
% Points will be defined in a clockwise direction from the first point,
% which is the top left most point of the geometry

% Copper Coil Geometry
P1 = [ro_tube,l_tube/2];
P2 = [P1(1)+dia_w,P1(2)];
P3 = [P2(1),-P1(2)];
P4 = [P1(1),P3(2)];
coords_coil = [P1;P2;P3;P4];

% Magnetorquer 3D-Printed Tube Geometry - U Shaped
P1 = [ri_tube,l_tube/2+t_tube];
P2 = [coords_coil(2,1),P1(2)];
P3 = [P2(1),P1(2)-t_tube];
P4 = [ro_tube,P3(2)];
P5 = [P4(1),-P4(2)];
P6 = [P2(1),P5(2)];
P7 = [P6(1),P6(2)-t_tube];
P8 = [P1(1),P7(2)];
coords_3dtube = [P1;P2;P3;P4;P5;P6;P7;P8];

% Working out a suitable radii for the Air Gap arcs and the BC arc
r_airgap_i = 1.2*sqrt(coords_coil(2,1)^2+coords_coil(2,2)^2);
r_airgap_o = 1.25*sqrt(coords_coil(2,1)^2+coords_coil(2,2)^2);
r_BC = 1.5*sqrt(coords_coil(2,1)^2+coords_coil(2,2)^2);

% Defining the Block Label Coordinates
% 3D Printed Tube
B1 = [(coords_3dtube(1,1)+coords_3dtube(4,1))/2,0];
% Copper Coil
B2 = [(coords_coil(1,1)+coords_coil(2,1))/2,0];
% Magnetic Field surrounding the Magnetorquer
B3 = [1.5*coords_coil(2,1),0];
% Air Gap surrounding the Magnetorquer
B4 = [(r_airgap_i+r_airgap_o)/2,0];
% Magnetic Field surrounding the Air Gap
B5 = [1.1*r_airgap_o,0];

%% Running FEMM 4.2
% Creating a filename for the created FEMM files
filename = 'MagnetorquerFEMM';

% Opening the FEMM 4.2 Program
openfemm;
% Setting up the Magnetostatics Document: 0 → Magnetics Problem
newdocument(0);
% Defining the Problem Type - 'axi' → Assymetrical Problem (ie rotate drawn
% geometry 360deg)
mi_probdef(0,'millimeters','axi',1e-8,0,-30);

% Drawing the Coil Using the Rectangle function to draw the Coil
mi_drawrectangle(coords_coil(1,1),coords_coil(1,2),...
    coords_coil(3,1),coords_coil(3,2));

% Drawing the Tube Drawing the Lines
mi_drawline(coords_3dtube(1,1),coords_3dtube(1,2),...
    coords_3dtube(2,1),coords_3dtube(2,2));
mi_drawline(coords_3dtube(2,1),coords_3dtube(2,2),...
    coords_3dtube(3,1),coords_3dtube(3,2));
mi_drawline(coords_3dtube(3,1),coords_3dtube(3,2),...
    coords_3dtube(4,1),coords_3dtube(4,2));
mi_drawline(coords_3dtube(4,1),coords_3dtube(4,2),...
    coords_3dtube(5,1),coords_3dtube(5,2));
mi_drawline(coords_3dtube(5,1),coords_3dtube(5,2),...
    coords_3dtube(6,1),coords_3dtube(6,2));
mi_drawline(coords_3dtube(6,1),coords_3dtube(6,2),...
    coords_3dtube(7,1),coords_3dtube(7,2));
mi_drawline(coords_3dtube(7,1),coords_3dtube(7,2),...
    coords_3dtube(8,1),coords_3dtube(8,2));
mi_drawline(coords_3dtube(8,1),coords_3dtube(8,2),...
    coords_3dtube(1,1),coords_3dtube(1,2));

% Drawing the Arcs bounding the Air Gap
% First Arc at r_BC
mi_drawarc(0,-r_airgap_i,0,r_airgap_i,180,1);
% Second Arc at r_BC*1.1
mi_drawarc(0,-r_airgap_o,0,r_airgap_o,180,1);
% Third Arc for the Periodic BC
mi_drawarc(0,-r_BC,0,r_BC,180,1);


% Adding the 5 Block Labels
mi_addblocklabel(B1(1),B1(2));
mi_addblocklabel(B2(1),B2(2));
mi_addblocklabel(B3(1),B3(2));
mi_addblocklabel(B4(1),B4(2));
mi_addblocklabel(B5(1),B5(2));

% 
% % Creating the Boundary Arcs (7 total) and defining the Boundary
% % Conditions: 0 → Dirichlet BCs
% mi_makeABC(7,r_BC,0,0,0);

% Defining the Circuit Properties for the Coil: 1 → Series Circuit
mi_addcircprop('Coil',I,1);

% Defining Material Properties to be used by the Block Label Nodes
mi_addmaterial('Air',1,1,0,0,0,0,0,1,0,0,0);
mi_addmaterial('MagField',1,1,Hc,0,0,0,0,1,0,0,0);
mi_addmaterial('Tube',mur_tube,mur_tube,0,0,0,0,0,1,0,0,0); 
mi_addmaterial('CuWire',1,1,0,0,58,0,0,1,3,0,0,1,dia_w);

% Defining each Block Label's corresponding material Enabling Automesh
automesh = 1;
% Defining the 3D Printed Tube Block Label Properties
mi_selectlabel([B1(1),B1(2)]);
mi_setblockprop('Tube',automesh,0,'<None>',0,0,0);
mi_clearselected();
% Defining the Coil Block Label Properties
mi_selectlabel([B2(1),B2(2)]);
mi_setblockprop('CuWire',automesh,0,'Coil',0,0,nturns);
mi_clearselected();
% Defining the Inner Magnetic Field Block Label Properties
mi_selectlabel([B3(1),B3(2)]);
mi_setblockprop('Air',automesh,0,'<None>',0,0,0);
% mi_setblockprop('MagField',automesh,0,'<None>',0,0,0);
mi_clearselected();
% Defining the Air Block Label Properties
mi_selectlabel([B4(1),B4(2)]);
mi_setblockprop('<No Mesh>',automesh,0,'<None>',0,0,0);
mi_clearselected();
% Defining the Outer Magnetic Field Block Label Properties
mi_selectlabel([B5(1),B5(2)]);
mi_setblockprop('Air',automesh,0,'<None>',0,0,0);
% mi_setblockprop('MagField',automesh,0,'<None>',0,0,0);
mi_clearselected();

% Defining the Periodic BCs
mi_addboundprop('ExtMag',0,0,0,0,0,0,0,0,4,0,0);
mi_addboundprop('AirGap',0,0,0,0,0,0,0,0,6,0,0);
% Setting up the External Magnetic Field Periodic BC
mi_selectarcsegment(0,r_BC);
mi_setarcsegmentprop(1,'ExtMag',0,0);
mi_clearselected();
% Setting up the Air Gap Periodic BCs
mi_selectarcsegment(0,r_airgap_i);
mi_setarcsegmentprop(5,'AirGap',0,0);
mi_clearselected();
mi_selectarcsegment(0,r_airgap_o);
mi_setarcsegmentprop(5,'AirGap',0,0);
mi_clearselected();

% Setting Zoom to show Final Geometry
mi_zoomnatural();

% Saving the Geometry
mi_saveas(strcat(filename,'.fem'));

%% Analysis of Geometry
% Running the analysis
mi_analyze();
% Loading the solution
mi_loadsolution();

% Simulating the Closing of the Valve Number of Intervals to discretise the
% simulation to
it_no = 10;
% Dividing the Stroke Length by the Number of Iterations
ds = s_sol/it_no;
% Creating a vector of Forces and Plunger Positions
svec = zeros(it_no+1,1);
Fvec = zeros(it_no+1,1);

% Creating a vector for the Coil Circuit Properties
props = zeros(it_no+1,3);

% Initialising Plot used in Loop
f1 = figure('Position',[100,100,1500,950]);

% Looping through the Plunger Positions
for i = 1:it_no+1
    if i ~= 1
        % Running the Analysis for the new Position
        mi_analyze();
        mi_loadsolution();
    end
    % Calculating the Block Integral for the Annular Plunger (Group 1)
    mo_groupselectblock(1);
    svec(i) = ds*(i-1);
    Fvec(i) = mo_blockintegral(19);
    mo_clearblock();
    % Translating the Plunger to the new Position
    mi_selectgroup(1);
    mi_movetranslate(0,ds);
    mi_clearselected();

    %% B Field Distribution Plot
    % Looking at the Flux Distribution over the Working Gap of the SV
    % Creating a Position Vector over the Gap
    rposvec = coords_VHAP(1,1):0.05:coords_VHAP(2,1);
    % Finding the vector's size
    [rowsr,colsr] = size(rposvec);
    % Initialising a y vector of equal size
    yposvec = zeros(rowsr,colsr);
    % Solving for the B Field vector - Note the outputted B vector is of
    the form Bvec=[Br(:,1),Bz(:,1)]
    Bvec = mo_getb(rposvec,yposvec);
    % Finding the magnitude of the B Field Vector
    Bvec = vecnorm(Bvec,2,2); % Second entry is 2 by default. Third entry specifies norm of each row
    % Saving the Coil Circuit Properties
    props(i,:) = mo_getcircuitproperties('Coil');
    % Extracting the maximum B field value for the first iteration for the
    % Density Contour Plot
    if i == 1
        Bmax = 1.3*max(Bvec);
    end

    % Subplot showing the change in the B Field Distribution over the
    % Contact Surface of the Annular Plunger as the Valve moves from Closed
    % to Open Plotting the B vector vs the Normalised Position Vector
    subplot(3,4,i);
    % Creating the Normalised Position Vector
    rposvec = rposvec-rposvec(1);
    plot(rposvec,Bvec);
    grid on;
    xlim([0,rposvec(end)]);
    
    % Resizing the FEMM Window Size
    main_resize(900,900);
    % Printing a Density Plot of the Magnetic Field
    mo_zoomnatural();
    mo_showdensityplot(1,0,Bmax,0,'mag');
    mo_savebitmap(['Bdistvalve',num2str(i),'_',VHMaterial,'.jpg']);
end

% Annotating Plot
figure(f1);
axs = axes;
xlbl = xlabel('Distance along the Contact Surface of the Plunger (L→R) [mm]');
ylbl = ylabel('B Field Distribution [T]');
titl = title('B distribution along the Contact Surface of the Plunger');
set(axs,'Visible','off');
set(xlbl,'Visible','on');
set(ylbl,'Visible','on');
set(titl,'Visible','on');

%% Force vs Plunger Displacement Plot
f2 = figure;
% Flipping Force Vector for Plot
Fvec = flip(Fvec);
% Finding Relative Position vector of Plunger
rposplunger = s_sol-svec;
plot(rposplunger,Fvec);
xlabel('Relative Position of Plunger from Valve Open Position [mm]');
ylabel('Induced Force [N]');
title('Plunger Force vs Displacement');
% Highlighting last point as being Valve Closed Position
hold on;
plot(rposplunger(1),Fvec(1),'ro');
text(rposplunger(1)*0.8,Fvec(1)+max(Fvec)*0.05,'Valve Closed');

% Creating the Plot Names
f1name = ['Bdistvss','_',VHMaterial];
f2name = ['Fvss','_',VHMaterial];

% Exporting the Plots
tofldr = '\Results\Images';
imgwd = GraphExport([f1,f2],[f1name,',',f2name],tofldr);
% Moving the Density Contour Plots to the images folder
movefile('*.jpg',imgwd);

%% Saving the Circuit Properties
% Calculating Flux/Current, Resistance and Power
B_I = props(:,3)./props(:,1);
V_I = props(:,2)./props(:,1);
IV = props(:,1).*props(:,2);

% Concatenating the 3 calculated parameters onto the Coil Circuit
% Properties
props = [props,B_I,V_I,IV];

% Defining the Text File Name
fname_res = ['CoilCircuitProps','_',VHMaterial,'.txt'];
% Opening the Text File
fid = fopen(fname_res,'w');
% Printing the props parameter into the Text File - Scientific Notation
% chosen for alignment reasons
fprintf(fid,'%-s %- s %- s %- s %- s %- s\n',"Current   ", "V Drop    ", "Flux Link ", "Self Ind  ", "Resist    ", "Power");
fprintf(fid,'%-.3E %- .3E %- .3E %- .3E %- .3E %- .3E\n', props);
% Closing the Text File
fclose(fid);

% Finding the Directory to move the Text File to
reswd = strrep(imgwd,'\Images','');
% Moving the Text File to the Results Folder
movefile(fname_res,reswd);

%% Saving the Data Structure
% NOTE: Only uncomment this section if using MATLAB Defining the XML File
% Name
fname_dstr = ['Parameters','_',VHMaterial,'.xml'];
writestruct(params,fname_dstr);
% Moving the XML File to the Results Folder
movefile(fname_dstr,reswd);
