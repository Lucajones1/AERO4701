% Function that simulates the Magnetorquer in FEMM
function mag = magfemm(mag,inp_fname)
% Adding the FEMM 4.2 Directories
addpath('C:\femm42\mfiles'); % MATLAB Functions
addpath('C:\femm42\bin'); % Executable

%% Design Parameters for the Magnetorquers
% Coil Number of Turns
nturns = mag.coil.nturns; %[Turns]
% Coil Number of Layers
nlayers = mag.coil.nlayers; %[Layers]
% Coil Current
I = mag.coil.I; %[A]
% Coil Voltage
V = mag.coil.V; %[VDC]
% Coil Resistance
R = mag.coil.R; %[Ω]
% Wire Diameter - 38 AWG
dia_w = mag.coil.dia_w*1000; %[mm]

% Magnetorquer Structure Design Tube Outer Diameter
OD_tube = mag.tube.OD*1000; %[mm]
% Tube Outer Radius
ro_tube = mag.tube.r_out*1000; %[mm]
% Tube Wall Thickness
t_tube = mag.tube.t*1000; %[mm]
% Tube Inner Diameter
ID_tube = mag.tube.ID*1000; %[mm]
ri_tube = mag.tube.r_in*1000; %[mm]
% Tube Length
L_tube = mag.tube.L*1000; %[mm]
% Tube Relative Permeability - Assumed as being Air
mur_tube = mag.tube.mur;

% Defining the External Magnetic Field for the simulation
B_e = 50e-6; %[T]
% Permeability of Free Space
mu0 = 1.2566e-6; %[N/A^2]
% Converting the Uniform External Magnetic Field to a Coercivity. Note this
% is done by using the relationship between the Magnetic Field and
% Coercivity for a Sphere.
Hc = 3/2*B_e/mu0;

%% Building the Geometry of the Magnetorquer
% Points will be defined in a clockwise direction from the first point,
% which is the top left most point of the geometry

% Copper Coil Geometry
P1 = [ro_tube,L_tube/2];
P2 = [P1(1)+(dia_w*nlayers),P1(2)];
P3 = [P2(1),-P1(2)];
P4 = [P1(1),P3(2)];
coords_coil = [P1;P2;P3;P4];

% Magnetorquer 3D-Printed Tube Geometry - U Shaped
P1 = [ri_tube,L_tube/2+t_tube];
P2 = [coords_coil(2,1),P1(2)];
P3 = [P2(1),P1(2)-t_tube];
P4 = [ro_tube,P3(2)];
P5 = [P4(1),-P4(2)];
P6 = [P2(1),P5(2)];
P7 = [P6(1),P6(2)-t_tube];
P8 = [P1(1),P7(2)];
coords_3dtube = [P1;P2;P3;P4;P5;P6;P7;P8];

% Defining the Block Label Coordinates Air surroudning the Magnetorquer
B1 = [2*coords_coil(2,1),0];
% 3D Printed Tube
B2 = [(coords_3dtube(1,1)+coords_3dtube(4,1))/2,0];
% Copper Coil
B3 = [(coords_coil(1,1)+coords_coil(2,1))/2,0.2*(L_tube/2-(t_tube+1))];

%% Running FEMM 4.2
% Creating a filename for the created FEMM files
filename = ['MagnetorquerFEMM',inp_fname];

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
% Drawing a Node for the Contour Integral across from the Tube
mi_addnode(0,coords_3dtube(1,2));

% Adding the 8 Block Labels
mi_addblocklabel(B1(1),B1(2));
mi_addblocklabel(B2(1),B2(2));
mi_addblocklabel(B3(1),B3(2));

% Defining the Circuit Properties for the Coil: 1 → Series Circuit
mi_addcircprop('Coil',I,1);

% Defining Material Properties to be used by the Block Label Nodes
mi_addmaterial('Air',1,1,0,0,0,0,0,1,0,0,0);
mi_addmaterial('Tube',mur_tube,mur_tube,0,0,0,0,0,1,0,0,0); 
mi_addmaterial('CuWire',1,1,0,0,58,0,0,1,3,0,0,1,dia_w);

% Defining each Block Label's corresponding material Enabling Automesh
automesh = 1;
% Defining the Air Block Label Properties
mi_selectlabel([B1(1),B1(2)]);
mi_setblockprop('Air',automesh,0,'<None>',0,0,0);
mi_clearselected();
% Defining the 3D Printed Tube Block Label Properties
mi_selectlabel([B2(1),B2(2)]);
mi_setblockprop('Tube',automesh,0,'<None>',0,0,0);
mi_clearselected();
% Defining the Coil Block Label Properties
mi_selectlabel([B3(1),B3(2)]);
mi_setblockprop('CuWire',automesh,0,'Coil',0,0,nturns*nlayers);
mi_clearselected();

% Calculating a radius to place the Improvised Asymptotic Boundary
% Conditions (IABCs) at
r_BC = 2*sqrt(coords_coil(2,1)^2+coords_coil(2,2)^2);
% Creating the Boundary Arcs (7 total) and defining the Boundary
% Conditions: 0 → Dirichlet BCs
mi_makeABC(7,r_BC,0,0,0);
% % Applying the calculated Coercivity above, for the External Magnetic
% % Field, to the BC materials u1→u7
% mi_modifymaterial('u1',3,Hc);
% mi_modifymaterial('u2',3,Hc);
% mi_modifymaterial('u3',3,Hc);
% mi_modifymaterial('u4',3,Hc);
% mi_modifymaterial('u5',3,Hc);
% mi_modifymaterial('u6',3,Hc);
% mi_modifymaterial('u7',3,Hc);
% % Rotating the Magnetic Field Direction for each of the BC materials
% mi_selectlabel(64,13);
% mi_setblockprop('u1',automesh,0,'<None>',90,0,0);
% mi_selectlabel(61,25.5);
% mi_setblockprop('u2',automesh,0,'<None>',90,0,0);
% mi_selectlabel(56,37);
% mi_setblockprop('u3',automesh,0,'<None>',90,0,0);
% mi_selectlabel(48,48);
% mi_setblockprop('u4',automesh,0,'<None>',90,0,0);
% mi_selectlabel(38,57);
% mi_setblockprop('u5',automesh,0,'<None>',90,0,0);
% mi_selectlabel(25.5,65);
% mi_setblockprop('u6',automesh,0,'<None>',90,0,0);
% mi_selectlabel(14,69);
% mi_setblockprop('u7',automesh,0,'<None>',90,0,0);

% Resizing the FEMM Window Size
main_resize(900,900);
% Setting Zoom to show Final Geometry
mi_zoomnatural();

% Saving the Geometry
mi_saveas(strcat(filename,'.fem'));
% Saving an Image of the Geometry
mi_savebitmap(['FEMMGeometry',inp_fname,'.jpg']);

%% Analysis of Geometry
% Running the analysis
mi_analyze();
% Loading the solution
mi_loadsolution();

% Creating the contour for the Magnetic Field Line Integral
mo_addcontour(0,coords_3dtube(1,2));
mo_addcontour(coords_3dtube(1,1),coords_3dtube(1,2));
% Integrating over the contour, using the B.n integral, to give the Normal
% Flux out of the Magnetorquer and the Average B Field
intres= mo_lineintegral(0);
% Normal Flux
phi_n = intres(1); %[Wb]
% Average B Field
B_avg = intres(2); %[T]

% Extracting Coil Circuit Properties
coilprops = mo_getcircuitproperties('Coil');
% Current carried by Coil
I_femm = coilprops(1); %[A]
% Voltage Drop over Coil
V_femm = coilprops(2); %[VDC]
% Coil's Flux Linkage
fluxlink = coilprops(2); %[Wb]
% Coil's Power Consumption
P_femm = I_femm*V_femm; %[W]

% Pulling out the Maximum B value
maxB = mo_getb(coords_3dtube(1,1)/2,0);
maxB = norm(maxB);

% Printing a Density Plot of the Magnetic Field
mo_zoomnatural();
mo_showdensityplot(1,0,round(1.1*maxB,5),0,'mag');
mo_savebitmap(['FEMMContour',inp_fname,'.jpg']);

% Working out the current Working Directory
wd = pwd;
% Folder to Erase from, from the Working Directory
reffrom = '\Subsystems\ADCS\MATLAB Scripts';
% Finding the index of the Code Folder in the Working Directory
idx = strfind(wd,reffrom);
% Erasing the Code Folder and Subfolders out of the Working Directory
wd = wd(1:idx-1);
% Navigating to the Images Folder
refto = '\Assignments\Critical Design Report\images';
% Creating the new Working Directory
wd = [wd,refto];
% Moving the Density Contour Plots to the images folder
movefile('*.jpg',wd);

% Closing FEMM
closefemm();

% Saving FEMM Design Parameters into mag structure
mag.femm.phi_n = phi_n;
mag.femm.B_avg = B_avg;
mag.femm.I = I_femm;
mag.femm.V = V_femm;
mag.femm.fluxlinkage = fluxlink;
mag.femm.P = P_femm;
end
