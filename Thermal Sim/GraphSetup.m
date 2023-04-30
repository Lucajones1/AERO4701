%% Graph Setup Code
% Author: 490432831

function colours = GraphSetup(fsz)
%% Determining provided Variables
% Checking if an input is provided
if nargin < 1
    fsz         = 11;                       % Fontsize
end
%% Colours
% List of Colours to use within the Graph
black         = [0 0 0]/255;
gray          = [89 89 89]/255;
red           = [216 30 49]/255;
green         = [0 128 0]/255;
blue          = [27 99 165]/255;
yellow        = [251 194 13]/255;
cyan          = [2 169 226]/255;
purple        = [126 47 142]/255;

% Setting the Colour Order
set(groot,'defaultAxesColorOrder',[black;blue;red;purple;green;yellow;cyan;gray]);

%% Plotting Parameters
% Defining Plotting Parameters
alw             = 1;                        % AxesLineWidth
lw              = 2;                        % LineWidth
msz             = 10;                       % MarkerSize

% Setting the Plotting Parameters
set(0,'defaultAxesLineWidth',alw);
set(0,'defaultAxesFontSize',fsz);
set(0,'defaultLineLineWidth',lw);
set(0,'defaultLineMarkerSize',msz);
set(0,'defaultfigurecolor','w');

% Defining and setting the font
font            = 'latex';
set(groot,'defaulttextinterpreter',font);
set(groot,'defaultAxesTickLabelInterpreter',font);
set(groot,'defaultLegendInterpreter',font);
set(groot,'defaultLegendEdgeColor',black);
set(groot,'defaultLegendFontSize',0.85*fsz);

%% Outputting a Colour Matrix
colours = [black;blue;red;green;yellow;purple;cyan;gray];
end