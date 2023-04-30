function GraphExport(fig,filename)
%% Creating the new Working Directory
% Working out the current Working Directory
wd = pwd;
% Folder to Erase from, from the Working Directory
reffrom = '\code';
% Finding the index of the Code Folder in the Working Directory
idx = strfind(wd,reffrom);
% Erasing the Code Folder and Subfolders out of the Working Directory
wd = wd(1:idx-1);
% Navigating to the Images Folder
refto = '\images';
% Creating the new Working Directory
wd = [wd,refto];

%% Exporting the Figure/s to the Images Folder
% Setting up Function to take an Array of Figures + Filenames
for i = 1:length(fig)
    % Printing to a .eps File
    print(fig(i),[wd,'\',convertStringsToChars(filename(i))],'-depsc','-r0');
end
end