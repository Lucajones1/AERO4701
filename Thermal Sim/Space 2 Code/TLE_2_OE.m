%% TLE to Orbital Elements Function
function [sat_name,TLE_line_1,TLE_line_2] = TLE_2_OE(filename)
%Opening the file and creating a fileID for the file so fgetl can extract
%the data from each row
file_number = fopen(filename, 'rb');

%Initialising counts for the different lines of the TLE codes
count0 = 1;
count1 = 1;
count2 = 1;

%Defining the standard gravitational parameter so the Mean Motion can be
%converted into the Semi-Major Axis
mu = 398600;

%Getting the three lines from the textfile for the first satellite
tle_line0 = fgetl(file_number);
tle_line1 = fgetl(file_number);
tle_line2 = fgetl(file_number);

%Initialise a while loop to go through the entire textfile and group the
%TLEs into different variables
while ischar(tle_line0)
    
    %Name of next object in the TLE codes
    next_name = convertCharsToStrings(tle_line0(~isspace(tle_line0)));
    %Creating an if statement that allows for the user to filter out any
    %debris or boosters that aren't part of the desired TLE data
    if strcmp(tle_line0(1),'I')
        %Provided the if statement if true, the data from the TLE code is
        %extracted for the given satellite
        
        %Extracting the Satellite Name from the first row
        sat_name(count0,1) = convertCharsToStrings(tle_line0(~isspace(tle_line0)));
        count0 = count0 + 1;
        
        %Extracting the Satellite Number, Epoch Year & Epoch Day
        sat_number(count1,1) = str2num(tle_line1(3:7));
        launch_year(count1,1) = str2num(tle_line1(10:11));
        launch_no(count1,1) = str2num(tle_line1(12:14));
        epoch_year(count1,1) = str2num(tle_line1(19:20));
        epoch_day(count1,1) = str2num(tle_line1(21:32));
        normalised_sat_no(count1,1) = count1;
        count1 = count1 + 1;
        
        %Extracting the Inclination (degrees), Right Ascension (degrees),
        %Eccentricity, Argument of Perigee (degrees), Mean Anomaly (degrees) &
        %Semimajor Axis (km). Note the semimajor axis (km) has been calculated
        %through conversion of the mean motion (rev/day)
        inclination(count2,1) = str2num(tle_line2(9:16));
        r_ascension(count2,1) = str2num(tle_line2(18:25));
        eccentricity(count2,1) = str2num(strcat('0.',tle_line2(27:33)));
        arg_perigee(count2,1) = str2num(tle_line2(35:42));
        mean_anomaly(count2,1) = str2num(tle_line2(44:51));
        sm_axis(count2,1) = (mu/(str2num(tle_line2(53:63))*2*pi/(24*60*60))^2)^(1/3);
        count2 = count2 + 1;
    else
    end
    
    %If there is more than one satellite, get the next 3 lines of TLE data
    %for the next satellite
    tle_line0 = fgetl(file_number);
    tle_line1 = fgetl(file_number);
    tle_line2 = fgetl(file_number);
end

%Compiling two matrices to export from the function. One for TLE line 1 
%and one for TLE line 2. The satellite names have been separated from the
%TLE matrices and exported as a separate variable
TLE_line_1 = [sat_number,launch_year,launch_no,epoch_year,epoch_day,normalised_sat_no];
TLE_line_2 = [inclination,r_ascension,eccentricity,arg_perigee,mean_anomaly,sm_axis];
end