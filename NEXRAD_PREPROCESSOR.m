%% Script to process NEXRAD netcdf scan files for use in fire perimeter tracking estimates
% Author: Neil P. Lareau, Ph.D. (nlareau@unr.edu)
% Univeristy of Nevada, Reno
% Department of Physics
% 

%%
close all
clear all
addpath('~/Dropbox/MATLAB/'); % this points the code to supporting MATLAB scripts and should be adjusted to point to your local MATLAB directory

%% Figure Flag
prompt = 'Display Figures?: yes or no (lower case)'; %this will tell the code which fire you're processing
str = input(prompt,'s');
if strcmp(str,'yes')
    pflag=1; %plot figures in code
else
    pflag=0; %don't plot figures in code
end

%% Constants
re=6.37*10^6;
daysec=24*60*60; %seconds per day
daymsec= daysec.*10^3;%millisecond per day (radar time is in msecs since midnight UTC
%%

prompt = 'Enter Fire Name: bear or camp (lower case):'; %this will tell the code which fire you're processing
str = input(prompt,'s');
disp(str)
if isempty(str)
    str = 'camp';
end
cases={str};

%% Start by loading the gridded radar reflecitivity data. This is the  reflectivity interpolated to an x,y,z grid defined by xvec, yvec,zvec.
if strcmp(cases,'camp')
    prompt = 'Navigate to Folder Containing Radar Data for the Camp Fire: hit "enter" to continue'; %this will tell the code which fire you're processing
    str = input(prompt,'s');%displays prompt to screen, requires user to hit enter
    dirname=uigetdir(); %user navigates to the data directory containing the radar data for this case
    cd(dirname); %change to the correct data director (needs to be changed for your machine)
    files=dir(strcat(dirname,'/KBBX_V06_2018*.nc')); %this directory pat will need to be changed to point to your local KBBX files for the camp fire time period
    fout='camp_fire_radar_for_perimeters.mat';

elseif  strcmp(cases,'bear')
    
    prompt = 'Navigate to Folder Containing Radar Data for the Bear Fire: hit "enter" to continue'; %this will tell the code which fire you're processing
    str = input(prompt,'s');%displays prompt to screen, requires user to hit enter
    dirname=uigetdir(); %user navigates to the data directory containing the radar data for this case
    cd(dirname); %change to the correct data director (needs to be changed for your machine)
    files=dir(strcat(dirname,'/KBBX_V06_2020*.nc'));%this directory pat will need to be changed to point to your local KBBX files for the bear fire time period
    fout='bear_fire_radar_for_perimeters.mat';
end
    
%% READ STATION METADATA FROM THE NETCDF FILE
StationLatitude = ncreadatt(files(1).name,'/','StationLatitude'); %read in the station latitude
StationLongitude = ncreadatt(files(1).name,'/','StationLongitude'); %read in the station longitude
StationElevationInMeters= ncreadatt(files(1).name,'/','StationElevationInMeters'); %read in the station elevation (in meters MSL)
scalefactor=(re.*cosd(StationLatitude)).*(pi./180); %used to convert longitude to x-distances
latfactor=111180; % meters per deg Latitude

% details of the Reflectivity data
%{
                       units                     = 'dBz'
                       long_name                 = 'Reflectivity_HI'
                       missing_value             = 0
                       signal_below_threshold    = 0
                       scale_factor              = 0.5
                       add_offset                = -33
                       _Unsigned                 = 'true'
                       SNR_threshold             = 16
                       range_folding_threshold   = 50
                       _CoordinateAxes           = 'timeR_HI elevationR_HI azimuthR_HI distanceR_HI'
                       scale_factor_metadata     = 0.5
                       add_offset_metadata       = -33
                       range_folded_value_packed = 1
                       range_folded_value        = -32.5
%}


%% INITIALIZE VARIABLES FOR STORING THE REFLECTIVITY DATA
aa=0; %counter variable
[dbZmaster,xmaster,ymaster,zmaster]=deal(nan(length(files)*2,1832,720)); %reflecitivyt matrix, note that we use 2x the # of files to accommodate the succesive scans as the same elevation (thus 2 times per file)
timemaster=nan(length(files)*2,1); %time vector for each scan
%% Construct a filter to remove islands of data
nf=zeros(9,9); %9x9 stencil
%set outer points to 1 
nf(1:3,:)=1; 
nf(7:9,:)=1;
nf(:,1:3)=1;
nf(:,7:9)=1;
nft=sum(nf(:)); %this is used to compute the fraction of the surrounding data that are NaNs.


%%
for ff=1:length(files)
    disp(ff./length(files));%prints to the screen what fraction complete the processing is (i.e., .3 = 30% of the way through the files
    fname=strcat(files(ff).folder,'/',files(ff).name); %create the full path to the input file
    
    % CONSTRUCT A TIME STAMP FROM THE FILE NAME
    YYYY=str2double(fname(end-17:end-14));
    MM=str2double(fname(end-13:end-12));
    DD=str2double(fname(end-11:end-10));
    HH=str2double(fname(end-8:end-7));
    MN=str2double(fname(end-6:end-5));
    timer=(datenum(YYYY,MM,DD,HH,MN,0));
    
    %% Read in High Resolution (SAILS) reflectivity data using low-level netcdf functions
    ncid = netcdf.open(fname,'NC_NOWRITE');
    %retreive SAILS reflectivity data (hi-resolution)
    varid = netcdf.inqVarID(ncid,'Reflectivity_HI');                                       
    dbZ_HI = double(netcdf.getVar(ncid,varid,'uint8'));                                      
    add_offset = netcdf.getAtt(ncid,varid,'add_offset');
    scale_factor = netcdf.getAtt(ncid,varid,'scale_factor');
    dbZ_HI(dbZ_HI==0)=nan;
    dbZ_HI = (double(dbZ_HI) * double(scale_factor)) + double(add_offset);
    netcdf.close(ncid);
    
    %% read in the scan metadata (azimuths, ranges, elevations) which can be read easily with MATLAB's high-level netcdf functions
    az_HI=ncread(fname,'azimuthR_HI'); 
    elev_HI=ncread(fname,'elevationR_HI');
    rdist_HI=ncread(fname,'distanceR_HI');
    time_HI=double(ncread(fname,'timeR_HI'));
    time_HI=time_HI./(daymsec)+floor(timer);

    
    %% FOR CAMP and BEAR FIRES we use the 1.5 deg beam which is the 3rd and 4th sweep of the NEXRAD data
    for ev=1:size(elev_HI,2)
        evnow=mode(elev_HI(:,ev)); % pull out the elevation for this scan.
        if evnow>1.3 & evnow<1.7        %if the VCP changes make sure it is still in this range
            aa=aa+1; %increment the counter
            time_now=squeeze(nanmean(time_HI(:,ev))); %mean time of the scan
            az1=az_HI(:,ev); %pull out the azimuths for that elevation
            [az1,sidx]=sort(az1); %sort the azimuths just to make sure they're allways in the same order
            dbZnow=squeeze(dbZ_HI(:,sidx,ev)); %now sort the Reflectivity in the same manner and parse out just this elevation scan
            evnow=mode(elev_HI(:,ev)); % pull out the elevation for this scan (scan bounces so I use the mode).
            ae=(4/3).*(re); %4/3 earth radius model for beam refraction
            z=(((rdist_HI.^2) + (ae.^2) + (2*rdist_HI.*ae.*(sind(evnow)))).^(1/2))-ae; %height of the radar beam above the radar site
            hdist=  ae.*(asin((rdist_HI.*cosd(evnow))./(ae+z)));% horizontal distance along the ground from the radar site
            [azmt,hdstmt]=meshgrid(az1,hdist); % create a mesh of azimuths and along ground distances
            [~,zdist]=meshgrid(az1,z); % z grid of data
            xdist=sind(azmt).*hdstmt; % x-component of the horizontal distance
            ydist=cosd(azmt).*hdstmt; % y-component of the horizontal distance
            
            %% Remove islands of data surrounded by NaNs
            dbzimg=dbZnow; %make a copy of the reflectivity
            dbzimg((dbZnow)>20)=1; %set points greater than 20dbz to 1
            dbzimg(isnan(dbZnow))=0; %everthing else is zero
            dbzimg=logical(dbzimg); %make sure it is a logical array (0s and 1s)
            BW2 = bwareafilt(dbzimg,[1500 inf]); %filter out small regions of connected pixels
            
            if pflag==1 
                %display the raw reflectivity
                figure(1);clf;
                subplot(1,2,1);cla;hold on;
                pcolor(xdist,ydist,dbZnow);shading flat;
                caxis([-20 50]);
                daspect([1 1 1]);
                xlim([-10 10]*10^4)
                ylim([-10 10]*10^4)
                box on; grid on;
                set(gca,'fontsize',15,'fontweight','bold','layer','top','linewidth',2);
            end    
            
            %perform the small object filtering and nanmask filtering
            dbZnow(BW2==0)=nan;%mask reflectivity data where the mask is zero
            nanmask=dbZnow; %now create a nan-mask
            nanmask(abs(nanmask)>0)=0; 
            nanmask(isnan(nanmask))=1;
            nanmaskf=filter2(nf,nanmask)./nft; %this should yield the percentage of adjacent data that are nans
            dbZnow(nanmaskf>.7)=nan; %if the local nan fraction in the 9x9 stencil of the boundary points is >70% toss this data point
            
            if pflag==1
                %plot the filtered data
                subplot(1,2,2);cla;
                pcolor(xdist,ydist,dbZnow);shading flat;
                caxis([-20 50]);
                daspect([1 1 1]);
                xlim([-10 10]*10^4)
                ylim([-10 10]*10^4)
                box on; grid on;
                set(gca,'fontsize',15,'fontweight','bold','layer','top','linewidth',2);
                pause(.1);
            end
            %%
            dbZnow(dbZnow<-20)=nan; %mask any very small values
            xdistlon=xdist./(scalefactor); %convert x-distance to longitude distances (not yet with the StationLongitude added in)
            ydistlat=ydist./(latfactor); %convert y-distance to latitude distances (not yet with the StationLatitude added in)
            dbZmaster(aa,:,:)=dbZnow; %store in dbZmaster for each sweep
            timemaster(aa)=time_now; %store time of each sweep
            xmaster(aa,:,:)=xdist; %store each sweep's x-distances
            ymaster(aa,:,:)=ydist; %store each sweep's y-distances
            zmaster(aa,:,:)=zdist; %store each sweep's z-heights; (note this doesn't include the station elevation, so it is height above the radar)
        end
    end
end

%% Optional Test Figure
if pflag==1
    dbzx=squeeze(nanmax(dbZmaster(:,:,:),[],1)); %plot the time-maximum radar reflectivity
    figure(1);clf;set(gcf,'pos',[100 100 1000 1000]);
    subplot(1,2,1);cla;hold on;
    pcolor(xdist,ydist,dbzx);shading flat;
    caxis([-20 50]);
    daspect([1 1 1]);
    xlim([-10 10]*10^4)
    ylim([-10 10]*10^4)
    box on; grid on;
    set(gca,'fontsize',15,'fontweight','bold','layer','top','linewidth',2);
    
    % we can use a use defined input polygon (xer,yer) to restrict data to
    % the domain of interest, which can speed the code up and remove
    % spurious detections outside of the fire area (e.g., terrain clutter)
    prompt='Would you like to restrict the data domain? yes or no';
    str = input(prompt,'s');
    if strcmp(str,'yes')
        [xer,yer]=ginput();
        inner=inpolygon(xdist,ydist,xer,yer);
        dbZmaster(:,inner==0)=nan; %set all points outside of the polygon (xer,yer) to be NaN
    end
    subplot(1,2,2);cla;hold on;
    dbzx=squeeze(nanmax(dbZmaster,[],1));  %compute the maximum reflectivity againt
    pcolor(xdist,ydist,dbzx);shading flat; %confirm that we've restricted the data as intended
    caxis([-20 50]);
    daspect([1 1 1]);
    xlim([-10 10]*10^4)
    ylim([-10 10]*10^4)
    box on; grid on;
    set(gca,'fontsize',15,'fontweight','bold','layer','top','linewidth',2);
    
end

%% OUTPUT THE RADAR DATA IN A '.MAT' FILE 

save(fout,'dbZmaster','timemaster','xmaster','ymaster','zmaster','StationElevationInMeters','StationLatitude','StationLongitude','-v7.3');
