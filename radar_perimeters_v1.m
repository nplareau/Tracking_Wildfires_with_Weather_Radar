%% FIRE PERIMETER TRACKING FROM NEXRAD RADAR REFLECTIVITY
% produced by Neil P. Lareau, Ph.D., University of Nevada, Reno (nlareau@unr.edu)
%
% OVERVIEW: 
% This script produces radar perimeter estimates from radar
% reflecitivity data. A description of the approach is contained in Lareau
% et al. 2022 (JGRA): "Tracking Wildfires with Weather Radar"
%
% INPUTS:
%  (1) *firename*_radar_for_perimeters.mat: This is a ".mat" MATLAB file containing
%  preprocessed radar observations for either the Camp, Bear, or King
%  Fires. These data are produced in an accompanying script ()
%     The contents of the file are as follows:
%     dbZmaster: a 3D array of reflectivity values structured as (time,xdist,ydist)
%     timemaster: a vector listing the times in UTC corresponding to the first dimension of dbZmaster
%     xmaster: a 2D array of the east-west distance in meters from the radar site (positive to the east)
%     ymaster: a 2D array of the north-south distance in meters from the radar site (positive to the north)
%     zmaster: a 2D array of the height of the radar point above the radar elevation. 
%     StationElevationInMeters: A scalar containing the radar site elevation above sea level (m MSL).
%     StationLatitude: a scalar containing the radar site's latitude
%     StationLongitude: a scalar containing the radar site's longitude

%
% OUTPUTS:
% (1) *firename*auto_perimeters_simp.mat : this contains the structure array "fire_cloud", which stores the fire perimeter information for each time step
% 
%     fire_cloud(aa).time: time stamp for current data in MATLAB serial time (fractional days since 1/1/1970) and in UTC
%     fire_cloud(aa).x: current or "active" x fire points
%     fire_cloud(aa).y: current or "active" y fire points
%     fire_cloud(aa).xall: current and all previous x-points-> this is the "point cloud"
%     fire_cloud(aa).yall: current and all previous y-points-> this is the "point cloud"
%     fire_cloud(aa).xspot: current x-spot fire points
%     fire_cloud(aa).yspot: current y-spot fire points
%     fire_cloud(aa).xback: current back edge 
%     fire_cloud(aa).yback: current back edge
%     fire_cloud(aa).xperim: smoothed perimeter polygon x-points fit to xall
%     fire_cloud(aa).yperim: smooth perimeter polygon y-points fit to yall

% (2) *firename*yyyymmddHHMM_progression.png : optional output of figures for each time step (yyyymmddHHMM)
%
%
% REQUIRED CODES AND PACKAGES: 
% (1) Image Processing Tool Box
% (https://www.mathworks.com/products/image.html), which contains some
% functions used including "medfilt2"
% (2) Signal Processing Tool Box
% (https://www.mathworks.com/help/signal/index.html?searchHighlight=signal%20processing%20toolbox&s_tid=srchtitle_signal%20processing%20toolbox_1)
% which includes "findpeaks"
% (3) Fit-Ellipse: Ohad Gal (2021). fit_ellipse
% (https://www.mathworks.com/matlabcentral/fileexchange/3215-fit_ellipse),
%  MATLAB Central File Exchange. Retrieved October 27, 2021.

%% 
close all
clear all

%% Add paths (will need to be changed based on your machine and the data directories)
addpath('~/Dropbox/MATLAB/');
addpath('~/Dropbox/CAMP_FIRE_RADAR/');
clrs=jet(400);

%% Select a Fire Case (camp, bear, king, or caldor)

prompt = 'Enter Fire Name: bear or camp (all lower case):';
str = input(prompt,'s');
disp(str)
if isempty(str)
    str = 'camp';
end
cases={str};

%% Start by loading the gridded radar reflecitivity data. This is the  reflectivity interpolated to an x,y,z grid defined by xvec, yvec,zvec.
if strcmp(cases,'camp')
    [file,path]=uigetfile();
    load(file,'dbZmaster','timemaster','xmaster','ymaster','zmaster','StationElevationInMeters','StationLatitude','StationLongitude');
    startf=18; %this is the starting time step for the analysis and is different for each case since the radar data may precede the fire ignition
    mdx=10^4; %initial search maximum distance
    xlimer=[-1.5 2.5]*10^4; %xlimits for plotting domain (in meters from radar)
    ylimer=[1.5 4.5]*10^4; %ylimits for plotting domain
    fout=strcat('camp_fire_auto_perimeters_',datestr(datenum(datetime('now')),'yyyymmdd'),'.mat'); %output file
    ofig='camp_progression'; %output figure names
    
elseif strcmp(cases,'bear')
    [file,path]=uigetfile();
    load(file,'dbZmaster','timemaster','xmaster','ymaster','zmaster','StationElevationInMeters','StationLatitude','StationLongitude');
    startf=19; %this is the starting time step for the analysis and is different for each case since the radar data may precede the fire ignition
    mdx=3*10^4;%initial search maximum distance
    xlimer=[0 6]*10^4;
    ylimer=[.5 4.2]*10^4;
    fout=strcat('bear_fire_auto_perimeters_',datestr(datenum(datetime('now')),'yyyymmdd'),'.mat'); %output file
    ofig='bear_progression';
    
end


%% Figure Flag
prompt = 'Display Figures?: yes or no (lower case)'; %this will tell the code which fire you're processing
str = input(prompt,'s');
if strcmp(str,'yes')
    pflag=1; %plot figures in code
else
    pflag=0; %don't plot figures in code
end
%%
close all; %make sure all figure windows are closed
[bb,cc]=butter(5,1/3); %spatial filtering coefficients (butterworth filter, 5th order, 1/n, where n= number of points in the filter
xall=[];yall=[]; %initialize x,y arrays for storing perimeter points

aa=0; %counter variable
mx1=0;mx2=0; %max search distances (may no longer be used)
%% Loop through the radar times and compute fire perimters
for tt=startf:2:size(dbZmaster,1)-1% for each time step of the radar data (moves in strides of 2 because we take the maximum for each r,azimuth point between the two PRF sweeps
    aa=aa+1; %log the number of increments
    xer=[];yer=[];xer2=[];yer2=[]; xback=[];yback=[]; %create some empy x and y point matrices that we'll push the perimter points into
    xspot=[]; yspot=[];%create empty arrays for potential spot fires
    zmax1=squeeze(nanmax(dbZmaster(tt-1:tt,:,:),[],1)); % define the column maximum reflectivity for each x, y point in the domain
    zmax1(zmax1<5)=nan; %remove non-smoke related returns (dbz threshold of 5 here is somewhat arbitrary, but based on experience. Certainly returns linked to active combustion tend to be much larger (e.g.>20 dbz)
    zmax1(isnan(zmax1))=0; %nan points = bad data but set to zero for smooth edges on radar returns
    zmax1=medfilt2(zmax1,[3 3]); % 3x3 median filter of the radar data to favor smooth results (range resolution is ~250 m, 0.5 deg azimutha resolution yields variable spatial scale as a function of distance)
    xmat=squeeze(xmaster(tt,:,:)); %parse out current sweep x-data
    ymat=squeeze(ymaster(tt,:,:)); %parse out currnet sweep y-data
    
    %% Create a gridded interpolant of the  maximum reflectivity. This interpolant will be used to identify the peaks along search radials later in the code
    tic
    F1=TriScatteredInterp(xmat(zmax1>5),ymat(zmax1>5),zmax1(zmax1>5)); %F1 is essentially a surface fit to the data which can be queeried for an x, y point. Restrict the surface to dbZ values above 5. 
    toc
    % to reduce the computational load we only use returns above 5 dbZ
    %% Plot the current time step radar data
    f=figure(11);clf; set(f,'pos',[100 100 800 800],'color','w');
    cmap=[[1 1 1]; jet(24)];
    colormap(cmap);
    %mask data outside the region of interest for computational efficiency and to avoid artificats far from the fire
    if aa>1
        zmax1(inner11)=nan;
    end
    pcolor(xmat,ymat,(zmax1));shading flat; %pcolor mesh display
    caxis([0 45]); 
    xlim(xlimer);% the x limits are provided above for each case
    ylim(ylimer);% the y limits are provided above for each case
    daspect([1 1 1]); %make sure it is data aspect 1 space (that is equal x and y spacing)
    hold on;
    grid on; box on; set(gca,'layer','top');
    title(datestr(timemaster(tt),'dd-mmm HH:MM'),'fontsize',15,'fontweight','bold'); %add time stamp as the title
    hold on;

    %% exclude any areas from the analysis that have known issues/terrain clutter
    if aa==1 %only for first time step of the data
        prompt = 'remove supurious region?: yes or no (lower case)'; %this will tell the code which fire you're processing
        str = input(prompt,'s');
        if strcmp(str,'yes')
            [xo1,yo1]=ginput();
            plot(xo1,yo1,'--k');
        end
    end
    if exist('xo1','var')
        inner11=inpolygon(xmat,ymat,xo1,yo1);
        zmax1(inner11)=nan;
    end
    %% Create polygons and search points
    if aa<3 %if it is the first timestep input a manually determined fire perimeter
        prompt = 'Use mouse to draw initial fire polygon estimate. Must include >=20 points. Preses "enter" to continue'; %this will tell the code which fire you're processing
        str = input(prompt,'s');
        [xer,yer]=ginput(); %this takes x,y inputs from the screen. Press "Enter" when you're done drawing the shape
        if numel(xer)>=20
            plot(xer,yer,'--k','linewidth',2)
        else
            prompt = 'Not enough points, try again: Use mouse to draw initial fire polygon estimate. Must include >=20 points. Preses "enter" to continue'; %this will tell the code which fire you're processing
            str = input(prompt,'s');
            [xer,yer]=ginput(); %this takes x,y inputs from the screen. Press "Enter" when you're done drawing the shape
        end
        
    elseif aa>=3 %if we're past the manual steps
        
        % loop over the polygon center, 1/4 and 3/4 search points which are contained in xc,yc 
        for it=1:(length(xc))% 
            x=xc(it);% first (through ith) x center;
            y=yc(it);% first (through ith) y center;
            azset=[0:.5:360];% azimuths in 0.5 degree increments for radial search
            xold=fire_cloud(aa-1).xperim; %read in the x-perimeter from the previous time step
            yold=fire_cloud(aa-1).yperim; %read in the y-perimeter from the previous time step
            maxer=min(sd(it),mdx);% set maximum search radius to be the minium of sd and mdx
            if it<length(xc) 
                z75=30;% this is the dbz_x threshold to be used in find peaks below
            elseif it==length(xc)
                z75=35; %use a slightly higher dbz_x for the final search point, this is used to focus on spot fires
                maxer=20*10^3; 
            end
            
            rl=[0:100:maxer]; %radial search distance vector (in 100 m increments) from the central x,y point
            for az=azset% loop over 360 degrees of azimuth in 0.5 degree increments
                xdist=cosd(az).*rl; %xdist along the radial
                ydist=sind(az).*rl; %ydist along the radial
                y2=y+ydist; %actual y points in cartesian space (i.e., added to the search center location)
                x2=x+xdist; %actual x points in cartesian space (i.e., added to the search center location)
                
                %find the part of the radial that is inside the old perimeter and trim the remaining part to 10 km (10*10^3) beyond the perimeter
                inner=inpolygon(x2,y2,xold,yold); %this locates the x2,y2 points falling inside of the perimeter defined by xold and yold
                dister=sqrt(xdist.^2 + ydist.^2); % along line distance
                dister(inner)=nan; %set distances to nan inside the polygon
                %find first point nearest edge of the polygon
                didx=find(isnan(dister)==0,1,'first'); %location of first non-nan data point (i.e., outside the boundary)
                %since we move in 100 m steps along the radial, 10 km from this point is 100 points
                if it<=3 %only trim this for the first two iterations, and allow the third sweep to search further
                    x2(didx+100:end)=[]; %remove search points far from the fire (i.e., set anything from the edge+100 points out to the end of the vector to NaN)
                    y2(didx+100:end)=[]; %remove search points far from the fire
                else %we can trim the third search radial too, but might opt to do so differently (commented out for now)
%                     x2(didx+200:end)=[]; %remove search points far from the fire
%                     y2(didx+200:end)=[];
                end
                
                
                hold on;
                if isreal(x2) %provided x2 contains real values (I think something got screwy onces, so added this as a check)
                    dbzline=F1(x2,y2); %interpolage the radar data (F1) to the x2,y2 points (this is the reflectivity along that azimuth)
                    dbzline(isnan(dbzline))=-5; %set nan values to -5 dbZ, to help in the find peaks routine
                    [bb,cc]=butter(5,1/5); %Butterworth filtering coefficients
                    dbzline=filtfilt(bb,cc,dbzline); %bidirectional filter to smooth line but perserve peak locations ()
                    dbzline(dbzline<0)=nan; %after smooth set any values <0 dbz to be nan (since we added these in before smoothing since filtfilt chokes on NaNs)
                    if numel(dbzline)>5 && nanmax(dbzline)>z75 %provide the line still has more than 5 data points and the peak is above our max dbz threshold
                        [pks,locs,widths,proms]=findpeaks(dbzline,'MinPeakHeight',z75,'MinPeakProminence',5,'MinPeakDistance',5,'Npeaks',2,'SortStr','descend'); %this uses Matlabs find peaks function (see text in manuscript for discussion)
                        dbznew=dbzline; %make a copy of the radar data
                                            
                        idx=[]; %blank the idx variable
                        if numel(locs)>0 %if there is 1 or more peak found by find peaks
                            if rl(locs(1))<max(rl) && it<length(xc) %provided the max is not the end point of the interpolated dbz line
                                idx1=locs(1); %set our "fire point" to be the first peak
                                nansummer=zeros(length(dbzline)); %set up a nan mask by creating a vector of zeros
                                nansummer(isnan(dbzline))=1; %set any value along the vector to 1 if it is NaN in the interpolated dbZ data
                                nancount=nansum(nansummer(1:idx1))./idx1; %this is the fraction of the data points that are NaN
                                if nancount<.1 %provide less than 10% of the data are NaN
                 
                                    if idx1 %if an index point was found
                                        xer=cat(1,xer,x2(idx1)); %add this point to a current list of active fire points (x points)
                                        yer=cat(1,yer,y2(idx1)); %add this point to a current list of active fire points (y points_
                                        
                                        %find likely back edge of combustion and short-range spotting zone
                                        dbznew2=dbznew; %create a copy of the interpolated dbz values
                                        dbznew2(1:locs(1))=nan; %remove any data past the maximum value (which in this case is the interior of the fire based on the way the interpolation and line segment is done).
                                        if pks(1)>40 %since deep flaming and spotting are likely linked to high dbz, only perform for regions were dbzmax>35 dbz
                                            idx2=find(dbznew2<(.9*pks(1)),1,'first'); %find where the value first drops below 90% of the peak value (e.g. if the peak is 45 dbz find the first point below 40.5 dbz)
                                            if numel(idx2)>0
                                                xback=cat(1,xback,x2(idx2)); %add these points to a "back edge"
                                                yback=cat(1,yback,y2(idx2)); %add these points to a "back edge" ->aka leading edge
                                            end
                                        end
                                        
                                    end
                                end
                            end
                            if it==length(xc) && numel(locs)>=2 && max(pks(2:end))>40 %provided there are at least two peaks and that the second peak has a large value (>35 dbz)
                                idx=locs(2:end);
                                xspot=cat(1,xspot,x2(idx)');%add these points to the spot fire potential list
                                yspot=cat(1,yspot,y2(idx)');
                            end 
                        end
                    end
                end
            end
        end
        
        %% remove any outlier points by examining the distribution of distances from each point to all other points in the point cloud
        for iii=1:length(xer) %loop over all the x,y points we've detected
            xdist=xer(iii)-[xer(:); xall(:)]; %distance in x, of the current point to all other points including previous points (stored in xall)
            ydist=yer(iii)-[yer(:); yall(:)]; %distance in y, of the current point to all other points including previous points (stored in yall)
            ddist=sqrt(xdist.^2 + ydist.^2); %total Euclidean Distance
            ddist(ddist==0)=nan; %replace zers distances with NaN (the zero point is the distance to itself and thus not of interest)
            nclose=nansum(ddist<.1*10^4); %this is the total number of points within 10 km
            if nclose<17 || min(ddist>7000)%if there are less than 10 proximal points or if the point is more than 5 km from the all other points
                %flag point as an outlier
                if pflag==1
                    plot(xer(iii),yer(iii),'*c');
                end
                %if it is outside the polygon add to a potential spot
                %fire list of pixels
                xspot=cat(1,xspot,xer(iii)); %add to spot fire list
                yspot=cat(1,yspot,yer(iii));
                xer(iii)=nan; %remove from perimeter list
                yer(iii)=nan; %remove from perimeter list
            end
            
        end
        
        xer(isnan(xer))=[]; %toss all nan x points
        yer(isnan(yer))=[]; %toss all nan y points
    end
    
    
    hold on;
    scatter(xer,yer,5,'*k');
    %add data to the master point cloud (xall, yall)
    xall=cat(1,xall,xer);
    yall=cat(1,yall,yer);
    
    %% add data points to the current point cloud
    
    %now use the MATLAB boundary function to determine the boundary of
    %the total point cloud. This includes all previous points to date,
    %so that it retains inactive older perimeter points and the fire perimeter thus grows in time to encompass all current and previous points
    k = boundary((xall),((yall)),.5); %the fit parameter (0.5 here) can be tuned. 0= convex hull, 1 = tightest possible polygon
    [bb,cc]=butter(5,1/3);%use a smoothing factor
    %smooth polygon
    xv=filtfilt(bb,cc,xall(k)); %smooth x
    yv=filtfilt(bb,cc,yall(k)); %smooth y
    pnow=polyshape(xv,yv,'Simplify',true); %convert to a polyshape 
    
    %% NOW WE STORE THE VALUES IN "fire_cloud" which is a structure array containg all the relevant fire data for a give time step (aa). 
    fire_cloud(aa).x=xer; %current or "active" x fire points
    fire_cloud(aa).y=yer; %current or "active" y fire points
    fire_cloud(aa).time=timemaster(tt); %time stamp for current data
    fire_cloud(aa).xall=xall; %current and all previous x-points-> this is the "point cloud"
    fire_cloud(aa).yall=yall; %current and all previous y-points-> this is the "point cloud"
    fire_cloud(aa).xspot=xspot; %current x-spot fire points
    fire_cloud(aa).yspot=yspot; %current y-spot fire points
    fire_cloud(aa).xback=xback; %current back edge 
    fire_cloud(aa).yback=yback; %current back edge
    fire_cloud(aa).xperim=xv; %smoothed perimeter polygon x-points fit to xall
    fire_cloud(aa).yperim=yv; %smooth perimeter polygon y-points fit to yall
    
    
    if pflag==1
        f=figure(11); hold on;
        plot(xall,yall,'*','markersize',1,'color','k');%clrs(tt,:));
        hold on;
        scatter(xer,yer,20,'*k');
        scatter(xspot,yspot,'sc');
        scatter(xback,yback,'*m');
        hold on;
        plot(xv,yv,'--k','linewidth',3);
        pause(.1)
    end
    
    %store xperim and yperim in temporary variables (xnow,ynow)
    xnow=fire_cloud(aa).xperim;
    ynow=fire_cloud(aa).yperim;
    
        
    %% FIT AN ELIPSE TO THE CURRENT PERIMETER
    test=fit_ellipse(xnow,ynow,gca); %this uses "fit_ellipse" available here: Ohad Gal (2021). fit_ellipse (https://www.mathworks.com/matlabcentral/fileexchange/3215-fit_ellipse), MATLAB Central File Exchange. Retrieved October 26, 2021.
    
    % use values returned from fit ellipse to determine the major axis
    R = [ cos(test.phi) sin(test.phi); -sin(test.phi) cos(test.phi)];
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = test.X0 + test.a*cos( theta_r );
    ellipse_y_r     = test.Y0 + test.b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];    
    edist=sqrt((rotated_ellipse(1,:)-test.X0_in).^2 + (rotated_ellipse(2,:)-test.Y0_in).^2);
    X0=test.X0;Y0=test.Y0;
    b=test.b;
    a=test.a;
    ver_line        = [ [X0 X0]; Y0+b*[-1 1] ];
    horz_line       = [ X0+a*[-1 1]; [Y0 Y0] ];
    new_ver_line    = R*ver_line;
    new_horz_line   = R*horz_line;
    
    
    %determine which is the major and minor axis between "horz" and "vert" lines.
    dx1=new_horz_line(1,2)-new_horz_line(1,1);
    dy1=new_horz_line(2,1)-new_horz_line(2,2);
    dister1=sqrt(dx1.^2 + dy1.^2); %length of "horz" line
    
    
    dx2=new_ver_line(1,2)-new_ver_line(1,1);
    dy2=new_ver_line(2,1)-new_ver_line(2,2);
    dister2=sqrt(dx2.^2 + dy2.^2); %length of "vert" line
    
    if dister2>dister1 %if the "Vertical" line is longer set it as the major axis
        major=new_ver_line;
        dx=dx2;
        dy=dy2;
    else %if the "horizontal line" is longer set it as the major axis. 
        major=new_horz_line;
        dx=dx1;
        dy=dy1;
    end
    
    %% CONSTRUCT AND TRIM MAJOR AXIS VECTOR
 
    % Now construct a line representing the major axis
    xp=major(1,:);xpi=linspace(xp(1),xp(2),100); 
    yp=major(2,:);ypi=linspace(yp(1),yp(2),100);
    
    % next we need to trim this line to only include points inside the current fire perimeter (stored as xnow,ynow)
    inner=inpolygon(xpi,ypi,xnow,ynow); %find points along the major axis inside the fire perimeter
    xpi=xpi(inner); %remove x points outside of the polygon
    ypi=ypi(inner); %remove y points outside of the polygon
   
    plot(major(1,:),major(2,:),'--k','linewidth',2); %plot the major axis on the figure
 
    
    %% LOCATE THREE SEARCH POINTS ALONG THE MAJOR AXIS
    
    sd=ones(3,1).*(1.5*nanmax(dister1,dister2)); %search distances set to be 1.5 times the length of the major axis
    lx=length(xpi);
    third=round(lx./3); %this is the 1/3 length of the major axis (rounded for integer value)
    twothird=2*third; %this is the 2/3 length of the major axis
    halfer=round(lx./2); % this is the center point (rounded for integer value)
    
    xc=[xpi(third) xpi(twothird) xpi(halfer)]; % these are the three search point xcenters
    yc=[ypi(third) ypi(twothird) ypi(halfer)]; % these are the three search point ycenters
    
    %note that the xc,yc values are stored from one time step to the next, such that they will be used above in the radial search
    %%
    if aa>3 %if we're past the manual estimation of the perimter
        xold=fire_cloud(aa-1).xperim; %select the previous time steps perimeter as "xold"
        yold=fire_cloud(aa-1).yperim; %select the previous time steps perimeter as "yold"
        inold=inpolygon(xc,yc,xold,yold);%check to see which polygon centers are in the original polygon
        xc=xc(inold==1); %only keep center points inside the old perimeter
        yc=yc(inold==1); %only keep center points inside the old perimter
        %note this step simply ensures that if there is an unphysical change in the fire size for some reason that the search remains centered within the previous polygon
    end
    
    
    if pflag==1
        plot(xc,yc,'sm','markersize',10,'linewidth',2,'markerfacecolor','m');
        grid on; box on;
        set(gca,'layer','top','linewidth',2,'fontsize',20,'fontweight','bold');
        xlabel('Distance [m]');
        ylabel('Distance [m]');
        cbh=colorbar;
        ylabel(cbh,'dbZ','fontsize',20,'fontweight','bold');
    end           
end


%% Save the perimeter data
save(fout,'fire_cloud');

