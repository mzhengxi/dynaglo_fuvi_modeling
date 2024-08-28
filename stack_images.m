% RUN THIS EVERY TIME
clear, close all, clc

% Load the different file names
%all_scenenums = ["00","01", "02","03", "04", "05", "06", "07", "08","09", "10"];
all_scenenums = ["00500","00501","00502","00503","00504","00505","00506","00507","00508","00509","00510","00511","00512","00513","00514","00515","00516","00517"]; 
%all_scenenums = ["00500"];

% Prepare the arrays for the image outputs and for the lat and lon for each
% image
allch1images = nan(256, 256, length(all_scenenums));
allch2images = nan(256, 256, length(all_scenenums));
allch3images = nan(256, 256, length(all_scenenums));
alllat = nan(256,256,length(all_scenenums));
alllon= nan(256,256,length(all_scenenums));
%% Loop through all scenes

% DON'T RUN THIS IF YOU HAVE THE IMAGES ALREADY SAVED

for i = 1:length(all_scenenums)
    close all
    %clearvars -except all_scenenums i
    scenenum = all_scenenums(i);
    
    % Set parameters
    filename = "\\lasp-store\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00500\0256\outputs\DAY\dayglow_scene_"+scenenum+"_0256_v1_r0.detrend_on2.sav";
    % filename = "\\lasp-store.lasp.colorado.edu\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\outputs\DAY\dayglow_scene_"+ scenenum +"_256_v1_r0.detrend_ampx100_on2.sav";
    IDLVar=restore_idl_modach(convertStringsToChars(filename));
    SPECT = IDLVar.SPECT;
    WL = IDLVar.LAM;
    Ch =[1,2,3];
    Int = 1.5*ones(size(Ch));
    binnum=8;
    
    % SZA correction, save lat and lon
    params = "\\lasp-store\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00500\0256\inputs\angles\geolocated_scene_params_"+scenenum+"_0256.sav";
    %params = "\\lasp-store\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\inputs\angles\geolocated_scene_params_"+scenenum+"_256.sav";
    IDLVar=restore_idl_modach(convertStringsToChars(params));
    lat = IDLVar.LAT;
    lon = IDLVar.LON;
    sza = IDLVar.SZA;
    DataCube = zeros(size(SPECT));
    for j = 1:size(SPECT,1)
        DataCube(j,:,:) = squeeze(SPECT(j,:,:))./cosd(sza);
    end

    % Generate Images
    PhotoEvents = DYNAGLO_ImageSIM2(DataCube,WL,Ch,Int,binnum);

    % Save images
    allch1images(:,:,i) = PhotoEvents(:,:,1);
    allch2images(:,:,i) = PhotoEvents(:,:,2);
    allch3images(:,:,i) = PhotoEvents(:,:,3);
    alllat(:,:,i) = reshape(lat,[256, 256, 1]);
    alllon(:,:,i) = reshape(lon,[256, 256, 1]);
end

% Save as a separate file so that it does not need to be run every time
save("fuvi_runs_081324.mat","all_scenenums", "allch1images", "allch2images", "allch3images", "alllat", "alllon");

%% Set up a common lat lon grid

% Load saved images and lat lon
load("fuvi_runs_081324.mat");

% Set up final lat lon for composite image interpolation later, 0.1 degree
% increments
commonlat = round(min(min(min(alllat))),1):0.1:round(max(max(max(alllat))),1);
commonlon = round(min(min(min(alllon))),1):0.1:round(max(max(max(alllon))),1);
[finallon, finallat] = meshgrid(commonlon, commonlat);
%% Correct and subtract images

%Add in FOV correction 
FOV = linspace(-20, 20, length(allch1images(:,:,1)));
[FX,FY]=meshgrid(FOV,FOV);
fov_over_image = (sqrt(FX.^2+FY.^2));
fov_correction = cosd(fov_over_image);
allch1images_corrected = allch1images./fov_correction;
allch2images_corrected = allch2images./fov_correction;
allch3images_corrected = allch3images./fov_correction;
ch1m2images_corrected = (allch1images-allch2images)./fov_correction;

%ch1mch2 = allch1images_corrected-allch2images_corrected;
%ch2mch3 = allch2images_corrected-allch3images_corrected;

% Take each image and divide by its mean
allch1images_corrected = allch1images_corrected(:,:,:)./mean(mean(allch1images_corrected(:,:,:)));
allch2images_corrected = allch2images_corrected(:,:,:)./mean(mean(allch2images_corrected(:,:,:)));
allch3images_corrected = allch3images_corrected(:,:,:)./mean(mean(allch3images_corrected(:,:,:)));
%ch1m2images_corrected = ch1m2images_corrected(:,:,:)./mean(mean(ch1m2images_corrected(:,:,:)));
%% Plot CONOPS
%ch1_vars = zeros(256, 256, 11);
%for k = 1:length(all_scenenums)
%    ch1_vars = 
%end
%ch2_vars = zeros(256, 256, 11);
%
figure 

for k = 1:length(all_scenenums)
    if mod(k,2)==1
        % For every other image, plot Ch1
        s = surf(alllat(:,:,k),alllon(:,:,k)-8*(k-1)*ones(256, 256, 1), normalize(allch1images_corrected(:,:,k), 'scale'));
        s.EdgeColor = 'none';
        view(0,90)
        text(alllat(1,1,k),alllon(100,100,k)-8*(k-1),"CH1 S" + all_scenenums(k))
        
    else
        % For every other image, plot Ch2
        s = surf(alllat(:,:,k), alllon(:,:,k)-8*(k-1)*ones(256, 256, 1), allch2images_corrected(:,:,k));
        s.EdgeColor = 'none';
        view(0,90)
        text(alllat(1,1,k),alllon(100,100,k)-8*(k-1),"CH2 S" + all_scenenums(k))
        
    end
    %axis square
    %pbaspect([2 1 1])
    hold all
    %s = surf(alllat(:,:,2),alllon(:,:,2)-8*ones(256, 256, 1), allch1images_corrected(:,:,2));
    %s.EdgeColor = 'none';
end
title('Concept of Operations')
set(gca,'ytick',[]);
xlabel('Latitude')
%test=min(min(allch1images_corrected(:,:,:)./mean(allch1images_corrected(:,:,:))));
%test2=max(max(allch2images_corrected(:,:,:)./mean(allch2images_corrected(:,:,:))));


%clim([.9 50*.006+.9])
%map = [0 0 1
%    0 1 1
%    0 1 0];
%colormap(map)
axis square

%% Plot all of a channel 1 overlapped
%ch1_vars = zeros(256, 256, 11);
%for k = 1:length(all_scenenums)
%    ch1_vars = 
%end
%ch2_vars = zeros(256, 256, 11);
%
figure 
% Every other image
for k = 1:2:length(all_scenenums)
    s = surf(alllat(:,:,k),alllon(:,:,k), allch1images_corrected(:,:,k));
    s.EdgeColor = 'none';
    view(0,90)
    %text(alllat(1,1,k),alllon(100,100,k)-8*(k-1),"CH1 S" + all_scenenums(k))
    %axis square
    %pbaspect([2 1 1])
    hold all
    %s = surf(alllat(:,:,2),alllon(:,:,2)-8*ones(256, 256, 1), allch1images_corrected(:,:,2));
    %s.EdgeColor = 'none';
end
title('All Channel 1 Images')
%set(gca,'ytick',[]);
xlabel('Latitude')
ylabel('Longitude')
%test=min(min(allch1images_corrected(:,:,:)./mean(allch1images_corrected(:,:,:))));
%test2=max(max(allch2images_corrected(:,:,:)./mean(allch2images_corrected(:,:,:))));


clim([.9 50*.006+.9])
map = [0 0 1
    0 1 1
    0 1 0];
colormap(map)
axis square



%% Plot all of a channel 2 overlapped
%ch1_vars = zeros(256, 256, 11);
%for k = 1:length(all_scenenums)
%    ch1_vars = 
%end
%ch2_vars = zeros(256, 256, 11);
%
figure 
% Every other image starting from 2
for k = 2:2:length(all_scenenums)
    s = surf(alllat(:,:,k),alllon(:,:,k), normalize(allch2images_corrected(:,:,k),'zscore'));
    s.EdgeColor = 'none';
    view(0,90)
    %text(alllat(1,1,k),alllon(100,100,k)-8*(k-1),"CH1 S" + all_scenenums(k))
    %axis square
    %pbaspect([2 1 1])
    hold all
    %s = surf(alllat(:,:,2),alllon(:,:,2)-8*ones(256, 256, 1), allch1images_corrected(:,:,2));
    %s.EdgeColor = 'none';
end
title('All Channel 2 Images')
%set(gca,'ytick',[]);
xlabel('Latitude')
ylabel('Longitude')
%test=min(min(allch1images_corrected(:,:,:)./mean(allch1images_corrected(:,:,:))));
%test2=max(max(allch2images_corrected(:,:,:)./mean(allch2images_corrected(:,:,:))));
%clim([.9 50*.006+.9])
%map = [0 0 1
%    0 0.25 1
%    0 0.5 1
%    0 0.75 1
%    0 1 1
%    0 1 0.75 
%    0 1 0.5
%    0 1 0.25
%    0 1 0];
%colormap(map)
axis square

%% Plot all of a channel 3 overlapped
%ch1_vars = zeros(256, 256, 11);
%for k = 1:length(all_scenenums)
%    ch1_vars = 
%end
%ch2_vars = zeros(256, 256, 11);
%
figure 

for k = 1:length(all_scenenums)
    s = surf(alllat(:,:,k),alllon(:,:,k), allch3images_corrected(:,:,k));
    s.EdgeColor = 'none';
    view(0,90)
    %text(alllat(1,1,k),alllon(100,100,k)-8*(k-1),"CH1 S" + all_scenenums(k))
    %axis square
    %pbaspect([2 1 1])
    hold all

end
title('All Channel 3 Images')
%set(gca,'ytick',[]);
xlabel('Latitude')
ylabel('Longitude')
%test=min(min(allch1images_corrected(:,:,:)./mean(allch1images_corrected(:,:,:))));
%test2=max(max(allch2images_corrected(:,:,:)./mean(allch2images_corrected(:,:,:))));
%clim([.9 50*.006+.9])
%map = [0 0 1
%    0 1 1
%    0 1 0];
%colormap(map)
axis square

%% Plot one image
figure
s = surf(alllat(:,:,1),alllon(:,:,1), allch3images_corrected(:,:,1));
s.EdgeColor = 'none';
title('Ch 3 Image 00')
view(0,90)
axis square
%clim([.9 50*.006+.9])
%map = [0 0 1
%    0 1 1
%    0 1 0];
%colormap(map)
%% Create interpolated images on common lat lon grid Channel 1

int_ch1 = zeros(length(finallat), size(finallon,2), length(all_scenenums));
%figure
for i = 1:length(all_scenenums)
    %subplot(11,1,i);
    %Vq = interp2(alllon(:,:,1), alllat(:,:,1), ch1mch2(:,:,1), finallon, finallat);
    image = allch1images_corrected(:,:,i);
    imagelon = alllon(:,:,i);
    imagelat = alllat(:,:,i);
    F = scatteredInterpolant(imagelat(:),imagelon(:), image(:),'linear','none');
    interp_image = F(finallat, finallon);
    
    
    %set extrapolated values to NaN
    DT = delaunayTriangulation(imagelat(:),imagelon(:));
    b = boundary(imagelat(:),imagelon(:));
    if ~inpolygon(finallat, finallon, DT.Points(b,1),DT.Points(b,2))
        interp_image = NaN;
    end

    int_ch1(:,:,i) = interp_image;
    
    %xq = griddata(imagelon(:),imagelat(:), image(:), finallon, finallat);
    %interp2(imagelon,imagelat,image,finallon,finallat);
    
    %figure
   %if mod(i,2)==1
   %    s = surf(finallat,finallon-8*(i-1)*ones(83, 83, 1), interp_image);
   %    s.EdgeColor = 'none';
   %    %title('Interpolated ' + all_scenenums(i))
   %    hold all
   %    view(0,90)
   %    %axis square
   %    %pbaspect([2 1 1])
   %    text(alllat(1,1,i),alllon(100,100,i)-8*(i-1),"CH1 S" + all_scenenums(i))
   %    set(gca,'ytick',[]);
   %    clim([.9 50*.006+.9])
   %    map = [0 0 1
   %            0 1 1
   %            0 1 0];
   %    colormap(map)
   %end
end

% Total composite
%j = 1;
%image1 = int_ch1(:,:,j);
%image2 = int_ch1(:,:,j+2);
%image3 = int_ch1(:,:,j+4);
%image4 = int_ch1(:,:,j+6);
%image5 = int_ch1(:,:,j+8);
%image6 = int_ch1(:,:,j+10);
%change nans to zeros to add together
%image1(isnan(image1))=0;
%image2(isnan(image2))=0;
%image3(isnan(image3))=0;
%image4(isnan(image1))=0;
%image5(isnan(image2))=0;
%image6(isnan(image3))=0;
%compositeimage1 = max(int_ch1(:,:,1:2:end), [],3, 'omitnan');
compositeimage1 = mean(int_ch1(:,:,1:2:end), 3, 'omitnan');

%set area outside thickest region to zero
%latmax = max(max(alllat(:,:,j)));
%latmin = min(min(alllat(:,:,j+10)));
%lonmin = min(min(alllon(:,:,j)));
%lonmax = max(max(alllon(:,:,j+10)));
%
%
%compositeimage(finallat>latmax) = NaN;
%compositeimage(finallat<latmin) = NaN;
%compositeimage(finallon>lonmax) = NaN;
%compositeimage(finallon<lonmin) = NaN;
%% Plot composite image
figure
s = surf(finallat,finallon, compositeimage1);
s.EdgeColor = 'none';
title('CH1 Composite')
view(0,90)
axis square
xlabel('Latitude')
ylabel('Longitude')
clim([0.8 50*.006+.9])

intervals = (0:49)*.006+.9;
clrs = winter(2*numel(intervals)-2);
colormap(clrs);
colorbar

map = [0 0 1
    0 0.125 1
    0 0.25 1
    0 0.325 1
    0 0.5 1
    0 0.625 1
    0 0.75 1
    0 0.825 1
    0 1 1
    0 1 1
    0 1 0.825 
    0 1 0.75 
    0 1 0.625 
    0 1 0.5
    0 1 0.325 
    0 1 0.25
    0 1 0.125 
    0 1 0];

colormap(map)
%colormap winter
%colorbar

colorbar;
%clim([.9 50*.006+.9])
map = [0 0 1
    0 .5 1
    0 1 1
    0 1 .5
    0 1 0];
%colormap(map)

%% Create interpolated images on common lat lon grid Channel 2

int_ch2 = zeros(length(finallat), size(finallon,2), length(all_scenenums));
figure
for i = 1:length(all_scenenums)
    image = allch2images_corrected(:,:,i);
    imagelon = alllon(:,:,i);
    imagelat = alllat(:,:,i);
    F = scatteredInterpolant(imagelat(:),imagelon(:), image(:),'linear','none');
    interp_image = F(finallat, finallon);
    
    
    %set extrapolated values to NaN
    DT = delaunayTriangulation(imagelat(:),imagelon(:));
    b = boundary(imagelat(:),imagelon(:));
    if ~inpolygon(finallat, finallon, DT.Points(b,1),DT.Points(b,2))
        interp_image = NaN;
    end

    int_ch2(:,:,i) = interp_image;
    
    %figure
    %if mod(i,2)==1
    %    s = surf(finallat,finallon-8*(i-1)*ones(83, 83, 1), interp_image);
    %    s.EdgeColor = 'none';
    %    %title('Interpolated ' + all_scenenums(i))
    %    hold all
    %    view(0,90)
    %    %axis square
    %    %pbaspect([2 1 1])
    %    text(alllat(1,1,i),alllon(100,100,i)-8*(i-1),"CH2 S" + all_scenenums(i))
    %    set(gca,'ytick',[]);
    %    clim([.9 50*.006+.9])
    %    map = [0 0 1
    %            0 1 1
    %            0 1 0];
    %    colormap(map)
    %end
end

% Total composite

%compositeimage2 = max(int_ch2(:,:,1:2:end),[],3, 'omitnan');
compositeimage2 = mean(int_ch2(:,:,2:2:end),3, 'omitnan');

%
figure
s = surf(finallat,finallon, compositeimage2);
s.EdgeColor = 'none';
title('CH2 Composite')
view(0,90)
axis square
xlabel('Latitude')
ylabel('Longitude')
clim([.6 50*.006+.9])
map = [0 0 1
   0 0.25 1
   0 0.5 1
   0 0.75 1
   0 1 1
   0 1 0.75 
   0 1 0.5
   0 1 0.25
   0 1 0];

colormap(map)
%colormap winter
colorbar
%%
%% Create interpolated images on common lat lon grid Channel 3

int_ch3 = zeros(length(finallat), length(finallon), length(all_scenenums));
figure
for i = 1:length(all_scenenums)
    image = allch3images_corrected(:,:,i);
    imagelon = alllon(:,:,i);
    imagelat = alllat(:,:,i);
    F = scatteredInterpolant(imagelat(:),imagelon(:), image(:),'linear','none');
    interp_image = F(finallat, finallon);
    
    
    %set extrapolated values to NaN
    DT = delaunayTriangulation(imagelat(:),imagelon(:));
    b = boundary(imagelat(:),imagelon(:));
    if ~inpolygon(finallat, finallon, DT.Points(b,1),DT.Points(b,2))
        interp_image = NaN;
    end

    int_ch3(:,:,i) = interp_image;
    
    %figure
    if mod(i,2)==1
        s = surf(finallat,finallon-8*(i-1)*ones(83, 83, 1), interp_image);
        s.EdgeColor = 'none';
        %title('Interpolated ' + all_scenenums(i))
        hold all
        view(0,90)
        %axis square
        %pbaspect([2 1 1])
        text(alllat(1,1,i),alllon(100,100,i)-8*(i-1),"CH3 S" + all_scenenums(i))
        set(gca,'ytick',[]);
        clim([.9 50*.006+.9])
        map = [0 0 1
                0 1 1
                0 1 0];
        colormap(map)
    end
end

% Total composite

%compositeimage3 = max(int_ch3(:,:,1:2:end),[],3, 'omitnan');
compositeimage3 = mean(int_ch3(:,:,1:2:end),3, 'omitnan');

%
figure
s = surf(finallat,finallon, compositeimage3);
s.EdgeColor = 'none';
title('CH3 Composite')
view(0,90)
axis square
xlabel('Latitude')
ylabel('Longitude')
clim([.9 50*.006+.9])
map = [0 0 1
    0 0.25 1
    0 0.5 1
    0 0.75 1
    0 1 1
    0 1 1
    0 1 0.75 
    0 1 0.5
    0 1 0.25
    0 1 0];

colormap(map)
%colormap winter
%colorbar

%%

%% Create interpolated images on common lat lon grid Channel 1 m Ch2

int_ch1m2 = zeros(length(finallat), length(finallon), length(all_scenenums));
figure
for i = 1:length(all_scenenums)
    image = ch1m2images_corrected(:,:,i);
    imagelon = alllon(:,:,i);
    imagelat = alllat(:,:,i);
    F = scatteredInterpolant(imagelat(:),imagelon(:), image(:),'linear','none');
    interp_image = F(finallat, finallon);
    
    
    %set extrapolated values to NaN
    DT = delaunayTriangulation(imagelat(:),imagelon(:));
    b = boundary(imagelat(:),imagelon(:));
    if ~inpolygon(finallat, finallon, DT.Points(b,1),DT.Points(b,2))
        interp_image = NaN;
    end

    int_ch1m2(:,:,i) = interp_image;
    
    %figure
    if mod(i,2)==1
        s = surf(finallat,finallon-8*(i-1)*ones(83, 83, 1), interp_image);
        s.EdgeColor = 'none';
        %title('Interpolated ' + all_scenenums(i))
        hold all
        view(0,90)
        %axis square
        %pbaspect([2 1 1])
        text(alllat(1,1,i),alllon(100,100,i)-8*(i-1),"CH3 S" + all_scenenums(i))
        set(gca,'ytick',[]);
        clim([.9 50*.006+.9])
        map = [0 0 1
                0 1 1
                0 1 0];
        colormap(map)
    end
end

% Total composite

%compositeimage3 = max(int_ch3(:,:,1:2:end),[],3, 'omitnan');
compositeimagech1m2 = mean(int_ch1m2(:,:,1:2:end),3, 'omitnan');

%
figure
s = surf(finallat,finallon, compositeimagech1m2);
s.EdgeColor = 'none';
title('CH1-Ch2 Composite')
view(0,90)
axis square
xlabel('Latitude')
ylabel('Longitude')
clim([.9 50*.006+.9])
map = [0 0 1
    0 0.25 1
    0 0.5 1
    0 0.75 1
    0 1 1
    0 1 0.75 
    0 1 0.5
    0 1 0.25
    0 1 0];

colormap(map)
%colormap winter
%colorbar

%% %% Channel subtractions
%
%figure
%s = surf(finallat,finallon, compositeimage1+mean(mean(mean(allch1images./fov_correction)))*ones(83, 83)-(compositeimage2+mean(mean(mean(allch2images./fov_correction)))*ones(83, 83)));
%s.EdgeColor = 'none';
%title('CH1-CH2 Composite')
%view(0,90)
%axis square
%xlabel('Latitude')
%ylabel('Longitude')
%%clim([.9 50*.006+.9])
%map = [0 0 1
%    0 0.1 1
%    0 0.2 1
%    0 0.25 1
%    0 0.3 1
%    0 0.4 1
%    0 0.5 1
%    0 0.75 1
%    0 1 1
%    0 1 1
%    0 1 0.75 
%    0 1 0.5
%    0 1 0.25
%    0 1 0];
%
%colormap(map)
%%colormap winter
%%colorbar
%%% O2 composite
%filename = "\\lasp-store.lasp.colorado.edu\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\outputs\DAY\dayglow_scene_"+ scenenum +"_256_v1_r0.detrend_ampx100_on2.sav";
%IDLVar=restore_idl_modach(convertStringsToChars(filename));

%% Channel subtractions

figure
s = surf(finallat,finallon, compositeimage1+mean(mean(mean(allch1images./fov_correction)))*ones(83, 83)-(compositeimage2+mean(mean(mean(allch2images./fov_correction)))*ones(83, 83)));
s.EdgeColor = 'none';
title('CH1-CH2 Composite')
view(0,90)
axis square
xlabel('Latitude')
ylabel('Longitude')
%clim([.9 50*.006+.9])
map = [0 0 1
    0 0.1 1
    0 0.2 1
    0 0.25 1
    0 0.3 1
    0 0.4 1
    0 0.5 1
    0 0.75 1
    0 1 1
    0 1 1
    0 1 0.75 
    0 1 0.5
    0 1 0.25
    0 1 0];

colormap(map)
%colormap winter
%colorbar
%% Find O and N2 at 150 km from IDL sav file
oden_all150km = zeros(256,256,11);
n2den_all150km = zeros(256,256,11);
on2 = nan(256, 256, length(all_scenenums));
for i = 1:length(all_scenenums)
    scenenum = all_scenenums(i);
    filename = "C:\Users\Maggie\Documents\DYNAGLO\Hoskins Matlab\Files\IDL Sav Files\d"+scenenum+".sav";
    
    %filename = "\\lasp-store\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00500\0256\outputs\DAY\dayglow_scene_"+scenenum+"_0256_v1_r0.detrend_on2.sav";
    IDLVar=restore_idl_modach(convertStringsToChars(filename));
    %oden_all150km(:,:,i) = IDLVar.ODEN(43,:,:)./mean(mean(IDLVar.ODEN(43,:,:)));
    oden_all150km(:,:,i) = IDLVar.ODEN(43,:,:);
    %n2den_all150km(:,:,i) = IDLVar.N2DEN(43,:,:)./mean(mean(IDLVar.N2DEN(43,:,:)));
    n2den_all150km(:,:,i) = IDLVar.N2DEN(43,:,:);
    %on2(:,:,i) = IDLVar.ON2;
end
%%


%% Create composite image for O
int_o2 = zeros(length(finallat), length(finallon), length(all_scenenums));
figure
for i = 1:length(all_scenenums)
    image = oden_all150km(:,:,i);
    imagelon = alllon(:,:,i);
    imagelat = alllat(:,:,i);
    F = scatteredInterpolant(imagelat(:),imagelon(:), image(:),'linear','none');
    interp_image = F(finallat, finallon);
    
    
    %set extrapolated values to NaN
    DT = delaunayTriangulation(imagelat(:),imagelon(:));
    b = boundary(imagelat(:),imagelon(:));
    if ~inpolygon(finallat, finallon, DT.Points(b,1),DT.Points(b,2))
        interp_image = NaN;
    end

    int_o2(:,:,i) = interp_image;
    
    %figure
    if mod(i,2)==1
        s = surf(finallat,finallon-8*(i-1)*ones(83, 83, 1), interp_image);
        s.EdgeColor = 'none';
        %title('Interpolated ' + all_scenenums(i))
        hold all
        view(0,90)
        %axis square
        %pbaspect([2 1 1])
        text(alllat(1,1,i),alllon(100,100,i)-8*(i-1),"O2 S" + all_scenenums(i))
        set(gca,'ytick',[]);
        clim([.9 50*.006+.9])
        map = [0 0 1
                0 1 1
                0 1 0];
        colormap(map)
    end
end

% Total composite
%%
%compositeimage3 = max(int_ch3(:,:,1:2:end),[],3, 'omitnan');
compositeimageo2 = mean(int_o2,3, 'omitnan');

%
figure
s = surf(finallat,finallon, compositeimageo2);
s.EdgeColor = 'none';
title('O Composite')
view(0,90)
axis square
xlabel('Latitude')
ylabel('Longitude')
%clim([0 1.8])
%colorbar
map = [0 0 1
    0 0.25 1
    0 0.5 1
    0 0.75 1
    0 1 1
    0 1 0.75 
    0 1 0.5];

colormap(map)
%%
%figure
%s = surf(finallat,finallon, compositeimage2-compositeimage3);
%s.EdgeColor = 'none';
%title('CH2-3 Composite')
%view(0,90)
%axis square
%xlabel('Latitude')
%ylabel('Longitude')
%%

%%% Create composite images
%composite_ch1 = zeros(length(finallon), length(finallat), 4);
%for j = 1:2:7
%    image1 = int_ch1(:,:,j);%./mean(mean(int_ch1(:,:,j), 'omitnan'), 'omitnan');
%    image2 = int_ch1(:,:,j+2);%./mean(mean(int_ch1(:,:,j+2), 'omitnan'), 'omitnan');
%    image3 = int_ch1(:,:,j+4);%./mean(mean(int_ch1(:,:,j+4), 'omitnan'), 'omitnan');
%
%    %change nans to zeros to add together
%    image1(isnan(image1))=0;
%    image2(isnan(image2))=0;
%    image3(isnan(image3))=0;
%    compositeimage = image1+image2+image3;
%    
%    %set area outside thickest region to zero
%    latmax = max(max(alllat(:,:,j)));
%    latmin = min(min(alllat(:,:,j+4)));
%    lonmin = min(min(alllon(:,:,j)));
%    lonmax = max(max(alllon(:,:,j+4)));
%    compositeimage(finallat>latmax) = NaN;
%    compositeimage(finallat<latmin) = NaN;
%    compositeimage(finallon>lonmax-.3) = NaN;
%    compositeimage(finallon<lonmin+.3) = NaN;
%    composite_ch1(:,:,ceil(j/2)) = compositeimage;
%end
%%% Plot trimmed composites
%titles = ["00 + 02 + 04", "02 + 04 + 06", "04 + 06 + 08", "06 + 08 + 10"];
%figure
%for num = 1:2
%    %figure
%    %subplot(4,1,num);
%    s = surf(finallat,finallon, composite_ch1(:,:,num), 'FaceColor', 'flat', 'FaceAlpha',1, 'FaceLighting','none');
%    s.EdgeColor = 'none';
%    title('Composite CH1 ' + titles(num))
%    view(0,90)
%    axis square
%    test = composite_ch1(:,:,num);
%    clim([mean(mean(test, 'omitnan'),'omitnan')-3*std(std(test, 'omitnan'),'omitnan') mean(mean(test, 'omitnan'),'omitnan')+3*std(std(test, 'omitnan'),'omitnan')])
%    %colorbar
%    %colormap gray
%    
%    hold all
%end
%
%%
%%

%%
%%%
%figure
%s = surf(finallon,finallat, composite_ch2mch3(:,:,2));
%s.EdgeColor = 'none';
%title('Composite ch2mch3')
%view(0,90)
%axis square
%%% trim composite
%
%%%
%test=composite_ch1mch2(:,:,1)./mean(mean(composite_ch1mch2(:,:,1), 'omitnan'), 'omitnan');
%figure
%s = surf(finallon,finallat, test);
%s.EdgeColor = 'none';
%title('Composite ch1mch2 (06+08+10)')
%view(0,90)
%axis square
%colormap gray
%colorbar
%%


%%%
%%%figure
%s = surf(Vq);
%s.EdgeColor = 'none';
%title('Ch 1-Ch 2')
%view(0,90)
%axis square
%
%%%
%%Add in FOV correction 
%FOV = linspace(-20, 20, length(PhotoEvents(:,:,1)));
%[FX,FY]=meshgrid(FOV,FOV);
%fov_over_image = (sqrt(FX.^2+FY.^2));
%fov_correction = cosd(fov_over_image);
%
%figure
%s = surf(PhotoEvents(:,:,1)./fov_correction);
%s.EdgeColor = 'none';
%title('no binning')
%view(0,90)
%axis square
%
%
%%%
%[X,Y,V] = peaks(10); 
%[Xq,Yq] = meshgrid(-3:.1:3,-3:.1:3);
%
%Vq = interp2(X,Y,V,Xq,Yq); 
%mesh(Xq,Yq,Vq)

%% Plot all images 
close all hidden
figure
for k = 1:length(all_scenenums)
    subplot(6,3,k)
    %flat = allch1images_corrected(:,:,1);
    %mean(mean(mean(allch1images_corrected(:,:,k))))*ones(size(allch1images_corrected(:,:,k))))

    smoothed_data= smooth2a(allch2images_corrected(:,:,k),3,3);
    %smoothed_data= smooth2a(normalize(allch1images_corrected(:,:,k),'zscore'),3,3);
    %smoothed_data = n2den_all150km(:,:,k);
    s = surf(alllat(:,:,k),alllon(:,:,k), smoothed_data);
    %s = surf(smoothed_data);
    s.EdgeColor = 'none';
    view(0,90)
    
    colorbar;
    %clim([.9 50*.006+.9])
    map = [0 0 1
        0 .5 1
        0 1 1
        0 1 .5
        0 1 0];
    colormap(map)
    title('Scene '+all_scenenums(k))
    sgtitle("Channel 2/mean smoothed")
    xlim([51 60]);
    ylim([-5 5]);
    clim([.9 1.1])
    %clim([23500 26000]);
    %text(alllat(1,1,k),alllon(100,100,k)-8*(k-1),"CH1 S" + all_scenenums(k))
    %axis square
    %pbaspect([2 1 1])
    %s = surf(alllat(:,:,2),alllon(:,:,2)-8*ones(256, 256, 1), allch1images_corrected(:,:,2));
    %s.EdgeColor = 'none';
    %mean(mean(mean(allch1images_corrected(:,:,k))))
end

%% Functions
function matrixOut = smooth2a(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
%
%function matrixOut = smooth2a(matrixIn,Nr,Nc)
% 
% This function smooths the data in matrixIn using a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
% element "i" by the mean of the rectange centered on "i".  Any NaN
% elements are ignored in the averaging.  If element "i" is a NaN, then it
% will be preserved as NaN in the output.  At the edges of the matrix,
% where you cannot build a full rectangle, as much of the rectangle that
% fits on your matrix is used (similar to the default on Matlab's builtin
% function "smooth").
% 
% "matrixIn": original matrix
% "Nr": number of points used to smooth rows
% "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
% 
% "matrixOut": smoothed version of original matrix
% 
% 
% 	Written by Greg Reeves, March 2009.
% 	Division of Biology
% 	Caltech
% 
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004
% 	Applied Research Laboratory
% 	Penn State University
% 
% 	Developed from code written by Olof Liungman, 1997
% 	Dept. of Oceanography, Earth Sciences Centre
% 	G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se
%
% Initial error statements and definitions
%
if nargin < 2, error('Not enough input arguments!'), end
N(1) = Nr; 
if nargin < 3, N(2) = N(1); else N(2) = Nc; end
if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end
%
% Building matrices that will compute running sums.  The left-matrix, eL,
% smooths along the rows.  The right-matrix, eR, smooths along the
% columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by- 
% (2*Nc+1) rectangle centered on element "i".
%
[row,col] = size(matrixIn);
eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);
%
% Setting all "NaN" elements of "matrixIn" to zero so that these will not
% affect the summation.  (If this isn't done, any sum that includes a NaN
% will also become NaN.)
%
A = isnan(matrixIn);
matrixIn(A) = 0;
%
% For each element, we have to count how many non-NaN elements went into
% the sums.  This is so we can divide by that number to get a mean.  We use
% the same matrices to do this (ie, "eL" and "eR").
%
nrmlize = eL*(~A)*eR;
nrmlize(A) = NaN;
%
% Actually taking the mean.
%
matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;

end
