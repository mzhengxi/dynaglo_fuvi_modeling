% RUN THIS EVERY TIME
clear, close all, clc

filename = "\\lasp-store\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\gold\2018_286\00002\0256\outputs\DAY\dayglow_scene_00002_256_v1_r0.sav";
IDLVar=restore_idl_modach(convertStringsToChars(filename));
DataCube = double(IDLVar.SPECT);
WL = double(IDLVar.LAM);
Ch =[1,2,3];
Int = 1.5*ones(size(Ch));
binnum=8;

% Generate Images
PhotoEvents = DYNAGLO_ImageSIM2(DataCube,WL,Ch,Int,binnum);
%%
%Add in FOV correction 
FOV = linspace(-20, 20, length(PhotoEvents(:,:,1)));
[FX,FY]=meshgrid(FOV,FOV);
fov_over_image = (sqrt(FX.^2+FY.^2));
fov_correction = cosd(fov_over_image);
allimages_corrected = PhotoEvents./fov_correction;

allimages_corrected = allimages_corrected./(mean(mean(allimages_corrected)));
%%
figure
channel_num = 3;
s = surf(allimages_corrected(:,:,channel_num));
s.EdgeColor = 'none';
title('GOLD Image Ch' + string(channel_num))
view(0,90)
axis square
map = [0 0 1
    0 .5 1
    0 1 1
    0 1 .5
    0 1 0];
colormap(map)
colorbar