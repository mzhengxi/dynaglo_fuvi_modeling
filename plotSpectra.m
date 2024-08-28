clear, close all, clc

scenenum = "00500";
filename = "\\lasp-store\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00500\0256\outputs\DAY\dayglow_scene_"+scenenum+"_0256_v1_r0.detrend_on2.sav";
IDLVar=restore_idl_modach(convertStringsToChars(filename));
SPECT = IDLVar.SPECT;
WL = IDLVar.LAM;
params = "\\lasp-store\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00500\0256\inputs\angles\geolocated_scene_params_"+scenenum+"_0256.sav";
IDLVar=restore_idl_modach(convertStringsToChars(params));
sza = IDLVar.SZA;
DataCube = zeros(size(SPECT));
for j = 1:size(SPECT,1)
    DataCube(j,:,:) = squeeze(SPECT(j,:,:))./cosd(sza);
end

%scenenum = "00";
%% Run for amplified file
%filename = "\\lasp-store.lasp.colorado.edu\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\outputs\DAY\dayglow_scene_"+ scenenum +"_256_v1_r0.detrend_ampx100_on2.sav";
%IDLVar=restore_idl_modach(convertStringsToChars(filename));
%% Set Parameters
%DataCube2 = IDLVar.SPECT/cosd(85);
%WL2 = IDLVar.LAM;

% Run for simulated GW
filename = "C:\Users\Maggie\Documents\DYNAGLO\FUVI_simulator_new\dayglow_scene_00012_0512_v1_r0.sav";
IDLVar=restore_idl_modach(convertStringsToChars(filename));
% Set Parameters
DataCube3 = IDLVar.SPEC;
WL3 = IDLVar.LAM;

% Run for flat field
filename = "C:\Users\Maggie\Documents\DYNAGLO\FUVI_simulator_new\dayglow_scene_00003_0512_v1_r0.sav";
IDLVar=restore_idl_modach(convertStringsToChars(filename));
% Set Parameters
DataCube4 = IDLVar.SPEC;
WL4 = IDLVar.LAM;

%%

scenenum = "00";
filename = "\\lasp-store.lasp.colorado.edu\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\outputs\DAY\dayglow_scene_"+ scenenum +"_256_v1_r0.detrend_on2.sav";
IDLVar=restore_idl_modach(convertStringsToChars(filename));
% Set Parameters
DataCube = IDLVar.SPECT/cosd(85);
WL = IDLVar.LAM;

% Run for amplified file
filename = "\\lasp-store.lasp.colorado.edu\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\outputs\DAY\dayglow_scene_"+ scenenum +"_256_v1_r0.detrend_ampx100_on2.sav";
IDLVar=restore_idl_modach(convertStringsToChars(filename));
% Set Parameters
DataCube2 = IDLVar.SPECT/cosd(85);
WL2 = IDLVar.LAM;

% Run for simulated GW
filename = "C:\Users\Maggie\Documents\DYNAGLO\FUVI_simulator_new\dayglow_scene_00012_0512_v1_r0.sav";
IDLVar=restore_idl_modach(convertStringsToChars(filename));
% Set Parameters
DataCube3 = IDLVar.SPEC;
WL3 = IDLVar.LAM;

% Run for flat field
filename = "C:\Users\Maggie\Documents\DYNAGLO\FUVI_simulator_new\dayglow_scene_00003_0512_v1_r0.sav";
IDLVar=restore_idl_modach(convertStringsToChars(filename));
% Set Parameters
DataCube4 = IDLVar.SPEC;
WL4 = IDLVar.LAM;

%%
clf
% Plot on one figure
figure(1)
semilogy(WL./10, DataCube(:,100,100), 'LineWidth',1.5)
hold all
%semilogy(WL2./10, DataCube2(:,100,100), 'r--','LineWidth',1.5)
semilogy(WL3./10, DataCube3(:,100,100), 'k','LineWidth',1.5)
%semilogy(WL3./10, (DataCube3(:,100,100)-DataCube4(:,100,100)), 'g.-','LineWidth',1.5)
legend('HIAMCM (scaled by SZA)','Simulated 100 km wave');
title('Spectra')
ylabel('Radiance R/nm')
xlabel('nm')

hold off
figure(2)
semilogy(WL./10, DataCube(:,100,100), 'LineWidth',1.5)
hold all
%semilogy(WL2./10, DataCube2(:,100,100), 'r--','LineWidth',1.5)
semilogy(WL3./10, DataCube3(:,100,100), 'k','LineWidth',1.5)
%semilogy(WL3./10, abs(DataCube3(:,100,100)-DataCube4(:,100,100)), 'g-','LineWidth',1.5)
%semilogy(WL4./10, DataCube4(:,100,100), 'y--','LineWidth',1.5)
%legend('No Amplification', 'x100', 'Simulated 100 km wave');
legend('HIAMCM (scaled by SZA)','Simulated 100 km wave');

title('Spectra')
ylabel('Radiance R/nm')
xlabel('nm')
xlim([100,300])
%%
figure()
semilogy(WL3./10, DataCube3(:,100,100), 'k','LineWidth',1.5)
hold all
semilogy(WL3./10, abs(DataCube3(:,100,100)-DataCube4(:,100,100)), 'g--','LineWidth',1.5)
semilogy(WL4./10, DataCube4(:,100,100), 'y--','LineWidth',1.5)
%%
figure(3)
diff = DataCube2-DataCube;
plot(WL./10, diff(:,100,100))
title('Difference (x100 - w/o Amp)')
ylabel('Difference (R/nm)')
xlabel('nm')
%%

figure(4)
diff = DataCube2-DataCube;
plot(WL./10, diff(:,100,100))
title('Difference (x100 - w/o Amp)')
ylabel('Difference (R/nm)')
xlabel('nm')
xlim([100,300])
yticks(-20:5)
%findpeaks(diff(:,100,100))
%[pks, locs] = findpeaks(diff(:,100,100));
%text(x,y,txt)


%% lat lon
filename = "C:\Users\Maggie\Documents\DYNAGLO\FUVI-add-sensor-output\src\n2den.sav";
IDLVar=restore_idl_modach1(convertStringsToChars(filename));
n2den = IDLVar.TESTDATA;
  
%%
scenenum = "10";
filename = "\\lasp-store\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\inputs\angles\geolocated_scene_params_"+scenenum+"_256.sav";

%filename = "\\lasp-store.lasp.colorado.edu\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\outputs\DAY\dayglow_scene_"+ scenenum +"_256_v1_r0.detrend_on2.sav";
IDLVar=restore_idl_modach1(convertStringsToChars(filename));
% Set Parameters
lat = IDLVar.LAT;
lon = IDLVar.LON;

%%
figure
s = surf(lat, lon, n2den);
s.EdgeColor = 'none';
title("N2 den (No Amplification)")
ylabel('Longitude')
xlabel('Latitude')
view(0,90)
axis square

figure
s = surf(lat, lon, n2denamp);
s.EdgeColor = 'none';
title("N2 den (with Amplification x100)")
ylabel('Longitude')
xlabel('Latitude')
view(0,90)
axis square
%savename = 'C:\Users\Maggie\Documents\DYNAGLO\Hoskins Matlab\Files\MATLAB Runs\No Amplify\'+scenenum+ '_ch1m2'+'_'+string(binnum) + '_cb';
%h=colorbar();
%ylabel(h, 'photoevents')
%saveas(gcf, savename+'.png')

