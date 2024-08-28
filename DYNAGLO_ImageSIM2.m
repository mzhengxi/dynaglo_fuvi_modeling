function PhotoEvents = DYNAGLO_ImageSIM2(DataCube,WL,Ch,Int,Binning)
%% DYNAGLO FUVI Image Simulator
% Author Alan Hoskins, Alan.Hoskins@LASP.COLORADO.edu
%
% DYNAGLO FUVO Image Simulator will simulate an FUVI image when passed a
% data cube containing the atmospheric flux in R/A over a 40 degree by 40
% degree FOV.
%
% Usage: [PhotoEvents,WLList] = DYNAGLO_ImageSIM(DataCube,WL,Ch,Int,Binning)
%
%       Where:
%         - DataCube is a LxMxN data cube of flux in R/A. The first
%           dimension corresponds with the wavelength of the light and the
%           second two dimensions correpond with the FOV.  Each resolution
%           element ("resel") will be of angular size 40/M x 40/N degrees
%           and the flux will be assumed to be uniformally distributed in
%           the resel.  The flux distribution vs. wavelength will be binned
%           up in gaussian spectral lines below 200 nm and uniform 20 nm
%           wide distributions above 200 nm.  These distributions are
%           treated as probability function for the purpose of the image
%           simulation.
%         - WL is a Lx1 vector of wavelengths in Angstroms corresonding to
%           the first dimension of the DataCube.
%         - Ch is a vector with S elements corresponding to the images to
%           create.  For a vector [1,1,2,3,4,4,4], the result will be 2
%           channel 1 images, 1 channel 2 images, 1 channel 3 images and
%           3 channel 4 images.
%         - Int is a vector with the same size as Ch corresponding to the
%           integration time for each image.
%         - Binning is the binning of the resultant image.  If binning = 0
%           the resulting images will be photoevents from the FUVI
%           MCP.  If Binning > 0, it will covert these images to DN
%           from the binned CMOS detector.  The output size will be
%           [2052x2052] / Binning x S (the size of Ch & Int).  While the CMOS is
%           2048 x 2048, I have expanded the array size to make it easier to
%           bin 38 x 38.  For a Binning value of 38, the output will be a
%           54x54xS array.
%
if ~nargin % This section gets the first of Scott's newest data cubes and produces 5 binned detector images from each channel with 1.5 second integration
    tic
    IDLVar=restore_idl_modach('\\lasp-store.lasp.colorado.edu\projects\Phase_Development\DYNAGLO\6000 FUVI\6300 Simulations\HIAMCM\Polar\1200\00200\256\outputs\DAY\dayglow_scene_10_256_v1_r0.detrend_ampx100_on2.sav');
    DataCube=double(IDLVar.SPECT)/cosd(85);   % Divide this data set by cosine(85 deg) per Aimee
    WL=double(IDLVar.LAM);
    clearvars IDLVar
    %Ch=repelem([1,2,3],1,5);
    Ch=[1,2,3];
    Int=1.5*ones(size(Ch));
    Binning=38;
    PlotResults=false;
    fprintf('Imported data cube in %.1f seconds. \n',toc)
else
    PlotResults=false;
end

SensitivityFactor=0.5; % This derates the instrument sensitivity by the indicated amount    

DataCube=double(DataCube);
WL=double(WL);
%% Bin up the atomic lines
% First let's bin wavelengths between 1356 and 2000 A using Gaussian functions

% The atomic lines are very narrow so we can treat them as monochromatic.
tic
MeanFlux=sum(DataCube,[2,3])/prod(size(DataCube,[2,3])); % Rayeighs per A over the full FOV

AtomicLines=[1356,1493,1641,1743]; % Note, there is another N I line at 1412, but I have not seen it in Scott's data cubes
ABins=zeros(size(DataCube,2),size(DataCube,3),length(AtomicLines));
AWLList=zeros(length(AtomicLines),2);  % Initialize a list of wavelengths for DYNAGLO Model

for I=1:length(AtomicLines)
    DI=find(WL==AtomicLines(I));
    LineRatio=1-mean(MeanFlux(DI+[-1,1]))/MeanFlux(DI);
    ABins(:,:,I)=LineRatio*squeeze(DataCube(DI,:,:));
    AWLList(I,:)=[AtomicLines(I),0.3]*1e-4; % Converted to um, first column is wavelength, second is std dev assuming normal dist.
    DataCube(DI,:,:)=(1-LineRatio)*DataCube(DI,:,:);
end

DI=squeeze(any(ABins,[1,2]));

ABins=ABins(:,:,DI);
AWLList=AWLList(DI,:);

if PlotResults
    figure
    Dock
    hold on
    DI=WL>1350&WL<1999;
    plot(WL(DI),MeanFlux(DI)-sum(DataCube(DI,:,:),[2,3])/prod(size(DataCube,[2,3])))
    plot(WL(DI),sum(DataCube(DI,:,:),[2,3])/prod(size(DataCube,[2,3])))
    xlabel('Wavelength (A)')
    ylabel('Flux (R/A)')
    legend('Atomic Spectral Features','LBH Spectral Features')
    axis('tight')
    box on
    title({'Dyanglo Data Cube Average Flux vs. Wavelength, 135 - 200 nm','Atomic  Lines and LBH Spectra'})
    set(gca,'YScale','Linear')
end
%% Now bin up the LBH spectra.
% To do this, I am going to low pass the spectra by convolving with a
% 2 A wide Gaussian and then isolate the LBH lines and fit them with gaussian functions.

DI=WL>1350&WL<1999;  % Data index for wavelengths of interest
MeanFlux=zeros(size(MeanFlux));
MeanFlux(DI)=LowPass(sum(DataCube(DI,:,:),[2,3])/prod(size(DataCube,[2,3])),2,2); % Rayeighs per A over the full FOV

% I am just using a arbitrary isolated LBH feature to limit the dimmer features.
% The LBH spectrum, being a function of electron energy, has constant ratios that
% only change in shape and in magnitude vs. the atomic lines so this method
% should always work.

LLim=0.5*max(MeanFlux(WL>1646&WL<1651));

% Now locate all the spectral features
[PkV,FI]=findpeaks(MeanFlux(DI),WL(DI),'MinPeakHeight',LLim);

if ~isempty(FI)
    PDiff=[diff(FI);10];
    
    LBHWLList=zeros(length(FI),2);  % Initialize a list of wavelengths for DYNAGLO Model
    LBHBins=zeros(size(DataCube,2),size(DataCube,3),length(FI));  % Intialize an array of photon fluxes
    GFits=zeros(length(FI),3);
    TempFlux=zeros(size(MeanFlux));
    
    Window=7;
    
    I=1;
    while I<=length(FI)
        if PDiff(I)<=5
            ExtPks=find(PDiff(I:end)>10,1)-1;
            if isempty(ExtPks)
                ExtPks=length(FI)-I;
            end
        else
            ExtPks=0;
        end
        StartPts=[];
        Lower=[];
        Upper=[];
        FitType=['gauss',num2str(ExtPks+1)];
        for J=0:ExtPks
            StartPts=[StartPts,MeanFlux(WL==FI(I+J)),FI(I+J),2];
            Lower=[Lower,0.4*MeanFlux(WL==FI(I+J)),FI(I+J)-1,1];
            Upper=[Upper,MeanFlux(WL==FI(I+J)),FI(I+J)+1,6];
        end
        GFit=fit(((FI(I)-Window):(FI(I+ExtPks)+Window))',MeanFlux(find(WL==(FI(I)-Window)):find(WL==FI(I+ExtPks)+Window)),FitType,'Start',StartPts,'Lower',Lower,'Upper',Upper); % Fit peak to a gaussian
        for J=0:ExtPks
            GFits(I+J,:)=[GFit.(['a',num2str(J+1)]),GFit.(['b',num2str(J+1)]),GFit.(['c',num2str(J+1)])];
            LBHWLList(I+J,:)=[GFit.(['b',num2str(J+1)]),GFit.(['c',num2str(J+1)])/sqrt(2)]*1e-4; % Converted to um, first column is wavelength, second is std dev assuming normal dist.
            TempFlux=TempFlux+GFits(I+J,1)*exp(-((WL-GFits(I+J,2))/GFits(I+J,3)).^2);
            if ExtPks
                LineWeights=GFits(I+J,1)*exp(-((((FI(I)-Window):(FI(I+ExtPks)+Window))'-GFits(I+J,2))/GFits(I+J,3)).^2)./MeanFlux(find(WL==(FI(I)-Window)):find(WL==FI(I+ExtPks)+Window));
            else
                LineWeights=ones(size(((FI(I)-Window):(FI(I+ExtPks)+Window))'));
            end
            LBHBins(:,:,I+J)=squeeze(sum(repmat(LineWeights,1,size(DataCube,2),size(DataCube,3)).*DataCube(find(WL==(FI(I)-Window)):find(WL==FI(I+ExtPks)+Window),:,:),1));
        end
        I=I+ExtPks+1;
    end
    
    if PlotResults
        figure
        Dock
        plot(WL(DI),MeanFlux(DI))
        hold on
        plot(WL(DI),TempFlux(DI))
        plot(FI,PkV,'or')
        xlabel('Wavelength (A)')
        ylabel(['Mean Flux (R/',char(197),')'])
        legend('Data Cube Mean Flux','Gaussian Fit','LBH Bin Wavelengths')
        title('Fitted LBH spectra')
        
        figure
        Dock
        plot(FI,squeeze(sum(LBHBins,[1,2]))/prod(size(LBHBins,[1,2])))
        hold on
        plot(FI,GFits(:,1).*GFits(:,3)*sqrt(pi))
        xlabel('Wavelength (A)')
        ylabel('Integrated Flux (R)')
        legend('Integrated Data Cube','Gaussian Fits')
        title(sprintf('Flux Gain / Loss = %0.2f%%',100*(sum(TempFlux)-sum(MeanFlux))/sum(MeanFlux)))
    end
else
    LBHWLList=[];
    LBHBins=[];
end
%% Now bin up the rest of the spectrum using uniform distributions
BinSize=200; % 20 nm bins
RTWL=(2000+BinSize/2):BinSize:(max(WL)-BinSize/2);  % This is the bin centers
RTWLList=zeros(length(RTWL),2);
RTBins=zeros(size(DataCube,2),size(DataCube,3),length(RTWL));

for I=1:length(RTWL)
    DI=abs(WL-RTWL(I))<=BinSize/2; % Data index within each bin
    TempCube=DataCube(DI,:,:); % Grab the relevant bins to integrate.
    % This is the integrated power in Rayleighs at each spectral feature (trapazoidal rule)
    RTBins(:,:,I)=squeeze(sum((TempCube(1:(end-1),:,:)+...
        TempCube(2:end,:,:))/2.*repmat(diff(WL(DI)),1,size(TempCube,2),size(TempCube,3)),1));
    RTWLList(I,:)=[RTWL(I),-BinSize]*1e-4;  % The negative sign indicates a uniform distribution
end

DI=squeeze(any(RTBins,[1,2]));

RTBins=RTBins(:,:,DI);
RTWLList=RTWLList(DI,:);
%% Now put it all together into a wavelength list and an array of fluxes to raytrace and dump the rest

% BinnedFlux=LBHBins;
% WLList=LBHWLList;

BinnedFlux=cat(3,ABins,LBHBins,RTBins);
WLList=cat(1,AWLList,LBHWLList,RTWLList);
[~,SI]=sort(WLList(:,1));
WLList=WLList(SI,:);
BinnedFlux=BinnedFlux(:,:,SI);

fprintf('Binned wavelengths in %.1f seconds. \n',toc)

% Clean up the variable space
clearvars -except WLList BinnedFlux Int Ch Binning SensitivityFactor
%% Ray trace each wavelength bins through optical model
load('C:\Users\Maggie\Documents\DYNAGLO\Hoskins Matlab\CsI_QE.mat','CsI_QE')

% In simulating the detector images, I am creating 2052 x 2052 images just because it is
% easier to bin up being a multiple of 38.  The variable OutputImages will contain the
% number of photoevents at the output of the intensifier.

OutputImages=zeros(2048,2048,numel(Int));

if any(Int>0) % If the Integration time is greater than 0, raytrace, else just use a transmission model (like Greg's)
    if true   % Extrapolate the the binned data so we fill to the detector edge, this is more realisitic
        tic
        ExtraDist=7;  % Extrapolated pixels of data to add to edges of Scott's data cube
        
        ExtraFlux=zeros(size(BinnedFlux,1)+2*ExtraDist,size(BinnedFlux,2)+2*ExtraDist,size(BinnedFlux,3));
        FOV1=linspace(-20,20,size(BinnedFlux,1))';
        FOV2=linspace(-20,20,size(BinnedFlux,2))';
        [FXi,FYi]=meshgrid(FOV2,FOV1); % This is the x and y components of the FOV.
        
        FOV1=mean(diff(FOV1));  % This is the angular size of each pixel in the data cube
        FOV2=mean(diff(FOV2));  % This is the angular size of each pixel in the data cube
        
        FOV1=linspace(-(20+2*ExtraDist*FOV1),(20+2*ExtraDist*FOV1),size(ExtraFlux,1))';
        FOV2=linspace(-(20+2*ExtraDist*FOV2),(20+2*ExtraDist*FOV2),size(ExtraFlux,2))';
        [FXf,FYf]=meshgrid(FOV2,FOV1); % This is the x and y components of the FOV.
        
        for I=1:size(BinnedFlux,3)
            Slice=BinnedFlux(:,:,I);
            % F=scatteredInterpolant(FYi(:),FXi(:),Slice(:));
            F=griddedInterpolant(FYi,FXi,Slice,'linear','nearest');
            ExtraFlux(:,:,I)=F(FYf,FXf);
        end
        
        FOV1=mean(diff(FOV1));  % This is the angular size of each pixel in the data cube
        FOV2=mean(diff(FOV2));  % This is the angular size of each pixel in the data cube
        fprintf('Extrapolated binned data in %.1f seconds. \n',toc)
    else
        FOV1=linspace(-20,20,size(BinnedFlux,1))';
        FOV2=linspace(-20,20,size(BinnedFlux,2))';
        [FXf,FYf]=meshgrid(FOV2,FOV1); % This is the x and y components of the FOV.
        
        FOV1=mean(diff(FOV1));  % This is the angular size of each pixel in the data cube
        FOV2=mean(diff(FOV2));  % This is the angular size of each pixel in the data cube
        
        ExtraFlux=BinnedFlux;
    end
    for ImCount=1:numel(Int)  % Ray Trace each Extroplated Image
        tic
        PhotoEvents=zeros(2048,2048,size(BinnedFlux,3)); % This is number of photoevents per wavelength
        for I=1:size(WLList,1) % Step through each wavelength slice
            % Transmission through DYNAGLO.  This could also be left to the
            % raytrace code rather than assuming normal incidence and uniform thickness.
            Trans=1;                                        % Transmission = 1
            [n,A]=BaF2(WLList(I,1));                        % Index of refraction and absorption of BaF2
            Trans=Trans*(1-Fresnel(n(1)))^2*exp(-A*6);         % Transmission reduced by absorptioon and Fresnel reflection
            [n,A]=CaF2(WLList(I,1));                        % Index of refraction and absorption of CaF2
            Trans=Trans*(1-Fresnel(n(1)))^2*exp(-A*11);        % Transmission reduced by absorptioon and Fresnel reflection
            
            switch Ch(ImCount)
                case 2
                    [n,A]=Quartz(WLList(I,1));             % Index of refraction and absorption of Quartz
                    Trans=Trans*(1-Fresnel(n(:,1))).^2.*exp(-A*3);
                case 3
                    [n,A]=Spectrosil_2000(WLList(I,1));    % Index of refraction and absorption of Spectrosil 2000
                    Trans=Trans*(1-Fresnel(n(:,1))).^2.*exp(-A*3);
                case 4
                    [n,A]=SUPRASIL_3001(WLList(I,1));    % Index of refraction and absorption of Suprasil 3001
                    Trans=Trans*(1-Fresnel(n(:,1))).^2.*exp(-A*8);
            end
            
            % Etendue (throughput) is the Solid Angle (Omega) * Aperture Area (mm) * Cosine of the aperture tilt with FOV
            Etendue=FOV1*FOV2*(pi/180)^2*(pi*10^2*1e-6)*cosd(sqrt(FXf.^2+FYf.^2));
            
            % Photons per second at the aperture for current wavelength slice and each datacube pixel

            PhotonRate=poissrnd(ExtraFlux(:,:,I).*Etendue*1e10/(4*pi));
            
            % To get the correct distribution of photo-events at the intesifier,
            % the number of photons traced should equal the photon rate times the
            % transmission, the CsI QE, and the integration time in seconds.

            % IncidentPhotons=PhotonRate*SensitivityFactor*Trans*CsI_QE(WLList(I,1)*1000)*Int(ImCount);

            % The above line does not have excess noise from the MCP, the code below does:
            
            if mean(PhotonRate(:))>100 % Estimate probability with a normal distribution
                IncidentPhotons=poissrnd(PhotonRate*SensitivityFactor*Trans*CsI_QE(WLList(I,1)*1000)*Int(ImCount));
            else % If there are very few photons use the more correct binomial distribution
                IncidentPhotons=binornd(PhotonRate,SensitivityFactor*Trans*CsI_QE(WLList(I,1)*1000)*Int(ImCount));

                % The generation of a random number of electrons a single photon with gain
                % is a binomial distribution, but as the number gets large it is
                % indistinguishable from a poisson distribution.  Cacluating the random
                % binomial distribution is very slow.
            end
            
            % Find the number of pixels in current wavelength bin that acutal has incident photons

            DI=IncidentPhotons(:)>0;
           
            %  I need to break the following line into chuncks.  Say
            %  sum(IncidentPhotons(DI)) < 50 million, not sure where limit is but 150
            %  million incident photons causes the ray tracer to write to disk as it
            %  runs out of memory.
            
            % Original Code
            % PhotoEvents(:,:,I)=PhotoEvents(:,:,I)+DynPhotonTrace([FXf(DI),FYf(DI),FOV1*ones(size(FXf(DI))),FOV2*ones(size(FXf(DI)))],...
            %     WLList(I,:)/1000,IncidentPhotons(DI));
           
            % Divided into 50 million photon batches:
            BatchNum=unique(floor(cumsum(IncidentPhotons(DI))/5e7))';
            BatchDI=false(length(DI),length(BatchNum));
            BatchIndex1=1;
            for J=BatchNum
                BatchIndex2=find(DI,find(floor(cumsum(IncidentPhotons(DI))/5e7)==J,1,'last'));
                BatchDI(BatchIndex1:BatchIndex2(end),J+1)=DI(BatchIndex1:BatchIndex2(end));
                BatchIndex1=BatchIndex2(end)+1;
            end
            for J=1:size(BatchDI,2)
                PhotoEvents(:,:,I)=PhotoEvents(:,:,I)+DynPhotonTrace(...
                    [FXf(BatchDI(:,J)),FYf(BatchDI(:,J)),FOV1*ones(size(FXf(BatchDI(:,J)))),FOV2*ones(size(FXf(BatchDI(:,J))))],...
                    WLList(I,:)/1000,IncidentPhotons(BatchDI(:,J)));
            end
            %----------------------------------------------------------------------%
            clearvars BatchNum BatchDI BatchIndex1 BatchIndex2
        end
        % The line below sums the Output photoevents into a single image.  If you want to
        % look at specific wavelengths, you should comment this out though it will break
        % the conversion to CMOS image code below.

        OutputImages(:,:,ImCount)=sum(PhotoEvents,3);  
        fprintf('Ray traced full image in %.1f seconds. \n',toc)
    end
else
    % Generates images just using statistics (This should be similar to Greg's method in IDL)
    tic
    ExtraFlux=zeros(2048,2048,size(BinnedFlux,3));
    
    FOV1=linspace(-20,20,size(BinnedFlux,1))';
    FOV2=linspace(-20,20,size(BinnedFlux,2))';
    [FXi,FYi]=meshgrid(FOV2,FOV1); % This is the x and y components of the FOV for the Binned Cube.
    
    FOV1=linspace(-20,20,size(ExtraFlux,1))';
    FOV2=linspace(-20,20,size(ExtraFlux,2))';
    
    [FXf,FYf]=meshgrid(FOV2,FOV1); % This is the x and y components of the FOV for the high resolution Cube
    
    FOV1=mean(diff(FOV1));  % This is the angular size of each pixel in the data cube
    FOV2=mean(diff(FOV2));  % This is the angular size of each pixel in the data cube
    
    for I=1:size(BinnedFlux,3)
        Slice=BinnedFlux(:,:,I);
        F=scatteredInterpolant(FYi(:),FXi(:),Slice(:));
        ExtraFlux(:,:,I)=F(FYf,FXf);
    end
    
    Etendue=FOV1*FOV2*(pi/180)^2*(pi*10^2*1e-6)*cosd(sqrt(FXf.^2+FYf.^2));
    fprintf('Interpolated binned data in %.1f seconds. \n',toc)
    for ImCount=1:numel(Int)
        tic
        PhotoEvents=zeros(2048,2048,size(BinnedFlux,3)); % This is number of photoevents per wavelength
        for I=1:size(ExtraFlux,3)
            
            Trans=1;                                        % Transmission = 1
            [n,A]=BaF2(WLList(I,1));                        % Index of refraction and absorption of BaF2
            Trans=Trans*(1-Fresnel(n(1)))^2*exp(-A*6);         % Transmission reduced by absorptioon and Fresnel reflection
            [n,A]=CaF2(WLList(I,1));                        % Index of refraction and absorption of CaF2
            Trans=Trans*(1-Fresnel(n(1)))^2*exp(-A*11);        % Transmission
            
            switch Ch(ImCount)
                case 2
                    [n,A]=Quartz(WLList(I,1));             % Index of refraction and absorption of Quartz
                    Trans=Trans*(1-Fresnel(n(:,1))).^2.*exp(-A*3);
                case 3
                    [n,A]=Spectrosil_2000(WLList(I,1));    % Index of refraction and absorption of Spectrosil 2000
                    Trans=Trans*(1-Fresnel(n(:,1))).^2.*exp(-A*3);
                case 4
                    [n,A]=SUPRASIL_3001(WLList(I,1));    % Index of refraction and absorption of Suprasil 3001
                    Trans=Trans*(1-Fresnel(n(:,1))).^2.*exp(-A*8);
            end
            
            PhotonRate=poissrnd(ExtraFlux(:,:,I).*Etendue*1e10/(4*pi));  % Photon rate at aperature
            
            % PhotoEvents(:,:,I)=binornd(PhotonRate,Trans*CsI_QE(WLList(I,1)*1000)*SensitivityFactor*abs(Int(ImCount))); % Photoevents at image intensifier
            PhotoEvents(:,:,I)=poissrnd(PhotonRate*Trans*CsI_QE(WLList(I,1)*1000)*SensitivityFactor*abs(Int(ImCount))); % Photoevents at image intensifier
            % PhotoEvents(:,:,I)=PhotonRate*Trans*CsI_QE(WLList(I,1)*1000)*SensitivityFactor*abs(Int(ImCount)); % Photoevents at image intensifier - No Noise
        end
        OutputImages(:,:,ImCount)=sum(PhotoEvents,3);
    end
    fprintf('Simulated full image in %.1f seconds. \n',toc)
end
%% If passed binning in the function call, turn the PhotoEvents into a binned CMOS image
if Binning>0 % Create a Binned CMOS Image
    % Detector settings.
    PhosphorGain = 70;  %   (Photon / e-)
    FiberTransmission = 0.7^3;  % (Photon / Photon)
    MCPGain=42.1; % round(200/(PhosphorGain*FiberTransmission));  %   (e- / e-)
    CMOS_QE = 0.9;  % (e- / Photon)
    ReadNoise=16; % e-
    e_per_DN=10.5263;
    % Note, I should really use a binomial distribution, but as the number increases a
    % poisson distribution is nearly identical so it does not matter.  This makes the code
    % much much faster.
    tic
    OutputImages=(OutputImages.*poissrnd(MCPGain*PhosphorGain*FiberTransmission*CMOS_QE,size(OutputImages))+poissrnd(ReadNoise,size(OutputImages)))/e_per_DN;
    if Binning>1
        DetImages=zeros(size(OutputImages,1)/Binning,size(OutputImages,2)/Binning,size(OutputImages,3));
        for I=1:size(OutputImages,3)
            DetImages(:,:,I)=LinearBinImage(OutputImages(:,:,I),[],Binning);
        end
    end
    PhotoEvents=DetImages;
    fprintf('Converted photoevents to binned CMOS images in %.1f seconds. \n',toc)
else
    PhotoEvents=OutputImages;
end
end
%% Suppporting Functions
function ImageData=LinearBinImage(ImageData,Dim,Num)

if isempty(Dim)
    for I=1:ndims(ImageData)
        ImageData=LinearBinImage(ImageData,I,Num);
    end
else
    if Dim==1
        TempData=zeros(size(ImageData,1)/Num,size(ImageData,2));
        for I=1:size(TempData,1)
            TempData(I,:)=sum(ImageData((I-1)*Num+(1:Num),:),1);
        end
    else
        TempData=zeros(size(ImageData,1),size(ImageData,2)/Num);
        for I=1:size(TempData,2)
            TempData(:,I)=sum(ImageData(:,(I-1)*Num+(1:Num)),2);
        end
    end
    ImageData=TempData;
end
end
function DynImage=DynPhotonTrace(FOVXY,Lambda,PhotonsPerFP,DetBinning)

% Author: Alan Hoskins (Alan.Hoskins@LASP.Colorado.edu)

% Input variables
% FOVXY: mx2 or mx4 array of x,y field points in degrees.  If FOVXY has 2
%   columns, it will be treated as a point source at the specified field point.
%   If FOVXY has 4 column elements, the last two are the angular width of a
%   rectangular source in x and y.
% Lambda: scalar or 1x2 array designating wavelength in mm.  If scalar, the
%   source will be monochromatic.  If 1x2, the sources will have either a gaussian
%   power spectra centered on column 1 and a standard deviation given by if
%   column 2 is positive.  If column 2 is negative, it indicates the full
%   width of a uniform spectrum.
% P: scalar or mx1 array of powers (AU) correponding to the m field points.
%   If scalar, it is assumed every field point has the same power P.
% RaysPerFP: number of rays to trace per field point.  Note, if the memory
%   usage gets large enough (~0.5 of total RAM) that it starts paging data
%   to disk, it will slow down significantly. Better to loop the program.

if numel(PhotonsPerFP)==1
    PhotonsPerFP=PhotonsPerFP*ones(size(FOVXY,1),1);
end
TotalPhotons=sum(PhotonsPerFP);

%% Detector Setup
if nargin<4
    DetBinning=1;
end
DetSize=2048/DetBinning;  % I am making the detector size 2052*2052 so it is easier to bin up to 54x54.
DetPitch=0.0137*DetBinning;
DetX=(1:(DetSize))*DetPitch;
DetX=DetX-mean(DetX);
Pix=polyfit(DetX,1:DetSize,1);
DynImage=zeros(DetSize);
%% Optical System Setup.
% Below are the origins of the surfaces in mm.  For planar surfaces, it is simply
% a point on the plane, for spherical, it is the center of the sphere.

% For Dynaglow, the first surface is the front of the BaF2 window at z=0.
% The second surface is the other side of the window at z=6 mm.

O = [0,0,0;
    0,0,6;
    0,0,127.6;
    0,0,-8.6;
    0,0,55.35;
    0,0,64.35]';

% Index of refraction of the next surface. Can be numeric of a function
% name.

n={'BaF2',1,'CaF2',1,'MgF2',1};

% Radius of Curvature for the surfaces - 0 indicates a planar surface

R=[0,0,114,-33.2,0,0]; % Radius of curvature

% Aperture Diameters
% D=2*[10,12.5,20,20,31.25,20.5]'; % Aperture diameter
D=2*[10,12.5,20,20,31.25,28]'; % Aperture diameter

% Aperture offset in x,y
DO=zeros(6,2);
%% Set up Ray Matrices

% Replicate the wavelenth vector - same length as number of rays.
if numel(Lambda)==1
    Lambda=Lambda*ones(1,TotalPhotons);                % Monochromatic
else
    if Lambda(2)>0
        Lambda=Lambda(1)+Lambda(2)*randn(1,TotalPhotons);  % Random Gaussian Distribution
    else
        Lambda=Lambda(1)-Lambda(2)*(rand(1,TotalPhotons)-0.5);  % Random uniform distribution
    end
end

% Here are the initial k vectors.  For a recangular source, the
% distribution of rays is random in the rectangle
if size(FOVXY,2)==2
    ki=sind(repelem(FOVXY,PhotonsPerFP,1))'; % Each field point it a point source
else
    % Each field point is a rectangular source
    ki=sind(repelem(FOVXY(:,1:2),PhotonsPerFP,1)+repelem(FOVXY(:,3:4),PhotonsPerFP,1).*(2*rand(TotalPhotons,2)-1)/2)';
end
ki(3,:)=sqrt(1-RSS(ki,1).^2);
ki=2*pi*ki./Lambda;  % The k vectors have a length of 2*pi*n/lambda in inverse mm

% The initial index of refraction for the rays
ni=ones(size(Lambda));

% Randomaly generate the photon positions in the aperture;

if TotalPhotons<1000
    x=D(1)*(2*rand(1,ceil(4*TotalPhotons))-1)/2;
    y=D(1)*(2*rand(1,ceil(4*TotalPhotons))-1)/2;
else
    x=D(1)*(2*rand(1,ceil(1.4*TotalPhotons))-1)/2;
    y=D(1)*(2*rand(1,ceil(1.4*TotalPhotons))-1)/2;
end

ri=[x(RSS([x;y],1)<=D(1)/2);y(RSS([x;y],1)<=D(1)/2)];
ri=[ri(:,1:TotalPhotons);zeros(1,TotalPhotons)];

clearvars x y

% Ray polarization, only matters if we use the fresnel reflection and
% transmission coefficients.
% Pol=logical(round(rand(1,size(ki,2)))); % Random polarization, S is true, P is false(1,size(ki,2))
%% Trace the rays.  Step through each surface
for I=1:length(R)
    % First find the parametric variable t to trace the ray to the next surface
    if R(I) % This finds t for a spherical surface
        ro=ri-repmat(O(:,I),1,size(ri,2));
        t=-(2*MatDot(ki,ro)+...
            sign(RSS(ro,1)-abs(R(I))).*sqrt((2*MatDot(ki,ro)).^2-...
            4*MatDot(ki,ki).*(MatDot(ro,ro)-R(I)^2)))./(2*MatDot(ki,ki));
    else
        % This finds t for a planar surface normal to the optical axis
        % (adding tilt is easy, but I stripped it out for Dynaglo.
        N=[zeros(2,size(ki,2));sign(ki(3,:))];
        t=MatDot(N,(repmat(O(:,I),1,size(ri,2))-ri))./MatDot(N,ki);
    end
    
    if ~isreal(t) % This removes rays that miss the next component - should never happen
        DI=abs(imag(t))==0;
        t=t(DI);
        ri=ri(:,DI);
        ki=ki(:,DI);
        ro=ro(:,DI);
        ni=ni(DI);
        Lambda=Lambda(DI);
    end
    
    r=ri+ki.*repmat(t,3,1); % The (x,y,z) intersection with the next surface for each ray
    
    if R(I)  % Find the normal vector for a spherical surface at each ray intersection.
        N=sign(RSS(ro,1)-abs(R(I))).*(repmat(O(:,I),1,size(r,2))-r);
        N=N./sqrt(repmat(MatDot(N,N),3,1));
    end
    
    % Get the index of refraction of the next material and ratio it with the current index
    if ischar(n{I})
        n2=eval([n{I},'(Lambda''*1e3)']);
        n12=ni./n2(:,1)';
    else
        n12=ni./repmat(n{I},1,length(Lambda));
    end
    
    % Find the transmitted k vector using Snell's law
    k=SnellsLaw(ki,N,n12);
    
    % Adjust the power using Fresnel equations. Note Pol=true is S polarization
    % RC=Fresnel(ki,k,n12,Pol,N);
    % T=1-RC; % Transmission Coefficient
    % P=T.*P; % We are only transmitting in Dynaglo
    
    % Calculate the reduced power due to absorption. Not used for Dynaglo
    % P=P.*exp(Alpha*RSS(ki.*repmat(t,3,1),1));
    
    % Find the rays that are vignetted and remove them
    DI=RSS(r(1:2,:)-repmat(DO(I,:)',1,size(r,2)),1)<=D(I)/2;
    
    r=r(:,DI);
    k=k(:,DI);
    ni=ni(DI);
    n12=n12(DI);
    Lambda=Lambda(DI);
    
    % Overwrite the rays positions, k vectors, and indices with the new values
    ri=r;
    ki=k;
    ni=ni./n12;
end
%% Image creation
% Find the pixel location the rays are striking
r=round(polyval(Pix,r(1:2,:)));

% Remove rays which miss the detector
DI=all(r>=1&r<=DetSize,1);
r=r(:,DI);

%Build up the image with the propagated ray powers
for I=1:size(r,2)
    DynImage(r(2,I),r(1,I))=DynImage(r(2,I),r(1,I))+1;
end
end
function k=SnellsLaw(k,N,n12)

% Author Alan Hoskins (Alan.Hoskins@Lasp.colorado.edu)
%
% Snells law for a ray entering a region of index2 from a region of index1
%
% Function Usage:  kout = SnellsLaw(kin,N,n12)
% Where kin is a 3xn array of column vectors representing the propagation
% vectors of different plane waves in region 1.
% N is a 3x1 or 3xn array of surface unit normal vectors pointing from the
% region of index 2 into the region of index 1.
% n12 is the ratio of the index of refraction in region 1 to the index in
% region 2 (n1/n2).
% kout is a 3xn array of column vectors representing the rays in region 2.
% The magnitude of kout vectors will by kin/n12 (i.e., assumes the vectors
% are k vectors and will obey consevation of momentum).
% Run without inputs for a drawing.

if ~nargin
    k=[sind(30),0,cosd(30)]';
    N=[0,0,1]';
    n12=1/1.5;
    kmag=1;
    figure
    Pos=get(gcf,'Position');
    set(gcf,'Position',[Pos(1:2)-[0,680-Pos(4)],680,680])
    xlim([-2 2])
    ylim([-2 2])
    daspect([1,1,1])
    box on
    hold on
    DrawRect([1,0],[2,4],0,'k','b',0.05)
    line([0;0],[-2,2],'color','k')
    DrawArrow([0;N(3)],[0,N(1)],'k',0.1*kmag)
    DrawArrow([-k(3);0],[-k(1);0],'r')
    DrawCircle([0,0],2*RSS(k,1),0.7*ones(1,3),'w',0)
    kR=SnellsLaw(k,N,0);
    DrawArrow([0;kR(3)],[0,kR(1)],'r',0.1*kmag/RSS(k,1))
    ko=SnellsLaw(k,N,n12);
    RC0=Fresnel(k,ko,n12,false,N);
    RC1=Fresnel(k,ko,n12,true,N);
    DrawArrow([0;ko(3)],[0,ko(1)],'r',0.1*kmag/RSS(ko,1))
    DrawCircle([0,0],2*RSS(ko,1),0.7*ones(1,3),'w',0)
    text(-1,1.75,'n=1')
    text(1,1.75,'n=1.5','HorizontalAlignment','right')
    text(-0.5,-0.4,'k_{in}')
    text(-0.5,0.4,'k_{ref}')
    text(0.3,0.4,'k_{out}')
    text(0.5,-0.15,'N')
    text(-1.9,-1.5,{'k = [0.5; 0; 0.866]','N=[0; 0; 1]','n12 = 1/1.5','SnellsLaw(k,N,n12) = [0.5; 0; 1.4142]'})
    text(-1.9,1.35,{'Reflection Coefficient',sprintf('S_{Pol} = %4.1f%%',RC1*100),sprintf('P_{Pol} = %4.1f%%',RC0*100)})
    xlabel('z')
    ylabel('x')
else
    if size(N,2)<size(k,2)
        N=repmat(N,1,size(k,2));
    end
    
    DI=n12>0; % Check for Reflection, occurs when n12=0, completely non-physical
    
    kT=k-repmat(MatDot(k,N),3,1).*N; % transverse ray vector
    if any(DI)
        k(:,DI)=kT(:,DI)+repmat(sqrt((RSS(k(:,DI),1)./n12(DI)).^2-RSS(kT(:,DI),1).^2),3,1).*N(:,DI); % transmitted rays
    end
    if any(~DI)
        k(:,~DI)=k(:,~DI)-2*repmat(MatDot(k(:,~DI),N(:,~DI)),3,1).*N(:,~DI); % reflected rays
    end
end
end
function [n,varargout]=BaF2(wl)

if nargin==0
    wl=0.58756;
end

A=[6.43356000E-001;
    5.06762000E-001;
    3.82610000E+000];
B=[3.34000000E-003;
    1.20300000E-002;
    2.15169810E+003];

n=Sellmeier(A,B,wl(:));

if nargout>1
    % Note: this absorbtion coefficient (mm^-1) is for -15 C LASP data pieced
    % together with Data thief from Korth website for 2 and 32 mm thick
    % Plus Thorlabs data for 10 mm thick sample
    
    Alpha=[0.118	2.4
        0.119	2.389560282
        0.12	2.35894925
        0.121	2.309229081
        0.122	2.241461958
        0.123	2.156710059
        0.124	2.056035566
        0.125	1.940500656
        0.126	1.811167512
        0.127	1.669098312
        0.128	1.515355238
        0.129	1.351000467
        0.13	1.177096182
        0.131	0.994704561
        0.132	0.804887785
        0.133	0.608661754
        0.134	0.289172876
        0.135	0.128372951
        0.136	0.070262342
        0.137	0.056356784
        0.138	0.052915741
        0.139	0.052233277
        0.14	0.051527824
        0.141	0.050802682
        0.142	0.050649251
        0.143	0.04989734
        0.144	0.048535327
        0.145	0.045951576
        0.146	0.044886235
        0.147	0.041315811
        0.148	0.037058586
        0.149	0.032040887
        0.15	0.026993769
        0.151	0.022636782
        0.152	0.019997655
        0.153	0.017145701
        0.154	0.016263552
        0.155	0.01606159
        0.156	0.015869389
        0.157	0.015686507
        0.158	0.015512507
        0.159	0.015346957
        0.16	0.015189425
        0.161	0.015039484
        0.162	0.014896707
        0.163	0.014760673
        0.164	0.014630961
        0.165	0.014507151
        0.166	0.014388827
        0.167	0.014275574
        0.168	0.014166978
        0.169	0.014062628
        0.17	0.013962112
        0.171	0.013865022
        0.172	0.013770949
        0.173	0.013679486
        0.174	0.013590227
        0.175	0.013502767
        0.176	0.013416702
        0.177	0.013331628
        0.178	0.013247142
        0.179	0.013162844
        0.18	0.013078332
        0.181	0.012993206
        0.182	0.012907068
        0.183	0.012819518
        0.184	0.012730159
        0.185	0.012638595
        0.186	0.01254443
        0.187	0.012447269
        0.188	1.23E-02
        0.189	1.22E-02
        0.19	1.21E-02
        0.191	1.20E-02
        0.192	1.19E-02
        0.193	1.18E-02
        0.194	1.17E-02
        0.195	1.15E-02
        0.196	1.14E-02
        0.197	1.12E-02
        0.198	1.11E-02
        0.199	1.09E-02
        0.2	1.07E-02
        0.201	1.06E-02
        0.202	1.04E-02
        0.203	1.02E-02
        0.204	9.95E-03
        0.205	9.53E-03
        0.206	9.12E-03
        0.207	8.57E-03
        0.208	7.57E-03
        0.209	6.40E-03
        0.21	5.48E-03
        0.211	4.56E-03
        0.212	4.06E-03
        0.213	3.58E-03
        0.214	3.28E-03
        0.215	3.01E-03
        0.216	2.82E-03
        0.217	2.60E-03
        0.218	2.46E-03
        0.219	2.18E-03
        0.22	2.14E-03
        0.221	1.98E-03
        0.222	2.00E-03
        0.223	1.86E-03
        0.224	1.75E-03
        0.225	1.58E-03
        0.226	1.51E-03
        0.227	1.45E-03
        0.228	1.41E-03
        0.229	1.36E-03
        0.23	1.32E-03
        0.231	1.22E-03
        0.232	1.18E-03
        0.233	1.14E-03
        0.234	1.13E-03
        0.235	1.10E-03
        0.236	9.85E-04
        0.237	9.22E-04
        0.238	8.57E-04
        0.239	9.57E-04
        0.24	8.37E-04
        0.241	9.00E-04
        0.242	9.44E-04
        0.243	8.41E-04
        0.244	9.21E-04
        0.245	8.72E-04
        0.246	7.86E-04
        0.247	8.10E-04
        0.248	8.52E-04
        0.249	7.65E-04
        0.25	8.06E-04
        0.251	7.73E-04
        0.252	7.95E-04
        0.253	8.17E-04
        0.254	8.01E-04
        0.255	6.76E-04
        0.256	6.97E-04
        0.257	7.54E-04
        0.258	7.01E-04
        0.259	7.20E-04
        0.26	7.76E-04
        0.261	6.86E-04
        0.262	7.23E-04
        0.263	7.24E-04
        0.264	6.69E-04
        0.265	6.87E-04
        0.266	7.23E-04
        0.267	5.95E-04
        0.268	6.49E-04
        0.269	6.66E-04
        0.27	5.92E-04
        0.271	6.09E-04
        0.272	6.61E-04
        0.273	6.23E-04
        0.274	6.75E-04
        0.275	6.37E-04
        0.276	6.52E-04
        0.277	6.86E-04
        0.278	7.55E-04
        0.279	7.70E-04
        0.28	8.76E-04
        0.281	9.64E-04
        0.282	1.03E-03
        0.283	1.14E-03
        0.284	1.28E-03
        0.285	1.33E-03
        0.286	1.44E-03
        0.287	1.47E-03
        0.288	1.50E-03
        0.289	1.57E-03
        0.29	1.65E-03
        0.291	1.50E-03
        0.292	1.51E-03
        0.293	1.42E-03
        0.294	1.27E-03
        0.295	1.17E-03
        0.296	9.80E-04
        0.297	8.10E-04
        0.298	7.31E-04
        0.299	6.70E-04
        0.3	5.72E-04
        0.301	6.20E-04
        0.302	5.23E-04
        0.303	4.07E-04
        0.304	4.18E-04
        0.305	3.92E-04
        0.306	4.58E-04
        0.307	4.68E-04
        0.308	3.88E-04
        0.309	3.08E-04
        0.31	3.91E-04
        0.311	3.83E-04
        0.312	3.75E-04
        0.313	3.48E-04
        0.314	3.58E-04
        0.315	3.68E-04
        0.316	3.96E-04
        0.317	3.69E-04
        0.318	3.42E-04
        0.319	4.06E-04
        0.32	4.15E-04
        0.321	3.88E-04
        0.322	4.51E-04
        0.323	4.60E-04
        0.324	3.97E-04
        0.325	4.23E-04
        0.326	4.32E-04
        0.327	4.22E-04
        0.328	4.13E-04
        0.329	4.03E-04
        0.33	3.93E-04
        0.331	4.20E-04
        0.332	4.10E-04
        0.333	3.63E-04
        0.334	4.07E-04
        0.335	3.79E-04
        0.336	3.33E-04
        0.337	3.95E-04
        0.338	3.66E-04
        0.339	3.37E-04
        0.34	3.45E-04
        0.341	3.70E-04
        0.342	3.78E-04
        0.343	4.03E-04
        0.344	3.74E-04
        0.345	3.27E-04
        0.346	3.70E-04
        0.347	3.59E-04
        0.348	3.84E-04
        0.349	4.27E-04
        0.35	3.80E-04
        0.351	3.32E-04
        0.352	3.57E-04
        0.353	3.45E-04
        0.354	3.34E-04
        0.355	3.76E-04
        0.356	3.47E-04
        0.357	3.17E-04
        0.358	3.05E-04
        0.359	3.29E-04
        0.36	3.36E-04
        0.361	3.78E-04
        0.362	3.48E-04
        0.363	2.63E-04
        0.364	2.69E-04
        0.365	3.11E-04
        0.366	3.35E-04
        0.367	3.77E-04
        0.368	3.29E-04
        0.369	2.44E-04
        0.37	2.50E-04
        0.371	2.73E-04
        0.372	2.61E-04
        0.373	2.84E-04
        0.374	2.71E-04
        0.375	2.76E-04
        0.376	2.82E-04
        0.377	2.88E-04
        0.378	2.97E-04
        0.379	3.11E-04
        0.38	3.22E-04
        0.381	3.31E-04
        0.382	3.41E-04
        0.383	3.47E-04
        0.384	3.48E-04
        0.385	3.47E-04
        0.386	3.46E-04
        0.387	3.45E-04
        0.388	3.44E-04
        0.389	3.41E-04
        0.39	3.39E-04
        0.391	3.35E-04
        0.392	3.27E-04
        0.393	3.22E-04
        0.394	3.19E-04
        0.395	3.17E-04
        0.396	3.17E-04
        0.397	3.21E-04
        0.398	3.23E-04
        0.399	3.24E-04
        0.4	3.23E-04
        0.401	3.20E-04
        0.402	3.18E-04
        0.403	3.16E-04
        0.404	3.14E-04
        0.405	3.15E-04
        0.406	3.17E-04
        0.407	3.17E-04
        0.408	3.17E-04
        0.409	3.20E-04
        0.41	3.20E-04
        0.411	3.19E-04
        0.412	3.15E-04
        0.413	3.10E-04
        0.414	3.04E-04
        0.415	3.01E-04
        0.416	2.98E-04
        0.417	3.00E-04
        0.418	3.03E-04
        0.419	3.05E-04
        0.42	3.06E-04
        0.421	3.09E-04
        0.422	3.11E-04
        0.423	3.11E-04
        0.424	3.11E-04
        0.425	3.08E-04
        0.426	3.03E-04
        0.427	3.00E-04
        0.428	3.02E-04
        0.429	3.03E-04
        0.43	3.06E-04
        0.431	3.08E-04
        0.432	3.05E-04
        0.433	3.00E-04
        0.434	2.98E-04
        0.435	2.99E-04
        0.436	3.01E-04
        0.437	3.05E-04
        0.438	3.09E-04
        0.439	3.12E-04
        0.44	3.12E-04
        0.441	3.15E-04
        0.442	3.17E-04
        0.443	3.20E-04
        0.444	3.21E-04
        0.445	3.21E-04
        0.446	3.19E-04
        0.447	3.20E-04
        0.448	3.18E-04
        0.449	3.16E-04
        0.45	3.14E-04
        0.451	3.14E-04
        0.452	3.11E-04
        0.453	3.09E-04
        0.454	3.06E-04
        0.455	2.99E-04
        0.456	2.92E-04
        0.457	2.93E-04
        0.458	2.95E-04
        0.459	2.96E-04
        0.46	3.00E-04
        0.461	3.04E-04
        0.462	2.99E-04
        0.463	2.95E-04
        0.464	2.94E-04
        0.465	2.95E-04
        0.466	2.96E-04
        0.467	2.99E-04
        0.468	3.02E-04
        0.469	3.08E-04
        0.47	3.11E-04
        0.471	3.15E-04
        0.472	3.21E-04
        0.473	3.20E-04
        0.474	3.13E-04
        0.475	3.07E-04
        0.476	2.98E-04
        0.477	2.92E-04
        0.478	2.86E-04
        0.479	2.86E-04
        0.48	2.87E-04
        0.481	2.89E-04
        0.482	2.92E-04
        0.483	2.99E-04
        0.484	3.03E-04
        0.485	3.04E-04
        0.486	3.06E-04
        0.487	3.08E-04
        0.488	3.08E-04
        0.489	3.11E-04
        0.49	3.14E-04
        0.491	3.17E-04
        0.492	3.20E-04
        0.493	3.24E-04
        0.494	3.25E-04
        0.495	3.28E-04
        0.496	3.26E-04
        0.497	3.17E-04
        0.498	3.10E-04
        0.499	3.08E-04
        0.5	3.02E-04
        0.501	2.99E-04
        0.502	3.02E-04
        0.503	3.02E-04
        0.504	2.94E-04
        0.505	2.90E-04
        0.506	2.88E-04
        0.507	2.80E-04
        0.508	2.74E-04
        0.509	2.72E-04
        0.51	2.76E-04
        0.511	2.83E-04
        0.512	2.95E-04
        0.513	3.03E-04
        0.514	3.08E-04
        0.515	3.01E-04
        0.516	2.90E-04
        0.517	2.80E-04
        0.518	2.74E-04
        0.519	2.68E-04
        0.52	2.66E-04
        0.521	2.64E-04
        0.522	2.58E-04
        0.523	2.56E-04
        0.524	2.57E-04
        0.525	2.56E-04
        0.526	2.55E-04
        0.527	2.51E-04
        0.528	2.44E-04
        0.529	2.39E-04
        0.53	2.38E-04
        0.531	2.36E-04
        0.532	2.36E-04
        0.533	2.33E-04
        0.534	2.27E-04
        0.535	2.23E-04
        0.536	2.25E-04
        0.537	2.25E-04
        0.538	2.29E-04
        0.539	2.29E-04
        0.54	2.25E-04
        0.541	2.26E-04
        0.542	2.27E-04
        0.543	2.25E-04
        0.544	2.26E-04
        0.545	2.24E-04
        0.546	2.17E-04
        0.547	2.14E-04
        0.548	2.15E-04
        0.549	2.17E-04
        0.55	2.18E-04
        0.551	2.20E-04
        0.552	2.23E-04
        0.553	2.29E-04
        0.554	2.31E-04
        0.555	2.34E-04
        0.556	2.36E-04
        0.557	2.30E-04
        0.558	2.15E-04
        0.559	2.06E-04
        0.56	2.04E-04
        0.561	2.01E-04
        0.562	2.01E-04
        0.563	2.08E-04
        0.564	2.07E-04
        0.565	1.96E-04
        0.566	1.93E-04
        0.567	1.94E-04
        0.568	1.96E-04
        0.569	2.02E-04
        0.57	2.15E-04
        0.571	2.17E-04
        0.572	2.16E-04
        0.573	2.14E-04
        0.574	2.13E-04
        0.575	2.08E-04
        0.576	2.12E-04
        0.577	2.14E-04
        0.578	2.14E-04
        0.579	2.17E-04
        0.58	2.22E-04
        0.581	2.24E-04
        0.582	2.24E-04
        0.583	2.25E-04
        0.584	2.23E-04
        0.585	2.23E-04
        0.586	2.22E-04
        0.587	2.24E-04
        0.588	2.21E-04
        0.589	2.18E-04
        0.59	2.13E-04
        0.591	2.10E-04
        0.592	2.04E-04
        0.593	2.02E-04
        0.594	1.98E-04
        0.595	1.91E-04
        0.596	1.86E-04
        0.597	1.87E-04
        0.598	1.87E-04
        0.599	1.86E-04
        0.6	1.88E-04
        0.601	1.89E-04
        0.602	1.85E-04
        0.603	1.85E-04
        0.604	1.88E-04
        0.605	1.91E-04
        0.606	1.90E-04
        0.607	1.92E-04
        0.608	1.95E-04
        0.609	1.98E-04
        0.61	1.98E-04
        0.611	2.01E-04
        0.612	1.96E-04
        0.613	1.90E-04
        0.614	1.85E-04
        0.615	1.84E-04
        0.616	1.83E-04
        0.617	1.87E-04
        0.618	1.85E-04
        0.619	1.79E-04
        0.62	1.73E-04
        0.621	1.69E-04
        0.622	1.66E-04
        0.623	1.69E-04
        0.624	1.69E-04
        0.625	1.69E-04
        0.626	1.67E-04
        0.627	1.64E-04
        0.628	1.58E-04
        0.629	1.58E-04
        0.63	1.58E-04
        0.631	1.59E-04
        0.632	1.59E-04
        0.633	1.64E-04
        0.634	1.63E-04
        0.635	1.63E-04
        0.636	1.58E-04
        0.637	1.55E-04
        0.638	1.54E-04
        0.639	1.55E-04
        0.64	1.56E-04
        0.641	1.59E-04
        0.642	1.58E-04
        0.643	1.53E-04
        0.644	1.49E-04
        0.645	1.45E-04
        0.646	1.44E-04
        0.647	1.44E-04
        0.648	1.43E-04
        0.649	1.41E-04
        0.65	1.36E-04
        0.651	1.30E-04
        0.652	1.21E-04
        0.653	1.09E-04
        0.654	9.87E-05
        0.655	9.18E-05
        0.656	8.62E-05
        0.657	8.29E-05
        0.658	8.52E-05
        0.659	8.33E-05
        0.66	7.99E-05
        0.661	7.72E-05
        0.662	7.59E-05
        0.663	7.40E-05
        0.664	7.78E-05
        0.665	7.94E-05
        0.666	8.03E-05
        0.667	8.26E-05
        0.668	8.56E-05
        0.669	8.65E-05
        0.67	8.74E-05
        0.671	8.47E-05
        0.672	7.62E-05
        0.673	6.85E-05
        0.674	6.29E-05
        0.675	5.87E-05
        0.676	5.59E-05
        0.677	5.39E-05
        0.678	4.90E-05
        0.679	4.34E-05
        0.68	3.85E-05
        0.681	3.36E-05
        0.682	3.16E-05
        0.683	2.52E-05
        0.684	1.38E-05
        0.685	3.16E-06
        0.686	0.00E+00
        0.687	0.00E+00
        0.688	4.17E-06
        0.689	1.58E-05
        0.69	2.59E-05
        0.691	2.89E-05
        0.692	3.04E-05
        0.693	2.76E-05
        0.694	2.19E-05
        0.695	1.48E-05
        0.696	9.12E-06
        0.697	2.72E-06
        0.698	1.35E-06
        0.699	5.00E-06
        0.7	1.08E-05
        0.701	1.66E-05
        0.702	2.31E-05
        0.703	2.24E-05
        0.704	1.60E-05
        0.705	9.59E-06
        0.706	1.01E-06
        0.707	0.00E+00
        0.708	0.00E+00
        0.709	0.00E+00
        0.71	0.00E+00
        0.711	0.00E+00
        0.712	0.00E+00
        0.713	0.00E+00
        0.714	0.00E+00
        0.715	0.00E+00
        0.716	0.00E+00
        0.717	0.00E+00
        0.718	0.00E+00
        0.719	0.00E+00
        0.72	0.00E+00
        0.721	0.00E+00
        0.722	0.00E+00
        0.723	0.00E+00
        0.724	0.00E+00
        0.725	0.00E+00
        0.726	0.00E+00
        0.727	0.00E+00
        0.728	0.00E+00
        0.729	0.00E+00
        0.73	0.00E+00
        0.731	0.00E+00
        0.732	0.00E+00
        0.733	0.00E+00
        0.734	0.00E+00
        0.735	0.00E+00
        0.736	0.00E+00
        0.737	0.00E+00
        0.738	0.00E+00
        0.739	0.00E+00
        0.74	0.00E+00
        0.741	0.00E+00
        0.742	0.00E+00
        0.743	0.00E+00
        0.744	0.00E+00
        0.745	0.00E+00
        0.746	0.00E+00
        0.747	0.00E+00
        0.748	0.00E+00
        0.749	0.00E+00
        0.75	0.00E+00
        0.751	0.00E+00
        0.752	0.00E+00
        0.753	0.00E+00
        0.754	0.00E+00
        0.755	0.00E+00
        0.756	0.00E+00
        0.757	0.00E+00
        0.758	0.00E+00
        0.759	0.00E+00
        0.76	0.00E+00
        0.761	0.00E+00
        0.762	0.00E+00
        0.763	0.00E+00
        0.764	0.00E+00
        0.765	0.00E+00
        0.766	0.00E+00
        0.767	0.00E+00
        0.768	0.00E+00
        0.769	0.00E+00
        0.77	0.00E+00
        0.771	0.00E+00
        0.772	0.00E+00
        0.773	0.00E+00
        0.774	0.00E+00
        0.775	0.00E+00
        0.776	0.00E+00
        0.777	0.00E+00
        0.778	0.00E+00
        0.779	0.00E+00
        0.78	0.00E+00
        0.781	0.00E+00
        0.782	0.00E+00
        0.783	0.00E+00
        0.784	0.00E+00
        0.785	0.00E+00
        0.786	0.00E+00
        0.787	0.00E+00
        0.788	0.00E+00
        0.789	0.00E+00
        0.79	0.00E+00
        0.791	0.00E+00
        0.792	0.00E+00
        0.793	0.00E+00
        0.794	0.00E+00
        0.795	0.00E+00
        0.796	0.00E+00
        0.797	0.00E+00
        0.798	0.00E+00
        0.799	0.00E+00
        0.8	0.00E+00
        0.801	0.00E+00
        0.802	0.00E+00
        0.803	0.00E+00
        0.804	0.00E+00
        0.805	0.00E+00
        0.806	0.00E+00
        0.807	0.00E+00
        0.808	0.00E+00
        0.809	0.00E+00
        0.81	0.00E+00
        0.811	0.00E+00
        0.812	0.00E+00
        0.813	0.00E+00
        0.814	0.00E+00
        0.815	0.00E+00
        0.816	0.00E+00
        0.817	0.00E+00
        0.818	0.00E+00
        0.819	0.00E+00
        0.82	0.00E+00
        0.821	0.00E+00
        0.822	0.00E+00
        0.823	0.00E+00
        0.824	0.00E+00
        0.825	0.00E+00
        0.826	0.00E+00
        0.827	0.00E+00
        0.828	0.00E+00
        0.829	0.00E+00
        0.83	0.00E+00
        0.831	0.00E+00
        0.832	0.00E+00
        0.833	0.00E+00
        0.834	0.00E+00
        0.835	0.00E+00
        0.836	0.00E+00
        0.837	0.00E+00
        0.838	0.00E+00
        0.839	0.00E+00
        0.84	0.00E+00
        0.841	0.00E+00
        0.842	0.00E+00
        0.843	0.00E+00
        0.844	0.00E+00
        0.845	0.00E+00
        0.846	0.00E+00
        0.847	0.00E+00
        0.848	0.00E+00
        0.849	0.00E+00
        0.85	0.00E+00
        0.851	0.00E+00
        0.852	0.00E+00
        0.853	0.00E+00
        0.854	0.00E+00
        0.855	0.00E+00
        0.856	0.00E+00
        0.857	0.00E+00
        0.858	0.00E+00
        0.859	0.00E+00
        0.86	0.00E+00
        0.861	0.00E+00
        0.862	0.00E+00
        0.863	0.00E+00
        0.864	0.00E+00
        0.865	0.00E+00
        0.866	0.00E+00
        0.867	0.00E+00
        0.868	0.00E+00
        0.869	0.00E+00
        0.87	0.00E+00
        0.871	0.00E+00
        0.872	0.00E+00
        0.873	0.00E+00
        0.874	0.00E+00
        0.875	0.00E+00
        0.876	0.00E+00
        0.877	0.00E+00
        0.878	0.00E+00
        0.879	0.00E+00
        0.88	0.00E+00
        0.881	0.00E+00
        0.882	0.00E+00
        0.883	0.00E+00
        0.884	0.00E+00
        0.885	0.00E+00
        0.886	0.00E+00
        0.887	0.00E+00
        0.888	0.00E+00
        0.889	0.00E+00
        0.89	0.00E+00
        0.891	0.00E+00
        0.892	0.00E+00
        0.893	0.00E+00
        0.894	0.00E+00
        0.895	0.00E+00
        0.896	0.00E+00
        0.897	0.00E+00
        0.898	0.00E+00
        0.899	0.00E+00
        0.9	0.00E+00
        0.901	0.00E+00
        0.902	0.00E+00
        0.903	0.00E+00
        0.904	0.00E+00
        0.905	0.00E+00
        0.906	0.00E+00
        0.907	0.00E+00
        0.908	0.00E+00
        0.909	0.00E+00
        0.91	0.00E+00
        0.911	0.00E+00
        0.912	0.00E+00
        0.913	0.00E+00
        0.914	0.00E+00
        0.915	0.00E+00
        0.916	0.00E+00
        0.917	0.00E+00
        0.918	0.00E+00
        0.919	0.00E+00
        0.92	0.00E+00
        0.921	0.00E+00
        0.922	0.00E+00
        0.923	0.00E+00
        0.924	0.00E+00
        0.925	0.00E+00
        0.926	0.00E+00
        0.927	0.00E+00
        0.928	0.00E+00
        0.929	0.00E+00
        0.93	0.00E+00
        0.931	0.00E+00
        0.932	0.00E+00
        0.933	0.00E+00
        0.934	0.00E+00
        0.935	0
        0.936	0
        0.937	0
        0.938	0
        0.939	0
        0.94	0
        0.941	0
        0.942	0
        0.943	0
        0.944	0
        0.945	0
        0.946	0
        0.947	0
        0.948	0.00E+00
        0.949	0.00E+00
        0.95	0.00E+00
        0.951	0.00E+00
        0.952	0.00E+00
        0.953	0.00E+00
        0.954	0.00E+00
        0.955	0.00E+00
        0.956	0.00E+00
        0.957	0.00E+00
        0.958	0.00E+00
        0.959	0.00E+00
        0.96	0.00E+00
        0.961	0.00E+00
        0.962	0.00E+00
        0.963	0.00E+00
        0.964	0.00E+00
        0.965	0.00E+00
        0.966	0.00E+00
        0.967	0.00E+00
        0.968	0.00E+00
        0.969	0.00E+00
        0.97	0.00E+00
        0.971	0.00E+00
        0.972	0.00E+00
        0.973	0.00E+00
        0.974	0.00E+00
        0.975	0.00E+00
        0.976	0.00E+00
        0.977	0.00E+00
        0.978	0.00E+00
        0.979	0.00E+00
        0.98	0.00E+00
        0.981	0.00E+00
        0.982	0.00E+00
        0.983	0.00E+00
        0.984	0.00E+00
        0.985	0.00E+00
        0.986	0.00E+00
        0.987	0.00E+00
        0.988	0.00E+00
        0.989	0.00E+00
        0.99	0.00E+00
        0.991	0.00E+00
        0.992	0.00E+00
        0.993	0.00E+00
        0.994	0
        0.995	0
        0.996	0
        0.997	0
        0.998	0
        0.999	0
        1	0];
    
    AlphaFit=fit(Alpha(:,1),Alpha(:,2),'linearinterp');
    varargout{1}=AlphaFit(wl);
end
end
function [n,alpha]=MgF2(wl)

if nargin==0
    wl=0.58756;
end

% These values are from M. J. Dodge, "Refractive properties of magnesium
% fluoride," APPL. OPT., v23, n12, pg 1980 (1984).  Valid for 200-700 nm.

Ao=[4.87551080e-1;3.98750310e-1;2.31203530e0];
Bo=[1.88217800e-3;8.95188847e-3;5.66135591e2];

Ae=[4.13440230e-1;5.04974990e-1;2.49048620e0];
Be=[1.35737865e-3;8.23767167e-3;5.65107755e2];

no=Sellmeier(Ao,Bo,wl);
ne=Sellmeier(Ae,Be,wl);

n=[no,ne];
alpha=0;


% These values are from M. W. Williams and E. T. Arakawa, "Optical
% properties of crystalline MgF 2 from 115 nm to 400 nm,"
% Appl. Opt. 18,1477-1478 (1979). - Note 22 degrees C
%
% lno=[115 1.714;
%     120 1.626;
%     125 1.588;
%     130 1.555;
%     140 1.513;
%     150 1.484;
%     160 1.464;
%     170 1.451;
%     180 1.440;
%     190 1.430;
%     200 1.422;
%     225 1.408;
%     250 1.401;
%     275 1.397;
%     300 1.393;
%     350 1.387;
%     400 1.382;
%     450 1.381;
%     500 1.381];
% lno(:,1)=lno(:,1)*1e-3;
%
% lnomne=[1759.0 13.91;
%     1724.3 13.99;
%     1690.6 14.06;
%     1657.8 14.12;
%     1626.0 14.18;
%     1594.0 14.22;
%     1561.9 14.25;
%     1529.5 14.27;
%     1494.6 14.24;
%     1453.2 14.14;
%     1350.5 13.14;
%     1322.8 12.61;
%     1315.5 12.27;
%     1305.9 11.91;
%     1297.1 11.57;
%     1289.3 11.24;
%     1282.6 10.92;
%     1277.2 10.62;
%     1272.2 10.32;
%     1267.4 10.02;
%     1263.4 9.73;
%     1259.5 9.45;
%     1255.8 9.17;
%     1252.4 8.89;
%     1249.3 8.61;
%     1246.3 8.34;
%     1243.6 8.07;
%     1240.9 7.80;
%     1238.4 7.53;
%     1235.0 7.26;
%     1233.8 7.00;
%     1231.5 6.74;
%     1229.5 6.48;
%     1227.7 6.22;
%     1225.6 5.96;
%     1223.8 5.70;
%     1222.1 5.45;
%     1220.2 5.20;
%     1218.7 4.94;
%     1217.2 4.69;
%     1216.0 4.44;
%     1214.4 4.18;
%     1213.2 3.93;
%     1212.0 3.68;
%     1210.0 3.43;
%     1208.7 3.18;
%     1207.3 2.93;
%     1205.9 2.69;
%     1205.0 2.44;
%     1203.7 2.19;
%     1202.8 1.95];
% lnomne(:,1)=lnomne(:,1)*1e4;
% lnomne(:,2)=lnomne(:,1)*1e-3;

%
% if nargin==0
%     wl=0.58756;
% end
%
% A=[2.85932450;
%     -1.02100270e-2;
%     3.19731620e-2;
%     2.02068690e-3;
%     -1.13993140e-4;
%     1.40129320e-5];
%
% n=Schott(A,wl);
end
function [n,varargout]=CaF2(wl)

if nargin==0
    wl=0.58756;
end

A=[5.67588800E-001;
    4.71091400E-001;
    3.84847230E+000];
B=[2.52643000E-003;
    1.00783330E-002;
    1.20055600E+003];

n=Sellmeier(A,B,wl(:));


if nargout>1
    % Absorption coefficient in mm^-1
    Alpha=[0.110000000000000	1.28430646846247
        0.111000000000000	1.07825006720856
        0.112000000000000	0.958205419245258
        0.113000000000000	0.873265954121463
        0.114000000000000	0.807504518466631
        0.115000000000000	0.753840747927856
        0.116000000000000	0.708505726511836
        0.117000000000000	0.669254502840269
        0.118000000000000	0.634641986223773
        0.119000000000000	0.574824068649236
        0.120000000000000	0.359311150478784
        0.121000000000000	0.296462794316046
        0.122000000000000	0.156421898119388
        0.123000000000000	0.0861145637439351
        0.124000000000000	0.0455028535011549
        0.125000000000000	0.0398952765512230
        0.126000000000000	0.0363524772500327
        0.127000000000000	0.0338977028305250
        0.128000000000000	0.0326874604361131
        0.129000000000000	0.0319782311289669
        0.130000000000000	0.0317418661689235
        0.131000000000000	0.0314913758885426
        0.132000000000000	0.0312282522846040
        0.133000000000000	0.0309538084432533
        0.134000000000000	0.0306557024791690
        0.135000000000000	0.0303366064475427
        0.136000000000000	0.0300095441261654
        0.137000000000000	0.0296753272362559
        0.138000000000000	0.0293346815504506
        0.139000000000000	0.0289882577262072
        0.140000000000000	0.0286366405470780
        0.141000000000000	0.0278850914089917
        0.142000000000000	0.0268879531804315
        0.143000000000000	0.0258972794677854
        0.144000000000000	0.0249132041068699
        0.145000000000000	0.0238404765025241
        0.146000000000000	0.0223672343639224
        0.147000000000000	0.0209139883022991
        0.148000000000000	0.0194802688793388
        0.149000000000000	0.0180726041656534
        0.150000000000000	0.0166983967931005
        0.151000000000000	0.0153417439055045
        0.152000000000000	0.0140022488721340
        0.153000000000000	0.0126825507692439
        0.154000000000000	0.0113943115326540
        0.155000000000000	0.0101216117142216
        0.156000000000000	0.00886411406924494
        0.157000000000000	0.00814953045020291
        0.158000000000000	0.00753491449142098
        0.159000000000000	0.00692255996601626
        0.160000000000000	0.00640477928486103
        0.161000000000000	0.00616542548071510
        0.162000000000000	0.00592492085834969
        0.163000000000000	0.00568335679824236
        0.164000000000000	0.00544081877082289
        0.165000000000000	0.00525753071568251
        0.166000000000000	0.00512295226330226
        0.167000000000000	0.00498709339113420
        0.168000000000000	0.00485002398175003
        0.169000000000000	0.00471180986003419
        0.170000000000000	0.00457251307961646
        0.171000000000000	0.00443219218537440
        0.172000000000000	0.00429090245430801
        0.173000000000000	0.00414869611683849
        0.174000000000000	0.00400562256036442
        0.175000000000000	0.00377942227594700
        0.176000000000000	0.00353995220766542
        0.177000000000000	0.00330018688495718
        0.178000000000000	0.00306016297288242
        0.179000000000000	0.00281991512611095
        0.180000000000000	0.00257947611475436
        0.181000000000000	0.00233887694095764
        0.182000000000000	0.00209814694703147
        0.183000000000000	0.00186846981284714
        0.184000000000000	0.00181814851610722
        0.185000000000000	0.00177575499114940
        0.186000000000000	0.00173264574725429
        0.187000000000000	0.00168884714989710
        0.188000000000000	0.00164438437858472
        0.189000000000000	0.00159928149262967
        0.190000000000000	0.00155356149258730
        0.191000000000000	0.00150724637768726
        0.192000000000000	0.00146035719956101
        0.193000000000000	0.00141291411254106
        0.194000000000000	0.00136493642078570
        0.195000000000000	0.00131644262246004
        0.196000000000000	0.00126745045118661
        0.197000000000000	0.00121797691496039
        0.198000000000000	0.00116803833270789
        0.199000000000000	0.00111765036865487
        0.200000000000000	0.00106682806465526
        0.201000000000000	0.00101558587062084
        0.202000000000000	0.000963937673180728
        0.203000000000000	0.000917234219813944
        0.204000000000000	0.000920168109200581
        0.205000000000000	0.000922731464610719
        0.206000000000000	0.000924888823799264
        0.207000000000000	0.000926651710094418
        0.208000000000000	0.000928031224235215
        0.209000000000000	0.000929038063518515
        0.210000000000000	0.000929682539913852
        0.211000000000000	0.000929974597209956
        0.212000000000000	0.000929923827253047
        0.213000000000000	0.000929539485332354
        0.214000000000000	0.000928830504765333
        0.215000000000000	0.000927805510730042
        0.216000000000000	0.000926472833390367
        0.217000000000000	0.000924840520355569
        0.218000000000000	0.000922916348513330
        0.219000000000000	0.000920707835273010
        0.220000000000000	0.000918222249252901
        0.221000000000000	0.000915466620443476
        0.222000000000000	0.000912447749876384
        0.223000000000000	0.000909172218827146
        0.224000000000000	0.000905646397577335
        0.225000000000000	0.000901876453761032
        0.226000000000000	0.000897868360318099
        0.227000000000000	0.000893627903075946
        0.228000000000000	0.000889160687979870
        0.229000000000000	0.000884472147990812
        0.230000000000000	0.000879567549668097
        0.231000000000000	0.000874451999454243
        0.232000000000000	0.000869130449677093
        0.233000000000000	0.000863607704284123
        0.234000000000000	0.000857888424322964
        0.235000000000000	0.000851977133180654
        0.236000000000000	0.000845878221594794
        0.237000000000000	0.000839595952447032
        0.238000000000000	0.000833134465350755
        0.239000000000000	0.000826497781042837
        0.240000000000000	0.000819757764605396
        0.241000000000000	0.000814101419488426
        0.242000000000000	0.000808311420497567
        0.243000000000000	0.000802360738222028
        0.244000000000000	0.000796252858921448
        0.245000000000000	0.000789991172721979
        0.246000000000000	0.000783578976879919
        0.247000000000000	0.000777019478913629
        0.248000000000000	0.000770315799609740
        0.249000000000000	0.000763470975909576
        0.250000000000000	0.000756487963681363
        0.251000000000000	0.000749369640383235
        0.252000000000000	0.000742118807622523
        0.253000000000000	0.000734738193615264
        0.254000000000000	0.000727230455551258
        0.255000000000000	0.000719598181868250
        0.256000000000000	0.000711843894439430
        0.257000000000000	0.000703970050678561
        0.258000000000000	0.000695979045565461
        0.259000000000000	0.000687873213596111
        0.260000000000000	0.000679654830660386
        0.261000000000000	0.000671326115850136
        0.262000000000000	0.000662889233201252
        0.263000000000000	0.000654346293372166
        0.264000000000000	0.000645821690250701
        0.265000000000000	0.000642158058365388
        0.266000000000000	0.000638725083522452
        0.267000000000000	0.000635192511587873
        0.268000000000000	0.000631562210082016
        0.269000000000000	0.000627836002026157
        0.270000000000000	0.000624015667244521
        0.271000000000000	0.000620102943621477
        0.272000000000000	0.000616099528314860
        0.273000000000000	0.000612007078927915
        0.274000000000000	0.000607827214641024
        0.275000000000000	0.000603561517304851
        0.276000000000000	0.000599211532496855
        0.277000000000000	0.000594778770541826
        0.278000000000000	0.000590264707498936
        0.279000000000000	0.000585670786115178
        0.280000000000000	0.000580998416748099
        0.281000000000000	0.000576248978257583
        0.282000000000000	0.000571423818868755
        0.283000000000000	0.000566524257006771
        0.284000000000000	0.000561551582104521
        0.285000000000000	0.000556507055384387
        0.286000000000000	0.000551391910615088
        0.287000000000000	0.000546207354844163
        0.288000000000000	0.000540954569107522
        0.289000000000000	0.000535634709116351
        0.290000000000000	0.000530248905922852
        0.291000000000000	0.000524798266564760
        0.292000000000000	0.000519283874690315
        0.293000000000000	0.000513706791163588
        0.294000000000000	0.000508068054651363
        0.295000000000000	0.000502368682192090
        0.296000000000000	0.000496609669747553
        0.297000000000000	0.000490791992737631
        0.298000000000000	0.000484916606559346
        0.299000000000000	0.000478984447089997
        0.300000000000000	0.000475182491286938
        0.301000000000000	0.000477979221230070
        0.302000000000000	0.000480720724385887
        0.303000000000000	0.000483407863551520
        0.304000000000000	0.000486041484305799
        0.305000000000000	0.000488622415429561
        0.306000000000000	0.000491151469314388
        0.307000000000000	0.000493629442358967
        0.308000000000000	0.000496057115354623
        0.309000000000000	0.000498435253859369
        0.310000000000000	0.000500764608561687
        0.311000000000000	0.000503045915633757
        0.312000000000000	0.000505279897074928
        0.313000000000000	0.000507467261045384
        0.314000000000000	0.000509608702190607
        0.315000000000000	0.000511704901956837
        0.316000000000000	0.000513756528897703
        0.317000000000000	0.000515764238972534
        0.318000000000000	0.000517728675836413
        0.319000000000000	0.000519650471122250
        0.320000000000000	0.000521530244715457
        0.321000000000000	0.000523368605020775
        0.322000000000000	0.000525166149222363
        0.323000000000000	0.000526923463536607
        0.324000000000000	0.000528641123458161
        0.325000000000000	0.000530319693999781
        0.326000000000000	0.000531959729925357
        0.327000000000000	0.000533561775977216
        0.328000000000000	0.000535126367097105
        0.329000000000000	0.000536654028641753
        0.330000000000000	0.000538145276592414
        0.331000000000000	0.000539600617759403
        0.332000000000000	0.000541020549981011
        0.333000000000000	0.000542405562317582
        0.334000000000000	0.000543756135240325
        0.335000000000000	0.000545072740815586
        0.336000000000000	0.000546355842884373
        0.337000000000000	0.000547605897237189
        0.338000000000000	0.000548823351784613
        0.339000000000000	0.000550008646723681
        0.340000000000000	0.000551162214699786
        0.341000000000000	0.000552284480965055
        0.342000000000000	0.000553375863532125
        0.343000000000000	0.000554436773324936
        0.344000000000000	0.000555467614325069
        0.345000000000000	0.000556468783714967
        0.346000000000000	0.000557440672017532
        0.347000000000000	0.000558383663232255
        0.348000000000000	0.000559298134968186
        0.349000000000000	0.000560184458573554
        0.350000000000000	0.000561042999262427
        0.351000000000000	0.000561874116238248
        0.352000000000000	0.000562678162814370
        0.353000000000000	0.000563455486532000
        0.354000000000000	0.000564206429274927
        0.355000000000000	0.000564931327381946
        0.356000000000000	0.000565630511756461
        0.357000000000000	0.000566304307973437
        0.358000000000000	0.000566953036384124
        0.359000000000000	0.000567577012218059
        0.360000000000000	0.000568176545683002
        0.361000000000000	0.000568751942062233
        0.362000000000000	0.000569303501809940
        0.363000000000000	0.000569831520644242
        0.364000000000000	0.000570336289638102
        0.365000000000000	0.000570818095308300
        0.366000000000000	0.000571277219702149
        0.367000000000000	0.000571713940482609
        0.368000000000000	0.000572128531011136
        0.369000000000000	0.000572521260428862
        0.370000000000000	0.000572892393736050
        0.371000000000000	0.000573242191869533
        0.372000000000000	0.000573570911778679
        0.373000000000000	0.000573878806499465
        0.374000000000000	0.000574166125227099
        0.375000000000000	0.000574433113387099
        0.376000000000000	0.000574680012704477
        0.377000000000000	0.000574907061271847
        0.378000000000000	0.000575114493615804
        0.379000000000000	0.000575302540761841
        0.380000000000000	0.000575471430298073
        0.381000000000000	0.000575621386437434
        0.382000000000000	0.000575752630078560
        0.383000000000000	0.000575865378865373
        0.384000000000000	0.000575959847245529
        0.385000000000000	0.000576036246527377
        0.386000000000000	0.000576094784936059
        0.387000000000000	0.000576135667667979
        0.388000000000000	0.000576159096944657
        0.389000000000000	0.000576165272065016
        0.390000000000000	0.000576154389456711
        0.391000000000000	0.000576126642726607
        0.392000000000000	0.000576082222709993
        0.393000000000000	0.000576021317518662
        0.394000000000000	0.000575944112588396
        0.395000000000000	0.000575850790725070
        0.396000000000000	0.000575741532150152
        0.397000000000000	0.000575616514544987
        0.398000000000000	0.000575475913094399
        0.399000000000000	0.000575319900529343
        0.400000000000000	0.000575148647168607
        0.401000000000000	0.000574962320959787
        0.402000000000000	0.000574761087519435
        0.403000000000000	0.000574545110172389
        0.404000000000000	0.000574314549990244
        0.405000000000000	0.000574069565829191
        0.406000000000000	0.000573810314367008
        0.407000000000000	0.000573536950139334
        0.408000000000000	0.000573249625575314
        0.409000000000000	0.000572948491032449
        0.410000000000000	0.000572633694830802
        0.411000000000000	0.000572305383286504
        0.412000000000000	0.000571963700744706
        0.413000000000000	0.000571608789611798
        0.414000000000000	0.000571240790387007
        0.415000000000000	0.000570859841693466
        0.416000000000000	0.000570466080308618
        0.417000000000000	0.000570059641194011
        0.418000000000000	0.000569640657524694
        0.419000000000000	0.000569209260717781
        0.420000000000000	0.000568765580460784
        0.421000000000000	0.000568309744739050
        0.422000000000000	0.000567841879863118
        0.423000000000000	0.000567362110495155
        0.424000000000000	0.000566870559675141
        0.425000000000000	0.000566367348846466
        0.426000000000000	0.000565852597881071
        0.427000000000000	0.000565326425104189
        0.428000000000000	0.000564788947318490
        0.429000000000000	0.000564240279827938
        0.430000000000000	0.000563680536461013
        0.431000000000000	0.000563109829593782
        0.432000000000000	0.000562528270172178
        0.433000000000000	0.000561935967734362
        0.434000000000000	0.000561333030432124
        0.435000000000000	0.000560719565052329
        0.436000000000000	0.000560095677037771
        0.437000000000000	0.000559461470507709
        0.438000000000000	0.000558817048278003
        0.439000000000000	0.000558162511880937
        0.440000000000000	0.000557497961584671
        0.441000000000000	0.000556823496412194
        0.442000000000000	0.000556139214160313
        0.443000000000000	0.000555445211417818
        0.444000000000000	0.000554741583583828
        0.445000000000000	0.000554028424885320
        0.446000000000000	0.000553305828394802
        0.447000000000000	0.000552573886047297
        0.448000000000000	0.000551832688657336
        0.449000000000000	0.000551082325935427
        0.450000000000000	0.000550322886504341
        0.451000000000000	0.000549554457915186
        0.452000000000000	0.000548777126662982
        0.453000000000000	0.000547990978202177
        0.454000000000000	0.000547196096961893
        0.455000000000000	0.000546392566360718
        0.456000000000000	0.000545580468821386
        0.457000000000000	0.000544759885785266
        0.458000000000000	0.000543930897726459
        0.459000000000000	0.000543093584165775
        0.460000000000000	0.000542248023684270
        0.461000000000000	0.000541394293936970
        0.462000000000000	0.000540532471665889
        0.463000000000000	0.000539662632713017
        0.464000000000000	0.000538784852033290
        0.465000000000000	0.000537899203707017
        0.466000000000000	0.000537005760952306
        0.467000000000000	0.000536104596137164
        0.468000000000000	0.000535195780791512
        0.469000000000000	0.000534279385618956
        0.470000000000000	0.000533355480508333
        0.471000000000000	0.000532424134545035
        0.472000000000000	0.000531485416022204
        0.473000000000000	0.000530539392451909
        0.474000000000000	0.000529586130575703
        0.475000000000000	0.000528625696375440
        0.476000000000000	0.000527658155083710
        0.477000000000000	0.000526683571194193
        0.478000000000000	0.000525702008471687
        0.479000000000000	0.000524713529962132
        0.480000000000000	0.000523718198002551
        0.481000000000000	0.000522716074230415
        0.482000000000000	0.000521707219593482
        0.483000000000000	0.000520691694358934
        0.484000000000000	0.000519669558122580
        0.485000000000000	0.000518640869818020
        0.486000000000000	0.000517605687725481
        0.487000000000000	0.000516564069480617
        0.488000000000000	0.000515516072083065
        0.489000000000000	0.000514461751905068
        0.490000000000000	0.000513401164699785
        0.491000000000000	0.000512334365609400
        0.492000000000000	0.000511261409173449
        0.493000000000000	0.000510182349336581
        0.494000000000000	0.000509097239456617
        0.495000000000000	0.000508006132312064
        0.496000000000000	0.000506909080109846
        0.497000000000000	0.000505806134492848
        0.498000000000000	0.000504697346547180
        0.499000000000000	0.000503582766809430
        0.500000000000000	0.000502462445273985
        0.501000000000000	0.000501336431399875
        0.502000000000000	0.000500204774117842
        0.503000000000000	0.000499067521837068
        0.504000000000000	0.000497924722451999
        0.505000000000000	0.000496776423348899
        0.506000000000000	0.000495622671412447
        0.507000000000000	0.000494463513032014
        0.508000000000000	0.000493298994108194
        0.509000000000000	0.000492129160058848
        0.510000000000000	0.000490954055825395
        0.511000000000000	0.000489773725878788
        0.512000000000000	0.000488588214225387
        0.513000000000000	0.000487397564413041
        0.514000000000000	0.000486201819536661
        0.515000000000000	0.000485001022244075
        0.516000000000000	0.000483795214741491
        0.517000000000000	0.000482584438799133
        0.518000000000000	0.000481368735756729
        0.519000000000000	0.000480148146528676
        0.520000000000000	0.000478922711609551
        0.521000000000000	0.000477692471079208
        0.522000000000000	0.000476457464608004
        0.523000000000000	0.000475217731461667
        0.524000000000000	0.000473973310506489
        0.525000000000000	0.000472724240214222
        0.526000000000000	0.000471470558666701
        0.527000000000000	0.000470212303560989
        0.528000000000000	0.000468949512213717
        0.529000000000000	0.000467682221565931
        0.530000000000000	0.000466410468187624
        0.531000000000000	0.000465134288282253
        0.532000000000000	0.000463853717691117
        0.533000000000000	0.000462568791897793
        0.534000000000000	0.000461279546032470
        0.535000000000000	0.000459986014876105
        0.536000000000000	0.000458688232864773
        0.537000000000000	0.000457386234093629
        0.538000000000000	0.000456080052321145
        0.539000000000000	0.000454769720973042
        0.540000000000000	0.000453455273146217
        0.541000000000000	0.000452136741612771
        0.542000000000000	0.000450814158823812
        0.543000000000000	0.000449487556913152
        0.544000000000000	0.000448156967701178
        0.545000000000000	0.000446822422698570
        0.546000000000000	0.000445483953109842
        0.547000000000000	0.000444141589836971
        0.548000000000000	0.000442795363482977
        0.549000000000000	0.000441445304355391
        0.550000000000000	0.000440091442469666
        0.551000000000000	0.000438733807552684
        0.552000000000000	0.000437372429045988
        0.553000000000000	0.000436007336109181
        0.554000000000000	0.000434638557623099
        0.555000000000000	0.000433266122193166
        0.556000000000000	0.000431890058152401
        0.557000000000000	0.000430510393564739
        0.558000000000000	0.000429127156227870
        0.559000000000000	0.000427740373676570
        0.560000000000000	0.000426350073185577
        0.561000000000000	0.000424956281772373
        0.562000000000000	0.000423559026200517
        0.563000000000000	0.000422158332982205
        0.564000000000000	0.000420754228381177
        0.565000000000000	0.000419346738415682
        0.566000000000000	0.000417935888861037
        0.567000000000000	0.000416521705252556
        0.568000000000000	0.000415104212888137
        0.569000000000000	0.000413683436830943
        0.570000000000000	0.000412259401912041
        0.571000000000000	0.000410832132733023
        0.572000000000000	0.000409401653668561
        0.573000000000000	0.000407967988868884
        0.574000000000000	0.000406531162262358
        0.575000000000000	0.000405091197557861
        0.576000000000000	0.000403648118247274
        0.577000000000000	0.000402201947607853
        0.578000000000000	0.000400752708704591
        0.579000000000000	0.000399300424392576
        0.580000000000000	0.000397845117319329
        0.581000000000000	0.000396386809927000
        0.582000000000000	0.000394925524454612
        0.583000000000000	0.000393461282940446
        0.584000000000000	0.000391994107224010
        0.585000000000000	0.000390524018948321
        0.586000000000000	0.000389051039562063
        0.587000000000000	0.000387575190321615
        0.588000000000000	0.000386096492293186
        0.589000000000000	0.000384614966354794
        0.590000000000000	0.000383130633198463
        0.591000000000000	0.000381643513332030
        0.592000000000000	0.000380153627081243
        0.593000000000000	0.000378660994591650
        0.594000000000000	0.000377165635830598
        0.595000000000000	0.000375667570589065
        0.596000000000000	0.000374166818483554
        0.597000000000000	0.000372663398958001
        0.598000000000000	0.000371157331285587
        0.599000000000000	0.000369648634570424
        0.600000000000000	0.000368137327749576
        0.601000000000000	0.000366623429594651
        0.602000000000000	0.000365106958713590
        0.603000000000000	0.000363587933552370
        0.604000000000000	0.000362066372396825
        0.605000000000000	0.000360542293374085
        0.606000000000000	0.000359015714454479
        0.607000000000000	0.000357486653453060
        0.608000000000000	0.000355955128031202
        0.609000000000000	0.000354421155698283
        0.610000000000000	0.000352884753813169
        0.611000000000000	0.000351345939585864
        0.612000000000000	0.000349804730078996
        0.613000000000000	0.000348261142209307
        0.614000000000000	0.000346715192749274
        0.615000000000000	0.000345166898328477
        0.616000000000000	0.000343616275435179
        0.617000000000000	0.000342063340417576
        0.618000000000000	0.000340508109485510
        0.619000000000000	0.000338950598711637
        0.620000000000000	0.000337390824033003
        0.621000000000000	0.000335828801252297
        0.622000000000000	0.000334264546039309
        0.623000000000000	0.000332698073932238
        0.624000000000000	0.000331129400338947
        0.625000000000000	0.000329558540538507
        0.626000000000000	0.000327985509682218
        0.627000000000000	0.000326410322795085
        0.628000000000000	0.000324832994777082
        0.629000000000000	0.000323253540404255
        0.630000000000000	0.000321671974330180
        0.631000000000000	0.000320088311086993
        0.632000000000000	0.000318502565086711
        0.633000000000000	0.000316914750622440
        0.634000000000000	0.000315324881869454
        0.635000000000000	0.000313732972886521
        0.636000000000000	0.000312139037616932
        0.637000000000000	0.000310543089889698
        0.638000000000000	0.000308945143420637
        0.639000000000000	0.000307345211813561
        0.640000000000000	0.000305743308561332
        0.641000000000000	0.000304139447046952
        0.642000000000000	0.000302533640544670
        0.643000000000000	0.000300925902221026
        0.644000000000000	0.000299316245135921
        0.645000000000000	0.000297704682243611
        0.646000000000000	0.000296091226393784
        0.647000000000000	0.000294475890332622
        0.648000000000000	0.000292858686703654
        0.649000000000000	0.000291239628048954
        0.650000000000000	0.000289618726809969
        0.651000000000000	0.000287995995328505
        0.652000000000000	0.000286371445847793
        0.653000000000000	0.000284745090513349
        0.654000000000000	0.000283116941373984
        0.655000000000000	0.000281487010382575
        0.656000000000000	0.000279855309397188
        0.657000000000000	0.000278221850181860
        0.658000000000000	0.000276586644407561
        0.659000000000000	0.000274949703653028
        0.660000000000000	0.000273311039405655
        0.661000000000000	0.000271670663062408
        0.662000000000000	0.000270028585930655
        0.663000000000000	0.000268384819228950
        0.664000000000000	0.000266739374088056
        0.665000000000000	0.000265092261551486
        0.666000000000000	0.000263443492576634
        0.667000000000000	0.000261793078035356
        0.668000000000000	0.000260141028714949
        0.669000000000000	0.000258487355318774
        0.670000000000000	0.000256832068467165
        0.671000000000000	0.000255175178698202
        0.672000000000000	0.000253516696468331
        0.673000000000000	0.000251856632153398
        0.674000000000000	0.000250194996049071
        0.675000000000000	0.000248531798371913
        0.676000000000000	0.000246867049259862
        0.677000000000000	0.000245200758773053
        0.678000000000000	0.000243532936894610
        0.679000000000000	0.000241863593531230
        0.680000000000000	0.000240192738513912
        0.681000000000000	0.000238520381598785
        0.682000000000000	0.000236846532467636
        0.683000000000000	0.000235171200728616
        0.684000000000000	0.000233494395917044
        0.685000000000000	0.000231816127495928
        0.686000000000000	0.000230136404856698
        0.687000000000000	0.000228455237319772
        0.688000000000000	0.000226772634135404
        0.689000000000000	0.000225088604484046
        0.690000000000000	0.000223403157477219
        0.691000000000000	0.000221716302157973
        0.692000000000000	0.000220028047501674
        0.693000000000000	0.000218338402416393
        0.694000000000000	0.000216647375743749
        0.695000000000000	0.000214954976259290
        0.696000000000000	0.000213261212673287
        0.697000000000000	0.000211566093631148
        0.698000000000000	0.000209869627714152
        0.699000000000000	0.000208171823439880
        0.700000000000000	0.000206472689262857
        0.701000000000000	0.000204772233575167
        0.702000000000000	0.000203070464706869
        0.703000000000000	0.000201367390926725
        0.704000000000000	0.000199663020442559
        0.705000000000000	0.000197957361401940
        0.706000000000000	0.000196250421892621
        0.707000000000000	0.000194542209943186
        0.708000000000000	0.000192832733523400
        0.709000000000000	0.000191122000544917
        0.710000000000000	0.000189410018861620
        0.711000000000000	0.000187696796270259
        0.712000000000000	0.000185982340510835
        0.713000000000000	0.000184266659267244
        0.714000000000000	0.000182549760167654
        0.715000000000000	0.000180831650784970
        0.716000000000000	0.000179112338637446
        0.717000000000000	0.000177391831189056
        0.718000000000000	0.000175670135849977
        0.719000000000000	0.000173947259977038
        0.720000000000000	0.000172223210874334
        0.721000000000000	0.000170497995793431
        0.722000000000000	0.000168771621934047
        0.723000000000000	0.000167044096444368
        0.724000000000000	0.000165315426421531
        0.725000000000000	0.000163585618912066
        0.726000000000000	0.000161854680912299
        0.727000000000000	0.000160122619368859
        0.728000000000000	0.000158389441179028
        0.729000000000000	0.000156655153191202
        0.730000000000000	0.000154919762205294
        0.731000000000000	0.000153183274973139
        0.732000000000000	0.000151445698198861
        0.733000000000000	0.000149707038539479
        0.734000000000000	0.000147967302605023
        0.735000000000000	0.000146226496959111
        0.736000000000000	0.000144484628119325
        0.737000000000000	0.000142741702557522
        0.738000000000000	0.000140997726700347
        0.739000000000000	0.000139252706929465
        0.740000000000000	0.000137506649582046
        0.741000000000000	0.000135759560951095
        0.742000000000000	0.000134011447285867
        0.743000000000000	0.000132262314792167
        0.744000000000000	0.000130512169632767
        0.745000000000000	0.000128761017927708
        0.746000000000000	0.000127008865754783
        0.747000000000000	0.000125255719149682
        0.748000000000000	0.000123501584106541
        0.749000000000000	0.000121746466578169
        0.750000000000000	0.000119990372476459
        0.751000000000000	0.000118233307672716
        0.752000000000000	0.000116475277997892
        0.753000000000000	0.000114716289243032
        0.754000000000000	0.000112956347159639
        0.755000000000000	0.000111195457459836
        0.756000000000000	0.000109433625816850
        0.757000000000000	0.000107670857865260
        0.758000000000000	0.000105907159201256
        0.759000000000000	0.000104142535383129
        0.760000000000000	0.000102376991931340
        0.761000000000000	0.000100610534329103
        0.762000000000000	9.88431680224325e-05
        0.763000000000000	9.70748984206203e-05
        0.764000000000000	9.53057308964202e-05
        0.765000000000000	9.35356707864820e-05
        0.766000000000000	9.17647233914324e-05
        0.767000000000000	8.99928939764114e-05
        0.768000000000000	8.82201877712080e-05
        0.769000000000000	8.64466099705516e-05
        0.770000000000000	8.46721657344151e-05
        0.771000000000000	8.28968601883619e-05
        0.772000000000000	8.11206984237368e-05
        0.773000000000000	7.93436854979027e-05
        0.774000000000000	7.75658264346874e-05
        0.775000000000000	7.57871262244744e-05
        0.776000000000000	7.40075898245618e-05
        0.777000000000000	7.22272221593535e-05
        0.778000000000000	7.04460281207056e-05
        0.779000000000000	6.86640125682184e-05
        0.780000000000000	6.68811803293049e-05
        0.781000000000000	6.50975361996384e-05
        0.782000000000000	6.33130849433206e-05
        0.783000000000000	6.15278312931294e-05
        0.784000000000000	5.97417799507983e-05
        0.785000000000000	5.79549355872635e-05
        0.786000000000000	5.61673028428106e-05
        0.787000000000000	5.43788863274880e-05
        0.788000000000000	5.25896906211868e-05
        0.789000000000000	5.07997202739656e-05
        0.790000000000000	4.90089798062637e-05
        0.791000000000000	4.72174737090919e-05
        0.792000000000000	4.54252064443684e-05
        0.793000000000000	4.36321824450210e-05
        0.794000000000000	4.18384061152779e-05
        0.795000000000000	4.00438818309144e-05
        0.796000000000000	3.82486139393995e-05
        0.797000000000000	3.64526067602312e-05
        0.798000000000000	3.46558645849832e-05
        0.799000000000000	3.28583916777181e-05
        0.800000000000000	3.10601922750448e-05
        0.801000000000000	2.92612705863540e-05
        0.802000000000000	2.74616307941647e-05
        0.803000000000000	2.56612770541706e-05
        0.804000000000000	2.38602134955198e-05
        0.805000000000000	2.20584442209836e-05
        0.806000000000000	2.02559733071914e-05
        0.807000000000000	1.84528048048444e-05
        0.808000000000000	1.66489427388725e-05
        0.809000000000000	1.48443911086142e-05
        0.810000000000000	1.30391538881073e-05
        0.811000000000000	1.12332350262020e-05
        0.812000000000000	9.42663844676275e-06
        0.813000000000000	7.61936804885880e-06
        0.814000000000000	5.81142770701066e-06
        0.815000000000000	4.00282127126961e-06
        0.816000000000000	2.19355256748619e-06
        0.817000000000000	3.83625397456410e-07
        0.818000000000000	0
        0.819000000000000	0
        0.820000000000000	0
        0.821000000000000	0
        0.822000000000000	0
        0.823000000000000	0
        0.824000000000000	0
        0.825000000000000	0
        0.826000000000000	0
        0.827000000000000	0
        0.828000000000000	0
        0.829000000000000	0
        0.830000000000000	0
        0.831000000000000	0
        0.832000000000000	0
        0.833000000000000	0
        0.834000000000000	0
        0.835000000000000	0
        0.836000000000000	0
        0.837000000000000	0
        0.838000000000000	0
        0.839000000000000	0
        0.840000000000000	0
        0.841000000000000	0
        0.842000000000000	0
        0.843000000000000	0
        0.844000000000000	0
        0.845000000000000	0
        0.846000000000000	0
        0.847000000000000	0
        0.848000000000000	0
        0.849000000000000	0
        0.850000000000000	0
        0.851000000000000	0
        0.852000000000000	0
        0.853000000000000	0
        0.854000000000000	0
        0.855000000000000	0
        0.856000000000000	0
        0.857000000000000	0
        0.858000000000000	0
        0.859000000000000	0
        0.860000000000000	0
        0.861000000000000	0
        0.862000000000000	0
        0.863000000000000	0
        0.864000000000000	0
        0.865000000000000	0
        0.866000000000000	0
        0.867000000000000	0
        0.868000000000000	0
        0.869000000000000	0
        0.870000000000000	0
        0.871000000000000	0
        0.872000000000000	0
        0.873000000000000	0
        0.874000000000000	0
        0.875000000000000	0
        0.876000000000000	0
        0.877000000000000	0
        0.878000000000000	0
        0.879000000000000	0
        0.880000000000000	0
        0.881000000000000	0
        0.882000000000000	0
        0.883000000000000	0
        0.884000000000000	0
        0.885000000000000	0
        0.886000000000000	0
        0.887000000000000	0
        0.888000000000000	0
        0.889000000000000	0
        0.890000000000000	0
        0.891000000000000	0
        0.892000000000000	0
        0.893000000000000	0
        0.894000000000000	0
        0.895000000000000	0
        0.896000000000000	0
        0.897000000000000	0
        0.898000000000000	0
        0.899000000000000	0
        0.900000000000000	0
        0.901000000000000	0
        0.902000000000000	0
        0.903000000000000	0
        0.904000000000000	0
        0.905000000000000	0
        0.906000000000000	0
        0.907000000000000	0
        0.908000000000000	0
        0.909000000000000	0
        0.910000000000000	0
        0.911000000000000	0
        0.912000000000000	0
        0.913000000000000	0
        0.914000000000000	0
        0.915000000000000	0
        0.916000000000000	0
        0.917000000000000	0
        0.918000000000000	0
        0.919000000000000	0
        0.920000000000000	0
        0.921000000000000	0
        0.922000000000000	0
        0.923000000000000	0
        0.924000000000000	0
        0.925000000000000	0
        0.926000000000000	0
        0.927000000000000	0
        0.928000000000000	0
        0.929000000000000	0
        0.930000000000000	0
        0.931000000000000	0
        0.932000000000000	0
        0.933000000000000	0
        0.934000000000000	0
        0.935000000000000	0
        0.936000000000000	0
        0.937000000000000	0
        0.938000000000000	0
        0.939000000000000	0
        0.940000000000000	0
        0.941000000000000	0
        0.942000000000000	0
        0.943000000000000	0
        0.944000000000000	0
        0.945000000000000	0
        0.946000000000000	0
        0.947000000000000	0
        0.948000000000000	0
        0.949000000000000	0
        0.950000000000000	0
        0.951000000000000	0
        0.952000000000000	0
        0.953000000000000	0
        0.954000000000000	0
        0.955000000000000	0
        0.956000000000000	0
        0.957000000000000	0
        0.958000000000000	0
        0.959000000000000	0
        0.960000000000000	0
        0.961000000000000	0
        0.962000000000000	0
        0.963000000000000	0
        0.964000000000000	0
        0.965000000000000	0
        0.966000000000000	0
        0.967000000000000	0
        0.968000000000000	0
        0.969000000000000	0
        0.970000000000000	0
        0.971000000000000	0
        0.972000000000000	0
        0.973000000000000	0
        0.974000000000000	0
        0.975000000000000	0
        0.976000000000000	0
        0.977000000000000	0
        0.978000000000000	0
        0.979000000000000	0
        0.980000000000000	0
        0.981000000000000	0
        0.982000000000000	0
        0.983000000000000	0
        0.984000000000000	0
        0.985000000000000	0
        0.986000000000000	0
        0.987000000000000	0
        0.988000000000000	0
        0.989000000000000	0
        0.990000000000000	0
        0.991000000000000	0
        0.992000000000000	0
        0.993000000000000	0
        0.994000000000000	0
        0.995000000000000	0
        0.996000000000000	0
        0.997000000000000	0
        0.998000000000000	0
        0.999000000000000	0
        1	0];
    
    AlphaFit=fit(Alpha(:,1),Alpha(:,2),'linearinterp');
    varargout{1}=AlphaFit(wl);
end
end
function data=RSS(data,dim)
% Simple function to do the square root of the sum of the squares of an array or vector
if nargin<1
    dim=0;
end

switch dim
    case 0
        data=sqrt(sum(data(:).^2));
    otherwise
        data=sqrt(sum(data.^2,dim));
end
end
function adotb=MatDot(a,b)
% Author Alan Hoskins (Alan.Hoskins@Lasp.colorado.edu)
%
% Function to get dot product of 2 arrays of vectors.
%
% Function Usage:  adotb=MatDot(a,b)
% Where a is a 3xn array of column vectors and b is either a single 3x1
% column vector or a 3xn array of vectors. In either case, adotb will be a
% 1xn array where each column is the dot product of the
% corresponding columns of a and b.

if nargin<2
    b=a;
end

if size(b,2) ~= size(a,2)
    b=repmat(b,1,sizea);
end

adotb=sum(a.*b,1);

end
function n=Sellmeier(A,B,lambda)
n=zeros(size(lambda));

for I=1:3
    n=n+A(I)*lambda.^2./(lambda.^2-B(I));
end

n=sqrt(1+n);

end
function RC=Fresnel(ki,k,n12,Pol,N)

% Fresnel equation for Reflection coefficient at a surface.
% ki (3xn) =  incident k vector(s)
% k (3xn) = transmitted k vector (k = SnellsLaw(ki,N,n12))
% n12 (1xn) = ratio of index of refraction at surface n1 / n2
% Pol (1xn) = boolean indicating polarization, true is s polarized
% N (3xn) = Surface Normal and incidence location in direction of ray.

% Run SnellsLaw.m without arguments for a diagram
if nargin==1 % normal incidence un-polarized
    RC=abs((1-ki)./(1+ki)).^2;
else
    DI=n12>0; % n12 = 0 is a mirror
    if any(DI)
        RC=zeros(size(Pol(DI)));    % Reflection Coefficient
        RC(Pol)=abs((MatDot(ki(:,Pol&DI),N(:,Pol&DI))-MatDot(k(:,Pol&DI),N(:,Pol&DI)))./...
            (MatDot(ki(:,Pol&DI),N(:,Pol&DI))+MatDot(k(:,Pol&DI),N(:,Pol&DI)))).^2;
        RC(~Pol)=abs((MatDot(ki(:,~Pol&DI),N(:,~Pol&DI))./n12(~Pol&DI)-MatDot(k(:,~Pol&DI),N(:,~Pol&DI)).*n12(~Pol&DI))./...
            (MatDot(ki(:,~Pol&DI),N(:,~Pol&DI))./n12(~Pol&DI)+MatDot(k(:,~Pol&DI),N(:,~Pol&DI)).*n12(~Pol&DI))).^2;
    end
    RC(~DI)=1;
end
end
function varargout=DrawCircle(Center,Diameter,EColor,FColor,FAlpha,LWidth)
% Author: Alan Hoskins (Alan.Hoskins@LASP.Colorado.edu)
%
% DrawRect draws a circle or oval on a plot.  Note, it is not an annotation, but
% an actual patch object in the axes object.
%
% Usage handle=DrawCircle(Center,Diameter,EColor,FColor,FAlpha)
%
% The output is the handle of the patch object.  The inputs are the center
% of the rectangle [x,y,z] or [x,y] for 2D plots, the size in
% [Widthx,Widthy], the EdgeColor, FaceColor, and the Transparency.  Only
% the first 2 arguments are required.
%
t = (1/100:1/100:1)'*2*pi;

if length(Diameter)>1
    x = Diameter(1)*sin(t)/2+Center(1);
    y = Diameter(2)*cos(t)/2+Center(2);
else
    x = Diameter*sin(t)/2+Center(1);
    y = Diameter*cos(t)/2+Center(2);
end
if numel(Center)>2
    z=Center(3)*ones(size(x));
else
    z=zeros(size(x));
end

switch nargin
    case 2
        RH=fill3(x,y,z,'w','EdgeColor','k','FaceAlpha',0);
    case 3
        RH=fill3(x,y,z,FColor,'EdgeColor','k','FaceAlpha',1);
    case 4
        RH=fill3(x,y,z,FColor,'EdgeColor',EColor,'FaceAlpha',1);
    case 5
        RH=fill3(x,y,z,FColor,'EdgeColor',EColor,'FaceAlpha',FAlpha);
    case 6
        RH=fill3(x,y,z,FColor,'EdgeColor',EColor,'FaceAlpha',FAlpha,'LineWidth',LWidth);
    otherwise
        disp('Incorrect parameters passed to DrawCircle');
        disp('Correct Syntax DrawCircle(Center,Diameter,Fill Color,Edge Color,Transparency)');
end

if nargout
    varargout{1}=RH;
end
end
function varargout=DrawRect(Center,Size,Angle,EColor,FColor,FAlpha,Lwidth)
% Author: Alan Hoskins (Alan.Hoskins@LASP.Colorado.edu)
%
% DrawRect draws a rectangle on a plot.  Note, it is not an annotation, but
% an actual patch object in the axes object.
%
% Usage handle=DrawRect(Center,Size,Angle,EColor,FColor,FAlpha,Lwidth)
%
% The output is the handle of the patch object.  The inputs are the center
% of the rectangle [x,y,z] or [x,y] for 2D plots, the size in [Widthx,Widthy],
% the angle to rotate the rectangle about its center (positive angle is
% counter clockwise), the EdgeColor, FaceColor, Transparency, and
% Linewidth.  Only the first 3 arguments are required.  Because the
% rectangle is a patch and not an annotation, it can create wierd results
% when used with different axis scales.

if numel(Center)==2
    Center=[Center 0];
end

if nargin > 7 || nargin < 2
    disp('Incorrect parameters passed to DrawRect');
    disp('Correct Syntax DrawRect(Center,Size,Rotation Angle,Edge Color,Fill Color,Transparency,LineWidth)');
    disp('Where Center and Size are vectors of form [x y] and 0 <= Transparency <= 1');
    return
end

if nargin >= 3
    RotMatrix=[ cosd(Angle),-sind(Angle);...
        sind(Angle),cosd(Angle)];
    
else
    RotMatrix = eye(2);
end

LineMatrix = [-Size(1)/2,-Size(2)/2;...
    Size(1)/2,-Size(2)/2;...
    Size(1)/2,Size(2)/2;...
    -Size(1)/2,Size(2)/2;...
    -Size(1)/2,-Size(2)/2];

for Coords=1:5
    LineMatrix(Coords,:)=(RotMatrix*LineMatrix(Coords,:)')';
end

LineMatrix(:,1) = LineMatrix(:,1) + Center(1);
LineMatrix(:,2) = LineMatrix(:,2) + Center(2);

LineMatrix=[LineMatrix Center(3)*ones(5,1)];

if nargin < 4
    RH=fill3(LineMatrix(:,1),LineMatrix(:,2),LineMatrix(:,3),...
        'w','EdgeColor','k');
elseif nargin == 4
    RH=fill3(LineMatrix(:,1),LineMatrix(:,2),LineMatrix(:,3),...
        'w','EdgeColor',EColor);
elseif nargin == 5
    RH=fill3(LineMatrix(:,1),LineMatrix(:,2),LineMatrix(:,3),...
        FColor,'EdgeColor',EColor);
elseif nargin == 6
    if FAlpha<0
        RH=fill3(LineMatrix(:,1),LineMatrix(:,2),LineMatrix(:,3),...
            'w','EdgeColor',EColor,'FaceColor','none','LineWidth',2);
    else
        RH=fill3(LineMatrix(:,1),LineMatrix(:,2),LineMatrix(:,3),...
            FColor,'EdgeColor',EColor,'FaceAlpha',FAlpha,'LineWidth',0.5);
    end
elseif nargin == 7
    RH=fill3(LineMatrix(:,1),LineMatrix(:,2),LineMatrix(:,3),...
        FColor,'EdgeColor',EColor,'FaceAlpha',FAlpha,'LineWidth',Lwidth);
end

if nargout
    varargout{1}=RH;
end
end
function [n,varargout]=Quartz(wl)

if nargin==0
    wl=0.58756;
end

% This is from Zemax Birefringent library and does not appear accurate
% below 180 nm

% Ao=[2.35676495e0,-1.13996924e-2,1.08741656e-2,3.32066914e-5,1.08609346e-5,-3.10123984e-7];
% Ae=[2.38421862e0,-1.20653449e-2,1.10138430e-2,1.28863130e-4,1.68747314e-7,4.92022338e-8];
% n=[Schott(Ao,wl(:)),Schott(Ae,wl(:))];

% This is from
% Gorachand Ghosh,
% Dispersion-equation coefficients for the refractive index and birefringence of calcite and quartz crystals,
% Optics Communications,
% Volume 163, Issues 1–3,
% 1999,
% Pages 95-102,

Ao=[1.28604141-1,1.07044083,1.10202242];
Bo=[0,1.00585997e-2,100];

Ae=[1.28851804-1,1.09509924,1.15662475];
Be=[0,1.02101864e-2,100];

no=Sellmeier(Ao,Bo,wl(:));
ne=Sellmeier(Ae,Be,wl(:));

n=[no(:),ne(:)];

if nargout>1 % Absorption Coefficient in mm^-1
    % This data is pieced together from:
    
    % Gorachand Ghosh,
    % Dispersion-equation coefficients for the refractive index and birefringence of calcite and quartz crystals,
    % Optics Communications,
    % Volume 163, Issues 1–3,
    % 1999,
    % Pages 95-102,
    
    % Data thief capture from Korth website, 100 mm thick c axis
    % tranmission
    
    % Data thief capture from Korth website, 2 mm thick c axis
    % tranmission
    
    Alpha=[0.12	2.4
        0.121	2.396299931
        0.122	2.385323337
        0.123	2.367255634
        0.124	2.34228224
        0.125	2.310588571
        0.126	2.272360046
        0.127	2.22778208
        0.128	2.177040091
        0.129	2.120319497
        0.13	2.057805714
        0.131	1.98968416
        0.132	1.916140251
        0.133	1.837359406
        0.134	1.75352704
        0.135	1.664828571
        0.136	1.571449417
        0.137	1.473574994
        0.138	1.37139072
        0.139	1.265082011
        0.14	1.154834286
        0.141	1.04083296
        0.142	0.923263451
        0.143	0.802311177
        0.144	0.678161554
        0.145	0.550995432
        0.146	0.365944471
        0.147	0.284664136
        0.148	0.237446314
        0.149	0.204068902
        0.15	0.176077345
        0.151	0.140558183
        0.152	0.105748572
        0.153	0.083534886
        0.154	0.074177333
        0.155	0.066430773
        0.156	0.059183699
        0.157	0.052375698
        0.158	0.045956658
        0.159	0.039884553
        0.16	0.034051885
        0.161	0.028414184
        0.162	0.023882215
        0.163	0.022167363
        0.164	0.02130311
        0.165	0.020443131
        0.166	0.019587471
        0.167	0.01873617
        0.168	0.017889258
        0.169	0.017046761
        0.17	0.016395604
        0.171	0.016034744
        0.172	0.015733962
        0.173	0.015431631
        0.174	0.015127862
        0.175	0.014822757
        0.176	0.014516415
        0.177	0.014208926
        0.178	0.01390038
        0.179	0.013590857
        0.18	0.013280436
        0.181	0.012967602
        0.182	0.012633644
        0.183	0.012280902
        0.184	0.011927817
        0.185	0.011574443
        0.186	0.011220835
        0.187	0.010867044
        0.188	0.010513117
        0.189	0.010159102
        0.19	9.81E-03
        0.191	9.45E-03
        0.192	9.16E-03
        0.193	8.98E-03
        0.194	8.79E-03
        0.195	8.60E-03
        0.196	8.42E-03
        0.197	8.23E-03
        0.198	8.04E-03
        0.199	7.85E-03
        0.2	7.66E-03
        0.201	7.47E-03
        0.202	7.27E-03
        0.203	7.08E-03
        0.204	6.89E-03
        0.205	6.70E-03
        0.206	6.50E-03
        0.207	6.31E-03
        0.208	6.11E-03
        0.209	5.92E-03
        0.21	5.73E-03
        0.211	5.64E-03
        0.212	5.57E-03
        0.213	5.49E-03
        0.214	5.41E-03
        0.215	5.34E-03
        0.216	5.26E-03
        0.217	5.18E-03
        0.218	5.10E-03
        0.219	5.02E-03
        0.22	4.94E-03
        0.221	4.86E-03
        0.222	4.78E-03
        0.223	4.69E-03
        0.224	4.61E-03
        0.225	4.53E-03
        0.226	4.45E-03
        0.227	4.36E-03
        0.228	4.28E-03
        0.229	4.19E-03
        0.23	4.11E-03
        0.231	4.02E-03
        0.232	3.94E-03
        0.233	3.85E-03
        0.234	3.77E-03
        0.235	3.68E-03
        0.236	3.59E-03
        0.237	3.51E-03
        0.238	3.44E-03
        0.239	3.40E-03
        0.24	3.36E-03
        0.241	3.32E-03
        0.242	3.28E-03
        0.243	3.24E-03
        0.244	3.20E-03
        0.245	3.16E-03
        0.246	3.12E-03
        0.247	3.08E-03
        0.248	3.04E-03
        0.249	3.00E-03
        0.25	2.96E-03
        0.251	2.92E-03
        0.252	2.88E-03
        0.253	2.84E-03
        0.254	2.79E-03
        0.255	2.75E-03
        0.256	2.71E-03
        0.257	2.67E-03
        0.258	2.62E-03
        0.259	2.58E-03
        0.26	2.54E-03
        0.261	2.49E-03
        0.262	2.45E-03
        0.263	2.40E-03
        0.264	2.36E-03
        0.265	2.31E-03
        0.266	2.27E-03
        0.267	2.22E-03
        0.268	2.18E-03
        0.269	2.13E-03
        0.27	2.09E-03
        0.271	2.04E-03
        0.272	2.00E-03
        0.273	1.95E-03
        0.274	1.91E-03
        0.275	1.90E-03
        0.276	1.90E-03
        0.277	1.89E-03
        0.278	1.89E-03
        0.279	1.89E-03
        0.28	1.89E-03
        0.281	1.88E-03
        0.282	1.88E-03
        0.283	1.88E-03
        0.284	1.87E-03
        0.285	1.87E-03
        0.286	1.87E-03
        0.287	1.86E-03
        0.288	1.86E-03
        0.289	1.85E-03
        0.29	1.85E-03
        0.291	1.84E-03
        0.292	1.84E-03
        0.293	1.83E-03
        0.294	1.83E-03
        0.295	1.82E-03
        0.296	1.82E-03
        0.297	1.81E-03
        0.298	1.81E-03
        0.299	1.80E-03
        0.3	1.80E-03
        0.301	1.79E-03
        0.302	1.79E-03
        0.303	1.78E-03
        0.304	1.77E-03
        0.305	1.77E-03
        0.306	1.76E-03
        0.307	1.75E-03
        0.308	1.75E-03
        0.309	1.74E-03
        0.31	1.73E-03
        0.311	1.73E-03
        0.312	1.72E-03
        0.313	1.71E-03
        0.314	1.71E-03
        0.315	1.70E-03
        0.316	1.69E-03
        0.317	1.68E-03
        0.318	1.68E-03
        0.319	1.67E-03
        0.32	1.66E-03
        0.321	1.65E-03
        0.322	1.64E-03
        0.323	1.64E-03
        0.324	1.63E-03
        0.325	1.63E-03
        0.326	1.63E-03
        0.327	1.64E-03
        0.328	1.64E-03
        0.329	1.65E-03
        0.33	1.65E-03
        0.331	1.66E-03
        0.332	1.66E-03
        0.333	1.67E-03
        0.334	1.67E-03
        0.335	1.68E-03
        0.336	1.68E-03
        0.337	1.69E-03
        0.338	1.69E-03
        0.339	1.70E-03
        0.34	1.70E-03
        0.341	1.71E-03
        0.342	1.71E-03
        0.343	1.72E-03
        0.344	1.72E-03
        0.345	1.72E-03
        0.346	1.73E-03
        0.347	1.73E-03
        0.348	1.74E-03
        0.349	1.74E-03
        0.35	1.74E-03
        0.351	1.75E-03
        0.352	1.75E-03
        0.353	1.76E-03
        0.354	1.76E-03
        0.355	1.76E-03
        0.356	1.77E-03
        0.357	1.77E-03
        0.358	1.77E-03
        0.359	1.78E-03
        0.36	1.78E-03
        0.361	1.78E-03
        0.362	1.79E-03
        0.363	1.79E-03
        0.364	1.79E-03
        0.365	1.80E-03
        0.366	1.80E-03
        0.367	1.80E-03
        0.368	1.81E-03
        0.369	1.81E-03
        0.37	1.81E-03
        0.371	1.81E-03
        0.372	1.82E-03
        0.373	1.82E-03
        0.374	1.82E-03
        0.375	1.82E-03
        0.376	1.83E-03
        0.377	1.83E-03
        0.378	1.83E-03
        0.379	1.84E-03
        0.38	1.84E-03
        0.381	1.84E-03
        0.382	1.84E-03
        0.383	1.83E-03
        0.384	1.82E-03
        0.385	1.80E-03
        0.386	1.79E-03
        0.387	1.78E-03
        0.388	1.77E-03
        0.389	1.76E-03
        0.39	1.74E-03
        0.391	1.73E-03
        0.392	1.72E-03
        0.393	1.71E-03
        0.394	1.69E-03
        0.395	1.68E-03
        0.396	1.67E-03
        0.397	1.66E-03
        0.398	1.64E-03
        0.399	1.63E-03
        0.4	1.62E-03
        0.401	1.61E-03
        0.402	1.59E-03
        0.403	1.58E-03
        0.404	1.57E-03
        0.405	1.56E-03
        0.406	1.54E-03
        0.407	1.53E-03
        0.408	1.52E-03
        0.409	1.51E-03
        0.41	1.49E-03
        0.411	1.48E-03
        0.412	1.47E-03
        0.413	1.45E-03
        0.414	1.44E-03
        0.415	1.43E-03
        0.416	1.42E-03
        0.417	1.40E-03
        0.418	1.39E-03
        0.419	1.38E-03
        0.42	1.36E-03
        0.421	1.35E-03
        0.422	1.34E-03
        0.423	1.32E-03
        0.424	1.31E-03
        0.425	1.30E-03
        0.426	1.28E-03
        0.427	1.27E-03
        0.428	1.26E-03
        0.429	1.24E-03
        0.43	1.23E-03
        0.431	1.22E-03
        0.432	1.20E-03
        0.433	1.19E-03
        0.434	1.18E-03
        0.435	1.16E-03
        0.436	1.15E-03
        0.437	1.14E-03
        0.438	1.12E-03
        0.439	1.11E-03
        0.44	1.10E-03
        0.441	1.08E-03
        0.442	1.06E-03
        0.443	1.04E-03
        0.444	1.03E-03
        0.445	1.01E-03
        0.446	9.91E-04
        0.447	9.73E-04
        0.448	9.56E-04
        0.449	9.38E-04
        0.45	9.20E-04
        0.451	9.03E-04
        0.452	8.85E-04
        0.453	8.68E-04
        0.454	8.50E-04
        0.455	8.32E-04
        0.456	8.14E-04
        0.457	7.97E-04
        0.458	7.79E-04
        0.459	7.61E-04
        0.46	7.43E-04
        0.461	7.26E-04
        0.462	7.08E-04
        0.463	6.90E-04
        0.464	6.72E-04
        0.465	6.54E-04
        0.466	6.36E-04
        0.467	6.19E-04
        0.468	6.01E-04
        0.469	5.83E-04
        0.47	5.65E-04
        0.471	5.47E-04
        0.472	5.29E-04
        0.473	5.11E-04
        0.474	4.93E-04
        0.475	4.75E-04
        0.476	4.57E-04
        0.477	4.39E-04
        0.478	4.21E-04
        0.479	4.03E-04
        0.48	3.85E-04
        0.481	3.67E-04
        0.482	3.49E-04
        0.483	3.31E-04
        0.484	3.13E-04
        0.485	2.95E-04
        0.486	2.78E-04
        0.487	2.74E-04
        0.488	2.75E-04
        0.489	2.76E-04
        0.49	2.77E-04
        0.491	2.77E-04
        0.492	2.78E-04
        0.493	2.79E-04
        0.494	2.80E-04
        0.495	2.81E-04
        0.496	2.81E-04
        0.497	2.82E-04
        0.498	2.83E-04
        0.499	2.84E-04
        0.5	2.84E-04
        0.501	2.85E-04
        0.502	2.86E-04
        0.503	2.86E-04
        0.504	2.87E-04
        0.505	2.88E-04
        0.506	2.88E-04
        0.507	2.89E-04
        0.508	2.90E-04
        0.509	2.90E-04
        0.51	2.91E-04
        0.511	2.91E-04
        0.512	2.92E-04
        0.513	2.93E-04
        0.514	2.93E-04
        0.515	2.94E-04
        0.516	2.94E-04
        0.517	2.95E-04
        0.518	2.95E-04
        0.519	2.96E-04
        0.52	2.96E-04
        0.521	2.97E-04
        0.522	2.97E-04
        0.523	2.98E-04
        0.524	2.98E-04
        0.525	2.99E-04
        0.526	2.99E-04
        0.527	3.00E-04
        0.528	3.00E-04
        0.529	3.00E-04
        0.53	3.01E-04
        0.531	3.01E-04
        0.532	3.01E-04
        0.533	3.02E-04
        0.534	3.02E-04
        0.535	3.03E-04
        0.536	3.03E-04
        0.537	3.03E-04
        0.538	3.04E-04
        0.539	3.04E-04
        0.54	3.04E-04
        0.541	3.04E-04
        0.542	3.05E-04
        0.543	3.05E-04
        0.544	3.05E-04
        0.545	3.06E-04
        0.546	3.06E-04
        0.547	3.06E-04
        0.548	3.06E-04
        0.549	3.07E-04
        0.55	3.07E-04
        0.551	3.07E-04
        0.552	3.07E-04
        0.553	3.07E-04
        0.554	3.08E-04
        0.555	3.08E-04
        0.556	3.08E-04
        0.557	3.08E-04
        0.558	3.08E-04
        0.559	3.08E-04
        0.56	3.09E-04
        0.561	3.09E-04
        0.562	3.09E-04
        0.563	3.09E-04
        0.564	3.09E-04
        0.565	3.09E-04
        0.566	3.09E-04
        0.567	3.09E-04
        0.568	3.10E-04
        0.569	3.10E-04
        0.57	3.10E-04
        0.571	3.10E-04
        0.572	3.10E-04
        0.573	3.10E-04
        0.574	3.10E-04
        0.575	3.10E-04
        0.576	3.10E-04
        0.577	3.10E-04
        0.578	3.10E-04
        0.579	3.10E-04
        0.58	3.10E-04
        0.581	3.10E-04
        0.582	3.10E-04
        0.583	3.10E-04
        0.584	3.10E-04
        0.585	3.10E-04
        0.586	3.10E-04
        0.587	3.10E-04
        0.588	3.09E-04
        0.589	3.09E-04
        0.59	3.09E-04
        0.591	3.08E-04
        0.592	3.08E-04
        0.593	3.07E-04
        0.594	3.07E-04
        0.595	3.07E-04
        0.596	3.06E-04
        0.597	3.06E-04
        0.598	3.05E-04
        0.599	3.05E-04
        0.6	3.04E-04
        0.601	3.04E-04
        0.602	3.04E-04
        0.603	3.03E-04
        0.604	3.03E-04
        0.605	3.02E-04
        0.606	3.02E-04
        0.607	3.01E-04
        0.608	3.01E-04
        0.609	3.00E-04
        0.61	3.00E-04
        0.611	2.99E-04
        0.612	2.99E-04
        0.613	2.98E-04
        0.614	2.98E-04
        0.615	2.97E-04
        0.616	2.97E-04
        0.617	2.96E-04
        0.618	2.96E-04
        0.619	2.95E-04
        0.62	2.95E-04
        0.621	2.94E-04
        0.622	2.94E-04
        0.623	2.93E-04
        0.624	2.92E-04
        0.625	2.92E-04
        0.626	2.91E-04
        0.627	2.91E-04
        0.628	2.90E-04
        0.629	2.90E-04
        0.63	2.89E-04
        0.631	2.88E-04
        0.632	2.88E-04
        0.633	2.87E-04
        0.634	2.87E-04
        0.635	2.86E-04
        0.636	2.85E-04
        0.637	2.85E-04
        0.638	2.84E-04
        0.639	2.84E-04
        0.64	2.83E-04
        0.641	2.82E-04
        0.642	2.82E-04
        0.643	2.81E-04
        0.644	2.80E-04
        0.645	2.80E-04
        0.646	2.79E-04
        0.647	2.78E-04
        0.648	2.78E-04
        0.649	2.77E-04
        0.65	2.76E-04
        0.651	2.76E-04
        0.652	2.75E-04
        0.653	2.74E-04
        0.654	2.74E-04
        0.655	2.73E-04
        0.656	2.72E-04
        0.657	2.72E-04
        0.658	2.71E-04
        0.659	2.70E-04
        0.66	2.70E-04
        0.661	2.69E-04
        0.662	2.68E-04
        0.663	2.67E-04
        0.664	2.67E-04
        0.665	2.66E-04
        0.666	2.65E-04
        0.667	2.65E-04
        0.668	2.64E-04
        0.669	2.63E-04
        0.67	2.62E-04
        0.671	2.62E-04
        0.672	2.61E-04
        0.673	2.60E-04
        0.674	2.59E-04
        0.675	2.59E-04
        0.676	2.58E-04
        0.677	2.57E-04
        0.678	2.56E-04
        0.679	2.56E-04
        0.68	2.55E-04
        0.681	2.54E-04
        0.682	2.53E-04
        0.683	2.52E-04
        0.684	2.52E-04
        0.685	2.51E-04
        0.686	2.50E-04
        0.687	2.49E-04
        0.688	2.48E-04
        0.689	2.48E-04
        0.69	2.47E-04
        0.691	2.46E-04
        0.692	2.45E-04
        0.693	2.44E-04
        0.694	2.44E-04
        0.695	2.43E-04
        0.696	2.42E-04
        0.697	2.41E-04
        0.698	2.40E-04
        0.699	2.40E-04
        0.7	2.39E-04
        0.701	2.38E-04
        0.702	2.37E-04
        0.703	2.36E-04
        0.704	2.35E-04
        0.705	2.35E-04
        0.706	2.34E-04
        0.707	2.33E-04
        0.708	2.32E-04
        0.709	2.31E-04
        0.71	2.30E-04
        0.711	2.29E-04
        0.712	2.29E-04
        0.713	2.28E-04
        0.714	2.27E-04
        0.715	2.26E-04
        0.716	2.25E-04
        0.717	2.24E-04
        0.718	2.23E-04
        0.719	2.22E-04
        0.72	2.22E-04
        0.721	2.21E-04
        0.722	2.20E-04
        0.723	2.19E-04
        0.724	2.18E-04
        0.725	2.17E-04
        0.726	2.16E-04
        0.727	2.15E-04
        0.728	2.14E-04
        0.729	2.14E-04
        0.73	2.13E-04
        0.731	2.12E-04
        0.732	2.11E-04
        0.733	2.10E-04
        0.734	2.09E-04
        0.735	2.08E-04
        0.736	2.07E-04
        0.737	2.06E-04
        0.738	2.05E-04
        0.739	2.04E-04
        0.74	2.03E-04
        0.741	2.03E-04
        0.742	2.02E-04
        0.743	2.01E-04
        0.744	2.00E-04
        0.745	1.99E-04
        0.746	1.98E-04
        0.747	1.97E-04
        0.748	1.96E-04
        0.749	1.95E-04
        0.75	1.94E-04
        0.751	1.93E-04
        0.752	1.92E-04
        0.753	1.91E-04
        0.754	1.90E-04
        0.755	1.89E-04
        0.756	1.88E-04
        0.757	1.87E-04
        0.758	1.86E-04
        0.759	1.85E-04
        0.76	1.85E-04
        0.761	1.84E-04
        0.762	1.83E-04
        0.763	1.82E-04
        0.764	1.81E-04
        0.765	1.80E-04
        0.766	1.79E-04
        0.767	1.78E-04
        0.768	1.77E-04
        0.769	1.76E-04
        0.77	1.75E-04
        0.771	1.74E-04
        0.772	1.73E-04
        0.773	1.72E-04
        0.774	1.71E-04
        0.775	1.70E-04
        0.776	1.69E-04
        0.777	1.68E-04
        0.778	1.67E-04
        0.779	1.66E-04
        0.78	1.65E-04
        0.781	1.64E-04
        0.782	1.63E-04
        0.783	1.62E-04
        0.784	1.61E-04
        0.785	1.60E-04
        0.786	1.59E-04
        0.787	1.58E-04
        0.788	1.57E-04
        0.789	1.56E-04
        0.79	1.55E-04
        0.791	1.54E-04
        0.792	1.53E-04
        0.793	1.52E-04
        0.794	1.51E-04
        0.795	1.50E-04
        0.796	1.49E-04
        0.797	1.47E-04
        0.798	1.46E-04
        0.799	1.45E-04
        0.8	1.44E-04
        0.801	1.43E-04
        0.802	1.42E-04
        0.803	1.41E-04
        0.804	1.40E-04
        0.805	1.39E-04
        0.806	1.38E-04
        0.807	1.37E-04
        0.808	1.36E-04
        0.809	1.35E-04
        0.81	1.34E-04
        0.811	1.33E-04
        0.812	1.32E-04
        0.813	1.31E-04
        0.814	1.30E-04
        0.815	1.29E-04
        0.816	1.28E-04
        0.817	1.27E-04
        0.818	1.26E-04
        0.819	1.24E-04
        0.82	1.23E-04
        0.821	1.22E-04
        0.822	1.21E-04
        0.823	1.20E-04
        0.824	1.19E-04
        0.825	1.18E-04
        0.826	1.17E-04
        0.827	1.16E-04
        0.828	1.15E-04
        0.829	1.14E-04
        0.83	1.13E-04
        0.831	1.12E-04
        0.832	1.11E-04
        0.833	1.10E-04
        0.834	1.08E-04
        0.835	1.07E-04
        0.836	1.06E-04
        0.837	1.05E-04
        0.838	1.04E-04
        0.839	1.03E-04
        0.84	1.02E-04
        0.841	1.01E-04
        0.842	9.98E-05
        0.843	9.87E-05
        0.844	9.76E-05
        0.845	9.65E-05
        0.846	9.54E-05
        0.847	9.43E-05
        0.848	9.32E-05
        0.849	9.21E-05
        0.85	9.10E-05
        0.851	8.99E-05
        0.852	8.88E-05
        0.853	8.77E-05
        0.854	8.66E-05
        0.855	8.55E-05
        0.856	8.44E-05
        0.857	8.33E-05
        0.858	8.22E-05
        0.859	8.11E-05
        0.86	8.00E-05
        0.861	7.89E-05
        0.862	7.77E-05
        0.863	7.66E-05
        0.864	7.55E-05
        0.865	7.44E-05
        0.866	7.33E-05
        0.867	7.22E-05
        0.868	7.11E-05
        0.869	7.00E-05
        0.87	6.88E-05
        0.871	6.77E-05
        0.872	6.66E-05
        0.873	6.55E-05
        0.874	6.43E-05
        0.875	6.32E-05
        0.876	6.21E-05
        0.877	6.10E-05
        0.878	5.99E-05
        0.879	5.87E-05
        0.88	5.76E-05
        0.881	5.65E-05
        0.882	5.53E-05
        0.883	5.42E-05
        0.884	5.31E-05
        0.885	5.19E-05
        0.886	5.08E-05
        0.887	4.97E-05
        0.888	4.85E-05
        0.889	4.74E-05
        0.89	4.63E-05
        0.891	4.51E-05
        0.892	4.40E-05
        0.893	4.29E-05
        0.894	4.17E-05
        0.895	4.06E-05
        0.896	3.94E-05
        0.897	3.83E-05
        0.898	3.72E-05
        0.899	3.60E-05
        0.9	3.49E-05
        0.901	3.37E-05
        0.902	3.26E-05
        0.903	3.14E-05
        0.904	3.03E-05
        0.905	2.91E-05
        0.906	2.80E-05
        0.907	2.68E-05
        0.908	2.57E-05
        0.909	2.45E-05
        0.91	2.34E-05
        0.911	2.22E-05
        0.912	2.11E-05
        0.913	1.99E-05
        0.914	1.88E-05
        0.915	1.76E-05
        0.916	1.64E-05
        0.917	1.53E-05
        0.918	1.41E-05
        0.919	1.30E-05
        0.92	1.18E-05
        0.921	1.06E-05
        0.922	9.48E-06
        0.923	8.32E-06
        0.924	7.16E-06
        0.925	6.00E-06
        0.926	4.83E-06
        0.927	3.67E-06
        0.928	2.50E-06
        0.929	1.33E-06
        0.93	1.68E-07
        0.931	0.00E+00
        0.932	0.00E+00
        0.933	0.00E+00
        0.934	0.00E+00
        0.935	0.00E+00
        0.936	0.00E+00
        0.937	0
        0.938	0
        0.939	0
        0.94	0
        0.941	0
        0.942	0
        0.943	0
        0.944	0
        0.945	0
        0.946	0
        0.947	0
        0.948	0
        0.949	0
        0.95	0.00E+00
        0.951	0.00E+00
        0.952	0.00E+00
        0.953	0.00E+00
        0.954	0.00E+00
        0.955	0.00E+00
        0.956	0.00E+00
        0.957	0.00E+00
        0.958	0.00E+00
        0.959	0.00E+00
        0.96	0.00E+00
        0.961	0.00E+00
        0.962	0.00E+00
        0.963	0.00E+00
        0.964	0.00E+00
        0.965	0.00E+00
        0.966	0.00E+00
        0.967	0.00E+00
        0.968	0.00E+00
        0.969	0.00E+00
        0.97	0.00E+00
        0.971	0.00E+00
        0.972	0.00E+00
        0.973	0.00E+00
        0.974	0.00E+00
        0.975	0.00E+00
        0.976	0.00E+00
        0.977	0.00E+00
        0.978	0.00E+00
        0.979	0.00E+00
        0.98	0.00E+00
        0.981	0.00E+00
        0.982	0.00E+00
        0.983	0.00E+00
        0.984	0.00E+00
        0.985	0.00E+00
        0.986	0.00E+00
        0.987	0.00E+00
        0.988	0.00E+00
        0.989	0.00E+00
        0.99	0.00E+00
        0.991	0.00E+00
        0.992	0.00E+00
        0.993	0.00E+00
        0.994	0.00E+00
        0.995	0.00E+00
        0.996	0
        0.997	0
        0.998	0
        0.999	0
        1	0];
    AlphaFit=fit(Alpha(:,1),Alpha(:,2),'linearinterp');
    varargout{1}=AlphaFit(wl);
end
end
function [n,varargout]=SUPRASIL_3001(wl)

if nargin==0
    wl=0.58756;
end

A=[4.31179025E-001;
    6.73687852E-001;
    9.05314180E-001];
B=[1.33349413E-002;
    4.50894737E-003;
    9.92209368E+001];

n=Sellmeier(A,B,wl(:));

if nargout>1  % This gives the absorption coefficient (alpha)
    % Absorption Coefficient in mm^-1
    % Pieced together from TSIS-2 SIM measurments, and the Heraeus
    % transmission calculator.
    
    % The data below is internal transmission for 10 mm
    % HeraeusData=[0.155	0.000001
    %     0.156	0.00001
    %     0.157	0.0001
    %     0.158	0.0003
    %     0.159	0.0007
    %     0.16	0.0013
    %     0.161	0.0019
    %     0.162	0.0025
    %     0.163	0.0032
    %     0.164	0.0041
    %     0.165	0.0056
    %     0.166	0.0085
    %     0.167	0.014
    %     0.168	0.024
    %     0.169	0.0428
    %     0.17	0.074
    %     0.171	0.1224
    %     0.172	0.1893
    %     0.173	0.2748
    %     0.174	0.3723
    %     0.175	0.4741
    %     0.176	0.5714
    %     0.177	0.6592
    %     0.178	0.7325
    %     0.179	0.7911
    %     0.18	0.8359
    %     0.181	0.8698
    %     0.182	0.8945
    %     0.183	0.9137
    %     0.184	0.9272
    %     0.185	0.9387
    %     0.186	0.9466
    %     0.187	0.9528
    %     0.188	0.9558
    %     0.189	0.9587
    %     0.19	0.9603
    %     0.191	0.9632
    %     0.192	0.9651
    %     0.193	0.9706
    %     0.194	0.9735
    %     0.195	0.978
    %     0.196	0.9807
    %     0.197	0.9829
    %     0.198	0.9828
    %     0.199	0.9835
    %     0.2	0.9827
    %     0.201	0.9825
    %     0.202	0.9868
    %     0.203	0.9912
    %     0.204	0.9929
    %     0.205	0.9979
    %     0.206	0.998
    %     0.207	0.998
    %     0.208	0.9981
    %     0.209	0.9982
    %     0.21	0.9983
    %     0.212	0.9984
    %     0.214	0.9985
    %     0.216	0.9986
    %     0.218	0.9987
    %     0.22	0.9987
    %     0.222	0.9988
    %     0.224	0.9989
    %     0.226	0.9989
    %     0.228	0.9989
    %     0.234	0.999
    %     0.24	0.9988
    %     0.246	0.9986
    %     0.252	0.9988
    %     0.258	0.9992
    %     0.264	0.9995
    %     0.27	0.9996
    %     0.276	0.9997
    %     0.282	0.9998
    %     0.288	0.9998
    %     0.294	0.9999
    %     0.3	0.9999
    %     0.324	0.9998
    %     0.36	0.9999
    %     0.384	0.9999
    %     0.474	0.9999
    %     0.504	0.9998
    %     0.516	0.9998
    %     0.552	0.9999
    %     0.6	0.9999
    %     0.678	1
    %     0.726	1
    %     0.792	0.9999
    %     0.858	0.9999
    %     0.894	0.9998
    %     0.948	0.9998
    %     0.975	0.9999
    %     1	0.9999
    %     1.05	0.9998
    %     1.3	0.9998
    %     1.35	0.9999
    %     1.65	0.9999
    %     1.69	1
    %     1.775	1
    %     1.78	0.9999
    %     1.926	0.9999
    %     1.95	0.9998
    %     2.16	0.9997
    %     2.196	0.9996
    %     2.214	0.9995
    %     2.232	0.9994
    %     2.244	0.9993
    %     2.268	0.9991
    %     2.286	0.9989
    %     2.298	0.9988
    %     2.31	0.9986
    %     2.346	0.9979
    %     2.394	0.9968
    %     2.418	0.9961];
    %
    %     AlphaFit=fit(ZemaxData(ZemaxData(:,2)>0,1),log(ZemaxData(ZemaxData(:,2)>0,2))/10,'pchipinterp');
    %     AlphaFit=fit(HeraeusData(:,1),log(HeraeusData(:,2))/10,'pchipinterp');
    % Extinction Coeff in mm^1 from Heraeus Website
    
    Alpha=[0.13	3.225201821
        0.131	3.231929479
        0.132	3.237523568
        0.133	3.241257231
        0.134	3.242403611
        0.135	3.24023585
        0.136	3.234027092
        0.137	3.223050478
        0.138	3.206579153
        0.139	3.183886259
        0.14	3.154244938
        0.141	3.116928333
        0.142	3.071209588
        0.143	3.016361844
        0.144	2.951658246
        0.145	2.876371935
        0.146	2.789776054
        0.147	2.691143746
        0.148	2.579748155
        0.149	2.454862422
        0.15	2.31575969
        0.151	2.161713104
        0.152	1.991995804
        0.153	1.805880934
        0.154	1.602641637
        0.155	1.381551056
        0.156	1.151292546
        0.157	0.921034037
        0.158	0.811172808
        0.159	0.726443022
        0.16	0.664539101
        0.161	0.626590139
        0.162	0.599146455
        0.163	0.574460447
        0.164	0.549676831
        0.165	0.518498868
        0.166	0.476768912
        0.167	0.426869795
        0.168	0.372970145
        0.169	0.315121718
        0.17	0.260369019
        0.171	0.210046091
        0.172	0.166442222
        0.173	0.129171172
        0.174	0.09880553
        0.175	0.074633701
        0.176	0.055966579
        0.177	0.04167283
        0.178	0.031129194
        0.179	0.02343309
        0.18	0.017924629
        0.181	0.013949198
        0.182	0.011149038
        0.183	0.009025299
        0.184	0.007558599
        0.185	0.006325934
        0.186	0.005487866
        0.187	0.004835026
        0.188	0.004520659
        0.189	0.004217708
        0.19	0.004050954
        0.191	0.00374942
        0.192	0.003552356
        0.193	0.002984084
        0.194	0.002685745
        0.195	0.002224561
        0.196	0.001948868
        0.197	0.001724789
        0.198	0.001734964
        0.199	0.001663764
        0.2	1.75E-03
        0.201	7.18E-04
        0.202	5.75E-04
        0.203	4.58E-04
        0.204	3.55E-04
        0.205	1.93E-04
        0.206	1.03E-04
        0.207	8.43E-05
        0.208	8.20E-05
        0.209	7.91E-05
        0.21	7.63E-05
        0.211	7.35E-05
        0.212	7.09E-05
        0.213	6.84E-05
        0.214	6.59E-05
        0.215	6.36E-05
        0.216	6.15E-05
        0.217	6.01E-05
        0.218	5.88E-05
        0.219	5.75E-05
        0.22	5.62E-05
        0.221	5.49E-05
        0.222	5.37E-05
        0.223	5.25E-05
        0.224	5.13E-05
        0.225	5.02E-05
        0.226	4.91E-05
        0.227	4.80E-05
        0.228	4.70E-05
        0.229	4.65E-05
        0.23	4.61E-05
        0.231	4.57E-05
        0.232	4.53E-05
        0.233	4.48E-05
        0.234	4.44E-05
        0.235	4.42E-05
        0.236	4.55E-05
        0.237	4.70E-05
        0.238	4.85E-05
        0.239	5.01E-05
        0.24	5.18E-05
        0.241	5.36E-05
        0.242	5.54E-05
        0.243	5.66E-05
        0.244	5.76E-05
        0.245	5.85E-05
        0.246	5.95E-05
        0.247	5.93E-05
        0.248	5.84E-05
        0.249	5.75E-05
        0.25	5.66E-05
        0.251	5.49E-05
        0.252	5.22E-05
        0.253	4.96E-05
        0.254	4.71E-05
        0.255	4.48E-05
        0.256	4.25E-05
        0.257	4.00E-05
        0.258	3.74E-05
        0.259	3.51E-05
        0.26	3.28E-05
        0.261	3.07E-05
        0.262	2.88E-05
        0.263	2.69E-05
        0.264	2.53E-05
        0.265	2.42E-05
        0.266	2.31E-05
        0.267	2.20E-05
        0.268	2.11E-05
        0.269	2.02E-05
        0.27	1.98E-05
        0.271	1.95E-05
        0.272	1.91E-05
        0.273	1.87E-05
        0.274	1.84E-05
        0.275	1.80E-05
        0.276	1.77E-05
        0.277	1.73E-05
        0.278	1.70E-05
        0.279	1.67E-05
        0.28	1.64E-05
        0.281	1.61E-05
        0.282	1.58E-05
        0.283	1.55E-05
        0.284	1.52E-05
        0.285	1.49E-05
        0.286	1.46E-05
        0.287	1.43E-05
        0.288	1.41E-05
        0.289	1.38E-05
        0.29	1.35E-05
        0.291	1.33E-05
        0.292	1.31E-05
        0.293	1.29E-05
        0.294	1.28E-05
        0.295	1.26E-05
        0.296	1.25E-05
        0.297	1.23E-05
        0.298	1.22E-05
        0.299	1.20E-05
        0.3	1.19E-05
        0.301	1.18E-05
        0.302	1.16E-05
        0.303	1.15E-05
        0.304	1.14E-05
        0.305	1.12E-05
        0.306	1.11E-05
        0.307	1.10E-05
        0.308	1.09E-05
        0.309	1.07E-05
        0.31	1.06E-05
        0.311	1.05E-05
        0.312	1.04E-05
        0.313	1.02E-05
        0.314	1.01E-05
        0.315	1.00E-05
        0.316	9.90E-06
        0.317	9.79E-06
        0.318	9.67E-06
        0.319	9.56E-06
        0.32	9.45E-06
        0.321	9.35E-06
        0.322	9.24E-06
        0.323	9.13E-06
        0.324	9.03E-06
        0.325	8.92E-06
        0.326	8.82E-06
        0.327	8.72E-06
        0.328	8.62E-06
        0.329	8.52E-06
        0.33	8.42E-06
        0.331	8.33E-06
        0.332	8.22E-06
        0.333	8.12E-06
        0.334	8.02E-06
        0.335	7.93E-06
        0.336	7.83E-06
        0.337	7.73E-06
        0.338	7.64E-06
        0.339	7.55E-06
        0.34	7.45E-06
        0.341	7.36E-06
        0.342	7.27E-06
        0.343	7.19E-06
        0.344	7.10E-06
        0.345	7.01E-06
        0.346	6.93E-06
        0.347	6.84E-06
        0.348	6.74E-06
        0.349	6.54E-06
        0.35	6.35E-06
        0.351	6.16E-06
        0.352	5.97E-06
        0.353	5.79E-06
        0.354	5.62E-06
        0.355	5.45E-06
        0.356	5.28E-06
        0.357	5.13E-06
        0.358	4.97E-06
        0.359	4.82E-06
        0.36	4.68E-06
        0.361	4.50E-06
        0.362	4.33E-06
        0.363	4.16E-06
        0.364	4.00E-06
        0.365	3.84E-06
        0.366	3.69E-06
        0.367	3.55E-06
        0.368	3.41E-06
        0.369	3.08E-06
        0.37	2.73E-06
        0.371	2.42E-06
        0.372	2.01E-06
        0.373	1.64E-06
        0.374	1.44E-06
        0.375	1.39E-06
        0.376	1.50E-06
        0.377	1.73E-06
        0.378	2.01E-06
        0.379	2.39E-06
        0.38	2.84E-06
        0.381	3.25E-06
        0.382	3.71E-06
        0.383	3.85E-06
        0.384	3.92E-06
        0.385	3.95E-06
        0.386	3.90E-06
        0.387	3.86E-06
        0.388	3.81E-06
        0.389	3.76E-06
        0.39	3.71E-06
        0.391	3.61E-06
        0.392	3.52E-06
        0.393	3.42E-06
        0.394	3.33E-06
        0.395	3.24E-06
        0.396	3.24E-06
        0.397	3.28E-06
        0.398	3.31E-06
        0.399	3.35E-06
        0.4	3.39E-06
        0.401	3.42E-06
        0.402	3.44E-06
        0.403	3.45E-06
        0.404	3.46E-06
        0.405	3.47E-06
        0.406	3.48E-06
        0.407	3.49E-06
        0.408	3.46E-06
        0.409	3.42E-06
        0.41	3.38E-06
        0.411	3.33E-06
        0.412	3.29E-06
        0.413	3.25E-06
        0.414	3.22E-06
        0.415	3.23E-06
        0.416	3.24E-06
        0.417	3.25E-06
        0.418	3.26E-06
        0.419	3.27E-06
        0.42	3.32E-06
        0.421	3.39E-06
        0.422	3.47E-06
        0.423	3.55E-06
        0.424	3.63E-06
        0.425	3.71E-06
        0.426	3.74E-06
        0.427	3.76E-06
        0.428	3.77E-06
        0.429	3.79E-06
        0.43	3.81E-06
        0.431	3.83E-06
        0.432	3.84E-06
        0.433	3.78E-06
        0.434	3.71E-06
        0.435	3.65E-06
        0.436	3.59E-06
        0.437	3.53E-06
        0.438	3.50E-06
        0.439	3.52E-06
        0.44	3.53E-06
        0.441	3.55E-06
        0.442	3.57E-06
        0.443	3.59E-06
        0.444	3.67E-06
        0.445	3.82E-06
        0.446	3.98E-06
        0.447	4.14E-06
        0.448	4.31E-06
        0.449	4.48E-06
        0.45	4.64E-06
        0.451	4.75E-06
        0.452	4.87E-06
        0.453	5.00E-06
        0.454	5.12E-06
        0.455	5.25E-06
        0.456	5.38E-06
        0.457	5.42E-06
        0.458	5.44E-06
        0.459	5.46E-06
        0.46	5.48E-06
        0.461	5.50E-06
        0.462	5.52E-06
        0.463	5.54E-06
        0.464	5.57E-06
        0.465	5.59E-06
        0.466	5.61E-06
        0.467	5.63E-06
        0.468	5.65E-06
        0.469	5.61E-06
        0.47	5.52E-06
        0.471	5.44E-06
        0.472	5.35E-06
        0.473	5.26E-06
        0.474	5.18E-06
        0.475	5.10E-06
        0.476	5.03E-06
        0.477	4.97E-06
        0.478	4.91E-06
        0.479	4.85E-06
        0.48	4.85E-06
        0.481	4.91E-06
        0.482	4.97E-06
        0.483	5.03E-06
        0.484	5.09E-06
        0.485	5.20E-06
        0.486	5.52E-06
        0.487	5.85E-06
        0.488	6.21E-06
        0.489	6.59E-06
        0.49	6.94E-06
        0.491	7.20E-06
        0.492	7.47E-06
        0.493	7.74E-06
        0.494	8.03E-06
        0.495	8.32E-06
        0.496	8.55E-06
        0.497	8.78E-06
        0.498	8.90E-06
        0.499	8.90E-06
        0.5	8.90E-06
        0.501	9.27E-06
        0.502	9.29E-06
        0.503	9.31E-06
        0.504	9.33E-06
        0.505	9.35E-06
        0.506	9.37E-06
        0.507	9.25E-06
        0.508	9.14E-06
        0.509	9.02E-06
        0.51	8.91E-06
        0.511	8.79E-06
        0.512	8.68E-06
        0.513	8.57E-06
        0.514	8.47E-06
        0.515	8.36E-06
        0.516	8.25E-06
        0.517	8.15E-06
        0.518	8.05E-06
        0.519	7.94E-06
        0.52	7.85E-06
        0.521	7.75E-06
        0.522	7.65E-06
        0.523	7.56E-06
        0.524	7.47E-06
        0.525	7.38E-06
        0.526	7.29E-06
        0.527	7.20E-06
        0.528	7.11E-06
        0.529	7.02E-06
        0.53	6.94E-06
        0.531	6.85E-06
        0.532	6.77E-06
        0.533	6.68E-06
        0.534	6.60E-06
        0.535	6.52E-06
        0.536	6.44E-06
        0.537	6.36E-06
        0.538	6.28E-06
        0.539	6.23E-06
        0.54	6.19E-06
        0.541	6.15E-06
        0.542	6.11E-06
        0.543	6.07E-06
        0.544	6.04E-06
        0.545	6.00E-06
        0.546	5.96E-06
        0.547	5.92E-06
        0.548	5.88E-06
        0.549	5.85E-06
        0.55	5.81E-06
        0.551	5.77E-06
        0.552	5.74E-06
        0.553	5.70E-06
        0.554	5.67E-06
        0.555	5.64E-06
        0.556	5.71E-06
        0.557	5.78E-06
        0.558	5.85E-06
        0.559	5.93E-06
        0.56	6.00E-06
        0.561	6.08E-06
        0.562	6.15E-06
        0.563	6.23E-06
        0.564	6.31E-06
        0.565	6.39E-06
        0.566	6.46E-06
        0.567	6.52E-06
        0.568	6.59E-06
        0.569	6.66E-06
        0.57	6.63E-06
        0.571	6.57E-06
        0.572	6.52E-06
        0.573	6.46E-06
        0.574	6.41E-06
        0.575	6.35E-06
        0.576	6.30E-06
        0.577	6.26E-06
        0.578	6.29E-06
        0.579	6.33E-06
        0.58	6.36E-06
        0.581	6.40E-06
        0.582	6.44E-06
        0.583	6.47E-06
        0.584	6.49E-06
        0.585	6.45E-06
        0.586	6.40E-06
        0.587	6.36E-06
        0.588	6.32E-06
        0.589	6.27E-06
        0.59	6.23E-06
        0.591	6.19E-06
        0.592	6.15E-06
        0.593	6.11E-06
        0.594	6.06E-06
        0.595	6.02E-06
        0.596	6.05E-06
        0.597	6.09E-06
        0.598	6.13E-06
        0.599	6.18E-06
        0.6	6.22E-06
        0.601	6.27E-06
        0.602	6.31E-06
        0.603	6.36E-06
        0.604	6.41E-06
        0.605	6.45E-06
        0.606	6.49E-06
        0.607	6.44E-06
        0.608	6.38E-06
        0.609	6.33E-06
        0.61	6.27E-06
        0.611	6.21E-06
        0.612	6.15E-06
        0.613	6.10E-06
        0.614	6.04E-06
        0.615	5.98E-06
        0.616	5.93E-06
        0.617	5.88E-06
        0.618	5.82E-06
        0.619	5.77E-06
        0.62	5.71E-06
        0.621	5.66E-06
        0.622	5.61E-06
        0.623	5.56E-06
        0.624	5.51E-06
        0.625	5.46E-06
        0.626	5.41E-06
        0.627	5.36E-06
        0.628	5.31E-06
        0.629	5.26E-06
        0.63	5.21E-06
        0.631	5.16E-06
        0.632	5.12E-06
        0.633	5.07E-06
        0.634	5.02E-06
        0.635	4.96E-06
        0.636	4.90E-06
        0.637	4.83E-06
        0.638	4.77E-06
        0.639	4.71E-06
        0.64	4.64E-06
        0.641	4.58E-06
        0.642	4.52E-06
        0.643	4.46E-06
        0.644	4.41E-06
        0.645	4.35E-06
        0.646	4.29E-06
        0.647	4.24E-06
        0.648	4.18E-06
        0.649	4.12E-06
        0.65	4.07E-06
        0.651	4.02E-06
        0.652	3.96E-06
        0.653	3.91E-06
        0.654	3.86E-06
        0.655	3.81E-06
        0.656	3.76E-06
        0.657	3.71E-06
        0.658	3.66E-06
        0.659	3.62E-06
        0.66	3.56E-06
        0.661	3.43E-06
        0.662	3.31E-06
        0.663	3.18E-06
        0.664	3.06E-06
        0.665	2.95E-06
        0.666	2.84E-06
        0.667	2.73E-06
        0.668	2.63E-06
        0.669	2.53E-06
        0.67	2.45E-06
        0.671	2.37E-06
        0.672	2.30E-06
        0.673	2.23E-06
        0.674	2.16E-06
        0.675	2.10E-06
        0.676	2.03E-06
        0.677	1.99E-06
        0.678	1.98E-06
        0.679	1.97E-06
        0.68	1.96E-06
        0.681	1.96E-06
        0.682	1.95E-06
        0.683	1.94E-06
        0.684	1.94E-06
        0.685	1.93E-06
        0.686	1.92E-06
        0.687	1.92E-06
        0.688	1.91E-06
        0.689	1.90E-06
        0.69	1.90E-06
        0.691	1.89E-06
        0.692	1.89E-06
        0.693	1.85E-06
        0.694	1.81E-06
        0.695	1.77E-06
        0.696	1.72E-06
        0.697	1.68E-06
        0.698	1.64E-06
        0.699	1.60E-06
        0.7	1.56E-06
        0.701	1.52E-06
        0.702	1.46E-06
        0.703	1.39E-06
        0.704	1.33E-06
        0.705	1.28E-06
        0.706	1.22E-06
        0.707	1.17E-06
        0.708	1.12E-06
        0.709	1.11E-06
        0.71	1.15E-06
        0.711	1.20E-06
        0.712	1.25E-06
        0.713	1.31E-06
        0.714	1.38E-06
        0.715	1.46E-06
        0.716	1.51E-06
        0.717	1.51E-06
        0.718	1.51E-06
        0.719	1.51E-06
        0.72	1.51E-06
        0.721	1.48E-06
        0.722	1.42E-06
        0.723	1.35E-06
        0.724	1.29E-06
        0.725	1.23E-06
        0.726	1.16E-06
        0.727	1.08E-06
        0.728	1.01E-06
        0.729	9.41E-07
        0.73	8.78E-07
        0.731	8.23E-07
        0.732	7.73E-07
        0.733	7.27E-07
        0.734	7.16E-07
        0.735	7.68E-07
        0.736	8.24E-07
        0.737	8.84E-07
        0.738	9.47E-07
        0.739	9.98E-07
        0.74	1.05E-06
        0.741	1.10E-06
        0.742	1.15E-06
        0.743	1.21E-06
        0.744	1.26E-06
        0.745	1.29E-06
        0.746	1.32E-06
        0.747	1.35E-06
        0.748	1.39E-06
        0.749	1.42E-06
        0.75	1.46E-06
        0.751	1.49E-06
        0.752	1.54E-06
        0.753	1.60E-06
        0.754	1.66E-06
        0.755	1.72E-06
        0.756	1.79E-06
        0.757	1.86E-06
        0.758	1.93E-06
        0.759	2.00E-06
        0.76	2.08E-06
        0.761	2.16E-06
        0.762	2.24E-06
        0.763	2.33E-06
        0.764	2.41E-06
        0.765	2.51E-06
        0.766	2.59E-06
        0.767	2.64E-06
        0.768	2.69E-06
        0.769	2.74E-06
        0.77	2.79E-06
        0.771	2.84E-06
        0.772	2.89E-06
        0.773	2.94E-06
        0.774	3.00E-06
        0.775	3.05E-06
        0.776	3.11E-06
        0.777	3.17E-06
        0.778	3.23E-06
        0.779	3.29E-06
        0.78	3.35E-06
        0.781	3.39E-06
        0.782	3.39E-06
        0.783	3.39E-06
        0.784	3.40E-06
        0.785	3.40E-06
        0.786	3.40E-06
        0.787	3.41E-06
        0.788	3.41E-06
        0.789	3.41E-06
        0.79	3.42E-06
        0.791	3.42E-06
        0.792	3.42E-06
        0.793	3.43E-06
        0.794	3.43E-06
        0.795	3.45E-06
        0.796	3.50E-06
        0.797	3.55E-06
        0.798	3.61E-06
        0.799	3.66E-06
        0.8	3.71E-06
        0.801	3.77E-06
        0.802	3.82E-06
        0.803	3.88E-06
        0.804	3.88E-06
        0.805	3.82E-06
        0.806	3.75E-06
        0.807	3.69E-06
        0.808	3.64E-06
        0.809	3.58E-06
        0.81	3.52E-06
        0.811	3.47E-06
        0.812	3.43E-06
        0.813	3.38E-06
        0.814	3.34E-06
        0.815	3.30E-06
        0.816	3.27E-06
        0.817	3.34E-06
        0.818	3.45E-06
        0.819	3.56E-06
        0.82	3.67E-06
        0.821	3.79E-06
        0.822	3.86E-06
        0.823	3.85E-06
        0.824	3.84E-06
        0.825	3.84E-06
        0.826	3.83E-06
        0.827	3.82E-06
        0.828	3.82E-06
        0.829	3.81E-06
        0.83	3.75E-06
        0.831	3.68E-06
        0.832	3.61E-06
        0.833	3.54E-06
        0.834	3.52E-06
        0.835	3.52E-06
        0.836	3.52E-06
        0.837	3.52E-06
        0.838	3.52E-06
        0.839	3.51E-06
        0.84	3.51E-06
        0.841	3.51E-06
        0.842	3.51E-06
        0.843	3.51E-06
        0.844	3.51E-06
        0.845	3.51E-06
        0.846	3.50E-06
        0.847	3.50E-06
        0.848	3.50E-06
        0.849	3.50E-06
        0.85	3.50E-06
        0.851	3.50E-06
        0.852	3.50E-06
        0.853	3.50E-06
        0.854	3.49E-06
        0.855	3.49E-06
        0.856	3.49E-06
        0.857	3.49E-06
        0.858	3.49E-06
        0.859	3.49E-06
        0.86	3.49E-06
        0.861	3.48E-06
        0.862	3.48E-06
        0.863	3.48E-06
        0.864	3.48E-06
        0.865	3.48E-06
        0.866	3.48E-06
        0.867	3.48E-06
        0.868	3.66E-06
        0.869	3.94E-06
        0.87	4.29E-06
        0.871	5.50E-06
        0.872	6.27E-06
        0.873	6.31E-06
        0.874	6.35E-06
        0.875	6.38E-06
        0.876	6.42E-06
        0.877	6.46E-06
        0.878	6.49E-06
        0.879	6.53E-06
        0.88	6.57E-06
        0.881	6.61E-06
        0.882	6.64E-06
        0.883	6.68E-06
        0.884	6.72E-06
        0.885	6.76E-06
        0.886	6.80E-06
        0.887	6.84E-06
        0.888	6.92E-06
        0.889	7.01E-06
        0.89	7.11E-06
        0.891	7.21E-06
        0.892	7.31E-06
        0.893	7.41E-06
        0.894	7.51E-06
        0.895	7.61E-06
        0.896	7.72E-06
        0.897	7.81E-06
        0.898	7.82E-06
        0.899	7.83E-06
        0.9	7.85E-06
        0.901	7.86E-06
        0.902	7.87E-06
        0.903	7.89E-06
        0.904	7.90E-06
        0.905	7.91E-06
        0.906	7.92E-06
        0.907	7.94E-06
        0.908	7.95E-06
        0.909	7.96E-06
        0.91	7.98E-06
        0.911	7.99E-06
        0.912	8.00E-06
        0.913	8.02E-06
        0.914	8.03E-06
        0.915	8.04E-06
        0.916	8.06E-06
        0.917	8.07E-06
        0.918	8.08E-06
        0.919	8.10E-06
        0.92	8.11E-06
        0.921	8.08E-06
        0.922	8.02E-06
        0.923	7.96E-06
        0.924	7.90E-06
        0.925	7.84E-06
        0.926	7.78E-06
        0.927	7.73E-06
        0.928	7.67E-06
        0.929	7.62E-06
        0.93	7.57E-06
        0.931	7.54E-06
        0.932	7.50E-06
        0.933	7.47E-06
        0.934	7.43E-06
        0.935	7.39E-06
        0.936	7.36E-06
        0.937	7.32E-06
        0.938	7.29E-06
        0.939	7.25E-06
        0.94	7.22E-06
        0.941	7.19E-06
        0.942	7.15E-06
        0.943	7.12E-06
        0.944	7.09E-06
        0.945	7.06E-06
        0.946	7.03E-06
        0.947	7.00E-06
        0.948	6.97E-06
        0.949	6.94E-06
        0.95	6.91E-06
        0.951	6.88E-06
        0.952	6.85E-06
        0.953	6.82E-06
        0.954	6.79E-06
        0.955	6.76E-06
        0.956	6.73E-06
        0.957	6.70E-06
        0.958	6.66E-06
        0.959	6.59E-06
        0.96	6.52E-06
        0.961	6.46E-06
        0.962	6.39E-06
        0.963	6.33E-06
        0.964	6.26E-06
        0.965	6.20E-06
        0.966	6.13E-06
        0.967	6.10E-06
        0.968	6.11E-06
        0.969	6.13E-06
        0.97	6.14E-06
        0.971	6.16E-06
        0.972	6.18E-06
        0.973	6.19E-06
        0.974	6.21E-06
        0.975	6.22E-06
        0.976	6.24E-06
        0.977	6.25E-06
        0.978	6.22E-06
        0.979	6.20E-06
        0.98	6.17E-06
        0.981	6.15E-06
        0.982	6.12E-06
        0.983	6.10E-06
        0.984	6.07E-06
        0.985	6.05E-06
        0.986	6.02E-06
        0.987	6.00E-06
        0.988	5.97E-06
        0.989	5.95E-06
        0.99	5.93E-06
        0.991	5.93E-06
        0.992	5.93E-06
        0.993	5.93E-06
        0.994	5.93E-06
        0.995	5.93E-06
        0.996	5.98E-06
        0.997	6.05E-06
        0.998	6.12E-06
        0.999	6.18E-06
        1	6.22E-06];
    AlphaFit=fit(Alpha(:,1),Alpha(:,2),'linearinterp');
    varargout{1}=AlphaFit(wl);
end
end
function [n,varargout]=Spectrosil_2000(wl)

if nargin==0
    wl=0.58756;
end

A=[4.2666235E-001;
    6.77390924E-001;
    8.74516460E-001];
B=[1.33495666E-002;
    4.54610259E-003;
    9.55752206E+001];

n=Sellmeier(A,B,wl(:));

if nargout>1
    % Absoprtion Coeff in mm^1 from Heraeus Website
    Alpha=[0.125	3.998631207
        0.126	3.99578315
        0.127	3.99120363
        0.128	3.984518351
        0.129	3.975353013
        0.13	3.963333319
        0.131	3.94808497
        0.132	3.929233669
        0.133	3.906405117
        0.134	3.879225017
        0.135	3.847319069
        0.136	3.810312977
        0.137	3.767832441
        0.138	3.719503164
        0.139	3.664950848
        0.14	3.603801194
        0.141	3.535679905
        0.142	3.460212681
        0.143	3.377025226
        0.144	3.285743241
        0.145	3.185992428
        0.146	3.077398488
        0.147	2.959587124
        0.148	2.832184038
        0.149	2.694814931
        0.15	2.547105505
        0.151	2.388681463
        0.152	2.219168505
        0.153	2.038192335
        0.154	1.845378653
        0.155	1.640353162
        0.156	1.422741563
        0.157	1.192169559
        0.158	0.954992586
        0.159	0.765000945
        0.16	0.606955529
        0.161	0.478090861
        0.162	0.376585862
        0.163	0.296631714
        0.164	0.223673061
        0.165	0.154380182
        0.166	0.106734969
        0.167	0.07411126
        0.168	0.049258877
        0.169	0.02970529
        0.17	0.018450246
        0.171	0.011672692
        0.172	0.007395881
        0.173	0.005264986
        0.174	0.003825343
        0.175	0.002779353
        0.176	0.002019375
        0.177	0.001468035
        0.178	0.001072721
        0.179	0.000786145
        0.18	0.000607152
        0.181	0.000468913
        0.182	0.000380626
        0.183	0.000320917
        0.184	0.000270574
        0.185	0.000235413
        0.186	0.000207387
        0.187	0.000182697
        0.188	0.000167467
        0.189	0.000156882
        0.19	0.000146966
        0.191	0.000137677
        0.192	0.000124939
        0.193	0.00011232
        0.194	0.000102308
        0.195	9.34E-05
        0.196	8.54E-05
        0.197	8.50E-05
        0.198	8.50E-05
        0.199	8.46E-05
        0.2	8.35E-05
        0.201	8.12E-05
        0.202	7.88E-05
        0.203	7.65E-05
        0.204	7.42E-05
        0.205	7.20E-05
        0.206	6.99E-05
        0.207	6.78E-05
        0.208	6.58E-05
        0.209	6.39E-05
        0.21	6.20E-05
        0.211	6.01E-05
        0.212	5.83E-05
        0.213	5.66E-05
        0.214	5.49E-05
        0.215	5.32E-05
        0.216	5.16E-05
        0.217	5.01E-05
        0.218	4.86E-05
        0.219	4.72E-05
        0.22	4.57E-05
        0.221	4.44E-05
        0.222	4.31E-05
        0.223	4.18E-05
        0.224	4.05E-05
        0.225	3.93E-05
        0.226	3.84E-05
        0.227	3.76E-05
        0.228	3.68E-05
        0.229	3.61E-05
        0.23	3.53E-05
        0.231	3.46E-05
        0.232	3.39E-05
        0.233	3.32E-05
        0.234	3.25E-05
        0.235	3.18E-05
        0.236	3.12E-05
        0.237	3.05E-05
        0.238	2.99E-05
        0.239	2.93E-05
        0.24	2.87E-05
        0.241	2.81E-05
        0.242	2.75E-05
        0.243	2.69E-05
        0.244	2.64E-05
        0.245	2.58E-05
        0.246	2.53E-05
        0.247	2.48E-05
        0.248	2.43E-05
        0.249	2.38E-05
        0.25	2.33E-05
        0.251	2.28E-05
        0.252	2.23E-05
        0.253	2.19E-05
        0.254	2.14E-05
        0.255	2.10E-05
        0.256	2.05E-05
        0.257	2.01E-05
        0.258	1.97E-05
        0.259	1.93E-05
        0.26	1.88E-05
        0.261	1.84E-05
        0.262	1.80E-05
        0.263	1.76E-05
        0.264	1.72E-05
        0.265	1.68E-05
        0.266	1.64E-05
        0.267	1.60E-05
        0.268	1.57E-05
        0.269	1.53E-05
        0.27	1.50E-05
        0.271	1.46E-05
        0.272	1.43E-05
        0.273	1.40E-05
        0.274	1.36E-05
        0.275	1.33E-05
        0.276	1.30E-05
        0.277	1.27E-05
        0.278	1.24E-05
        0.279	1.22E-05
        0.28	1.19E-05
        0.281	1.16E-05
        0.282	1.14E-05
        0.283	1.12E-05
        0.284	1.10E-05
        0.285	1.08E-05
        0.286	1.06E-05
        0.287	1.04E-05
        0.288	1.02E-05
        0.289	9.98E-06
        0.29	9.79E-06
        0.291	9.60E-06
        0.292	9.42E-06
        0.293	9.24E-06
        0.294	9.07E-06
        0.295	8.90E-06
        0.296	8.79E-06
        0.297	8.70E-06
        0.298	8.61E-06
        0.299	8.52E-06
        0.3	8.43E-06
        0.301	8.34E-06
        0.302	8.25E-06
        0.303	8.16E-06
        0.304	8.07E-06
        0.305	7.99E-06
        0.306	7.90E-06
        0.307	7.82E-06
        0.308	7.74E-06
        0.309	7.65E-06
        0.31	7.57E-06
        0.311	7.49E-06
        0.312	7.41E-06
        0.313	7.33E-06
        0.314	7.25E-06
        0.315	7.17E-06
        0.316	7.09E-06
        0.317	7.70E-06
        0.318	8.60E-06
        0.319	7.93E-06
        0.32	7.13E-06
        0.321	6.94E-06
        0.322	6.81E-06
        0.323	6.68E-06
        0.324	6.55E-06
        0.325	6.42E-06
        0.326	6.31E-06
        0.327	6.28E-06
        0.328	6.26E-06
        0.329	6.23E-06
        0.33	6.21E-06
        0.331	6.18E-06
        0.332	6.16E-06
        0.333	6.13E-06
        0.334	6.11E-06
        0.335	6.08E-06
        0.336	6.06E-06
        0.337	6.03E-06
        0.338	6.01E-06
        0.339	5.98E-06
        0.34	5.96E-06
        0.341	5.93E-06
        0.342	5.91E-06
        0.343	5.89E-06
        0.344	5.86E-06
        0.345	5.84E-06
        0.346	5.81E-06
        0.347	5.77E-06
        0.348	5.73E-06
        0.349	5.70E-06
        0.35	5.66E-06
        0.351	5.62E-06
        0.352	5.58E-06
        0.353	5.54E-06
        0.354	5.50E-06
        0.355	5.47E-06
        0.356	5.43E-06
        0.357	5.39E-06
        0.358	5.36E-06
        0.359	5.32E-06
        0.36	5.28E-06
        0.361	5.11E-06
        0.362	4.91E-06
        0.363	4.72E-06
        0.364	4.53E-06
        0.365	4.36E-06
        0.366	4.19E-06
        0.367	4.03E-06
        0.368	3.87E-06
        0.369	3.72E-06
        0.37	3.57E-06
        0.371	3.40E-06
        0.372	3.10E-06
        0.373	2.83E-06
        0.374	2.58E-06
        0.375	2.35E-06
        0.376	2.01E-06
        0.377	2.59E-06
        0.378	3.34E-06
        0.379	3.72E-06
        0.38	4.09E-06
        0.381	4.51E-06
        0.382	4.50E-06
        0.383	4.31E-06
        0.384	4.14E-06
        0.385	4.08E-06
        0.386	4.03E-06
        0.387	3.97E-06
        0.388	3.92E-06
        0.389	3.87E-06
        0.39	3.82E-06
        0.391	3.77E-06
        0.392	3.82E-06
        0.393	3.94E-06
        0.394	4.06E-06
        0.395	4.18E-06
        0.396	4.31E-06
        0.397	4.44E-06
        0.398	4.57E-06
        0.399	4.68E-06
        0.4	4.38E-06
        0.401	4.08E-06
        0.402	3.81E-06
        0.403	3.55E-06
        0.404	3.48E-06
        0.405	3.52E-06
        0.406	3.56E-06
        0.407	3.60E-06
        0.408	3.65E-06
        0.409	3.69E-06
        0.41	3.73E-06
        0.411	3.77E-06
        0.412	3.81E-06
        0.413	3.83E-06
        0.414	3.85E-06
        0.415	3.87E-06
        0.416	3.89E-06
        0.417	3.91E-06
        0.418	3.93E-06
        0.419	3.94E-06
        0.42	3.96E-06
        0.421	3.98E-06
        0.422	4.00E-06
        0.423	4.05E-06
        0.424	4.10E-06
        0.425	4.14E-06
        0.426	4.19E-06
        0.427	4.24E-06
        0.428	4.29E-06
        0.429	4.34E-06
        0.43	4.39E-06
        0.431	4.45E-06
        0.432	4.50E-06
        0.433	4.55E-06
        0.434	4.60E-06
        0.435	4.63E-06
        0.436	4.62E-06
        0.437	4.60E-06
        0.438	4.59E-06
        0.439	4.58E-06
        0.44	4.57E-06
        0.441	4.57E-06
        0.442	4.64E-06
        0.443	4.72E-06
        0.444	4.79E-06
        0.445	4.87E-06
        0.446	4.94E-06
        0.447	5.02E-06
        0.448	5.10E-06
        0.449	5.18E-06
        0.45	5.27E-06
        0.451	5.35E-06
        0.452	5.43E-06
        0.453	5.52E-06
        0.454	5.61E-06
        0.455	5.70E-06
        0.456	5.74E-06
        0.457	5.76E-06
        0.458	5.78E-06
        0.459	5.81E-06
        0.46	5.83E-06
        0.461	5.85E-06
        0.462	5.88E-06
        0.463	5.90E-06
        0.464	5.92E-06
        0.465	5.95E-06
        0.466	5.97E-06
        0.467	5.99E-06
        0.468	6.02E-06
        0.469	6.04E-06
        0.47	6.07E-06
        0.471	6.09E-06
        0.472	6.11E-06
        0.473	6.14E-06
        0.474	6.16E-06
        0.475	6.19E-06
        0.476	6.24E-06
        0.477	6.43E-06
        0.478	6.64E-06
        0.479	6.84E-06
        0.48	7.06E-06
        0.481	7.28E-06
        0.482	7.51E-06
        0.483	7.75E-06
        0.484	7.99E-06
        0.485	8.22E-06
        0.486	8.36E-06
        0.487	8.51E-06
        0.488	8.66E-06
        0.489	8.82E-06
        0.49	8.97E-06
        0.491	9.12E-06
        0.492	9.25E-06
        0.493	9.39E-06
        0.494	9.53E-06
        0.495	9.50E-06
        0.496	9.13E-06
        0.497	8.77E-06
        0.498	8.48E-06
        0.499	8.67E-06
        0.5	8.81E-06
        0.501	8.91E-06
        0.502	8.93E-06
        0.503	8.94E-06
        0.504	8.95E-06
        0.505	8.97E-06
        0.506	8.98E-06
        0.507	8.99E-06
        0.508	9.01E-06
        0.509	8.92E-06
        0.51	8.68E-06
        0.511	8.46E-06
        0.512	8.24E-06
        0.513	8.02E-06
        0.514	7.81E-06
        0.515	7.61E-06
        0.516	7.41E-06
        0.517	7.22E-06
        0.518	7.03E-06
        0.519	6.85E-06
        0.52	6.79E-06
        0.521	6.81E-06
        0.522	6.83E-06
        0.523	6.85E-06
        0.524	6.87E-06
        0.525	6.90E-06
        0.526	6.92E-06
        0.527	6.94E-06
        0.528	6.96E-06
        0.529	6.99E-06
        0.53	7.02E-06
        0.531	7.05E-06
        0.532	7.07E-06
        0.533	7.10E-06
        0.534	7.13E-06
        0.535	7.16E-06
        0.536	7.19E-06
        0.537	7.22E-06
        0.538	7.28E-06
        0.539	7.40E-06
        0.54	7.51E-06
        0.541	7.62E-06
        0.542	7.74E-06
        0.543	7.86E-06
        0.544	7.98E-06
        0.545	8.10E-06
        0.546	8.22E-06
        0.547	8.35E-06
        0.548	8.48E-06
        0.549	8.52E-06
        0.55	8.41E-06
        0.551	8.30E-06
        0.552	8.20E-06
        0.553	8.09E-06
        0.554	8.06E-06
        0.555	8.18E-06
        0.556	8.31E-06
        0.557	8.44E-06
        0.558	8.57E-06
        0.559	8.71E-06
        0.56	8.84E-06
        0.561	8.98E-06
        0.562	9.12E-06
        0.563	9.27E-06
        0.564	9.41E-06
        0.565	9.56E-06
        0.566	9.71E-06
        0.567	9.86E-06
        0.568	9.84E-06
        0.569	9.81E-06
        0.57	9.78E-06
        0.571	9.75E-06
        0.572	9.72E-06
        0.573	9.69E-06
        0.574	9.66E-06
        0.575	9.63E-06
        0.576	9.60E-06
        0.577	9.58E-06
        0.578	9.55E-06
        0.579	9.52E-06
        0.58	9.51E-06
        0.581	9.55E-06
        0.582	9.60E-06
        0.583	9.65E-06
        0.584	9.70E-06
        0.585	9.75E-06
        0.586	9.80E-06
        0.587	9.85E-06
        0.588	9.90E-06
        0.589	9.95E-06
        0.59	1.00E-05
        0.591	1.01E-05
        0.592	1.01E-05
        0.593	1.02E-05
        0.594	1.02E-05
        0.595	1.03E-05
        0.596	1.03E-05
        0.597	1.02E-05
        0.598	1.02E-05
        0.599	1.02E-05
        0.6	1.02E-05
        0.601	1.02E-05
        0.602	1.02E-05
        0.603	1.02E-05
        0.604	1.01E-05
        0.605	1.01E-05
        0.606	1.01E-05
        0.607	1.01E-05
        0.608	1.01E-05
        0.609	1.01E-05
        0.61	1.01E-05
        0.611	1.01E-05
        0.612	1.00E-05
        0.613	1.00E-05
        0.614	1.00E-05
        0.615	1.00E-05
        0.616	9.94E-06
        0.617	9.87E-06
        0.618	9.81E-06
        0.619	9.74E-06
        0.62	9.68E-06
        0.621	9.61E-06
        0.622	9.55E-06
        0.623	9.49E-06
        0.624	9.42E-06
        0.625	9.36E-06
        0.626	9.30E-06
        0.627	9.24E-06
        0.628	9.18E-06
        0.629	9.08E-06
        0.63	8.90E-06
        0.631	8.72E-06
        0.632	8.54E-06
        0.633	8.37E-06
        0.634	8.20E-06
        0.635	8.04E-06
        0.636	7.87E-06
        0.637	7.71E-06
        0.638	7.67E-06
        0.639	7.80E-06
        0.64	7.94E-06
        0.641	8.09E-06
        0.642	8.23E-06
        0.643	8.32E-06
        0.644	8.26E-06
        0.645	8.19E-06
        0.646	8.13E-06
        0.647	8.07E-06
        0.648	8.01E-06
        0.649	7.94E-06
        0.65	7.88E-06
        0.651	7.82E-06
        0.652	7.76E-06
        0.653	7.70E-06
        0.654	7.64E-06
        0.655	7.59E-06
        0.656	7.54E-06
        0.657	7.50E-06
        0.658	7.45E-06
        0.659	7.41E-06
        0.66	7.36E-06
        0.661	7.32E-06
        0.662	7.27E-06
        0.663	7.23E-06
        0.664	7.18E-06
        0.665	7.14E-06
        0.666	7.10E-06
        0.667	7.05E-06
        0.668	7.03E-06
        0.669	7.02E-06
        0.67	7.00E-06
        0.671	6.99E-06
        0.672	6.98E-06
        0.673	6.96E-06
        0.674	6.95E-06
        0.675	6.93E-06
        0.676	6.92E-06
        0.677	6.91E-06
        0.678	6.89E-06
        0.679	6.88E-06
        0.68	6.86E-06
        0.681	6.85E-06
        0.682	6.84E-06
        0.683	6.82E-06
        0.684	6.81E-06
        0.685	6.79E-06
        0.686	6.78E-06
        0.687	6.77E-06
        0.688	6.75E-06
        0.689	6.74E-06
        0.69	6.72E-06
        0.691	6.71E-06
        0.692	6.70E-06
        0.693	6.68E-06
        0.694	6.67E-06
        0.695	6.66E-06
        0.696	6.65E-06
        0.697	6.63E-06
        0.698	6.62E-06
        0.699	6.61E-06
        0.7	6.60E-06
        0.701	6.58E-06
        0.702	6.57E-06
        0.703	6.56E-06
        0.704	6.55E-06
        0.705	6.53E-06
        0.706	6.52E-06
        0.707	6.51E-06
        0.708	6.50E-06
        0.709	6.49E-06
        0.71	6.47E-06
        0.711	6.46E-06
        0.712	6.45E-06
        0.713	6.44E-06
        0.714	6.64E-06
        0.715	6.94E-06
        0.716	7.26E-06
        0.717	7.59E-06
        0.718	7.93E-06
        0.719	8.29E-06
        0.72	8.67E-06
        0.721	9.41E-06
        0.722	1.08E-05
        0.723	1.21E-05
        0.724	1.25E-05
        0.725	1.27E-05
        0.726	1.30E-05
        0.727	1.32E-05
        0.728	1.27E-05
        0.729	1.22E-05
        0.73	1.16E-05
        0.731	1.12E-05
        0.732	1.08E-05
        0.733	1.05E-05
        0.734	1.02E-05
        0.735	9.96E-06
        0.736	9.69E-06
        0.737	9.43E-06
        0.738	9.17E-06
        0.739	8.92E-06
        0.74	8.68E-06
        0.741	8.44E-06
        0.742	8.21E-06
        0.743	7.99E-06
        0.744	7.76E-06
        0.745	7.51E-06
        0.746	7.27E-06
        0.747	7.03E-06
        0.748	6.80E-06
        0.749	6.58E-06
        0.75	6.36E-06
        0.751	6.16E-06
        0.752	5.96E-06
        0.753	5.76E-06
        0.754	5.59E-06
        0.755	5.52E-06
        0.756	5.47E-06
        0.757	5.42E-06
        0.758	5.37E-06
        0.759	5.32E-06
        0.76	5.27E-06
        0.761	5.22E-06
        0.762	5.17E-06
        0.763	5.12E-06
        0.764	5.08E-06
        0.765	5.06E-06
        0.766	5.04E-06
        0.767	5.03E-06
        0.768	5.01E-06
        0.769	4.99E-06
        0.77	4.95E-06
        0.771	4.70E-06
        0.772	4.39E-06
        0.773	4.13E-06
        0.774	4.15E-06
        0.775	4.23E-06
        0.776	4.30E-06
        0.777	4.38E-06
        0.778	4.46E-06
        0.779	4.53E-06
        0.78	4.61E-06
        0.781	4.70E-06
        0.782	4.78E-06
        0.783	4.86E-06
        0.784	4.94E-06
        0.785	4.92E-06
        0.786	4.88E-06
        0.787	4.83E-06
        0.788	4.79E-06
        0.789	4.74E-06
        0.79	4.70E-06
        0.791	4.65E-06
        0.792	4.61E-06
        0.793	4.56E-06
        0.794	4.43E-06
        0.795	4.11E-06
        0.796	3.81E-06
        0.797	3.83E-06
        0.798	3.97E-06
        0.799	4.12E-06
        0.8	4.27E-06
        0.801	4.43E-06
        0.802	4.59E-06
        0.803	4.75E-06
        0.804	4.91E-06
        0.805	4.80E-06
        0.806	4.61E-06
        0.807	4.43E-06
        0.808	4.25E-06
        0.809	4.08E-06
        0.81	3.92E-06
        0.811	3.82E-06
        0.812	3.80E-06
        0.813	3.77E-06
        0.814	3.75E-06
        0.815	3.73E-06
        0.816	3.71E-06
        0.817	3.69E-06
        0.818	3.67E-06
        0.819	3.65E-06
        0.82	3.62E-06
        0.821	3.60E-06
        0.822	3.58E-06
        0.823	3.56E-06
        0.824	3.52E-06
        0.825	3.36E-06
        0.826	3.20E-06
        0.827	3.04E-06
        0.828	2.90E-06
        0.829	2.76E-06
        0.83	2.62E-06
        0.831	2.36E-06
        0.832	2.09E-06
        0.833	2.13E-06
        0.834	2.26E-06
        0.835	2.40E-06
        0.836	2.54E-06
        0.837	2.70E-06
        0.838	2.86E-06
        0.839	3.01E-06
        0.84	3.14E-06
        0.841	3.23E-06
        0.842	2.96E-06
        0.843	2.65E-06
        0.844	2.38E-06
        0.845	2.25E-06
        0.846	2.16E-06
        0.847	2.08E-06
        0.848	2.02E-06
        0.849	1.99E-06
        0.85	1.96E-06
        0.851	1.93E-06
        0.852	1.90E-06
        0.853	1.87E-06
        0.854	1.85E-06
        0.855	1.82E-06
        0.856	1.79E-06
        0.857	1.77E-06
        0.858	1.73E-06
        0.859	1.69E-06
        0.86	1.65E-06
        0.861	1.61E-06
        0.862	1.66E-06
        0.863	1.96E-06
        0.864	2.31E-06
        0.865	2.73E-06
        0.866	3.10E-06
        0.867	2.73E-06
        0.868	2.40E-06
        0.869	2.31E-06
        0.87	2.37E-06
        0.871	2.44E-06
        0.872	2.50E-06
        0.873	2.57E-06
        0.874	2.63E-06
        0.875	2.70E-06
        0.876	2.77E-06
        0.877	2.84E-06
        0.878	2.92E-06
        0.879	2.99E-06
        0.88	3.07E-06
        0.881	3.10E-06
        0.882	3.08E-06
        0.883	3.06E-06
        0.884	3.03E-06
        0.885	3.01E-06
        0.886	2.99E-06
        0.887	2.97E-06
        0.888	2.92E-06
        0.889	2.74E-06
        0.89	2.57E-06
        0.891	2.42E-06
        0.892	2.27E-06
        0.893	2.13E-06
        0.894	2.12E-06
        0.895	2.25E-06
        0.896	2.40E-06
        0.897	2.55E-06
        0.898	2.72E-06
        0.899	2.88E-06
        0.9	2.98E-06
        0.901	3.07E-06
        0.902	3.04E-06
        0.903	2.57E-06
        0.904	2.17E-06
        0.905	2.04E-06
        0.906	2.02E-06
        0.907	2.01E-06
        0.908	1.99E-06
        0.909	1.97E-06
        0.91	1.95E-06
        0.911	1.93E-06
        0.912	1.90E-06
        0.913	1.87E-06
        0.914	1.85E-06
        0.915	1.82E-06
        0.916	1.80E-06
        0.917	1.77E-06
        0.918	1.74E-06
        0.919	1.67E-06
        0.92	1.61E-06
        0.921	1.58E-06
        0.922	2.03E-06
        0.923	2.81E-06
        0.924	3.88E-06
        0.925	5.35E-06
        0.926	7.39E-06
        0.927	1.02E-05
        0.928	1.41E-05
        0.929	1.94E-05
        0.93	2.42E-05
        0.931	2.78E-05
        0.932	3.20E-05
        0.933	3.68E-05
        0.934	4.23E-05
        0.935	4.87E-05
        0.936	5.60E-05
        0.937	6.44E-05
        0.938	7.40E-05
        0.939	8.21E-05
        0.94	8.93E-05
        0.941	9.71E-05
        0.942	0.00010556
        0.943	0.000114783
        0.944	0.000124736
        0.945	0.000129254
        0.946	0.00013001
        0.947	0.00013077
        0.948	0.000131535
        0.949	0.000131065
        0.95	0.000124971
        0.951	0.00011916
        0.952	0.000113619
        0.953	0.000108336
        0.954	0.000103299
        0.955	9.85E-05
        0.956	9.39E-05
        0.957	8.95E-05
        0.958	8.54E-05
        0.959	8.14E-05
        0.96	7.76E-05
        0.961	7.40E-05
        0.962	7.06E-05
        0.963	6.73E-05
        0.964	6.41E-05
        0.965	6.06E-05
        0.966	5.72E-05
        0.967	5.40E-05
        0.968	5.10E-05
        0.969	4.81E-05
        0.97	4.54E-05
        0.971	4.29E-05
        0.972	4.05E-05
        0.973	3.82E-05
        0.974	3.63E-05
        0.975	3.47E-05
        0.976	3.31E-05
        0.977	3.15E-05
        0.978	3.01E-05
        0.979	2.87E-05
        0.98	2.74E-05
        0.981	2.61E-05
        0.982	2.49E-05
        0.983	2.38E-05
        0.984	2.30E-05
        0.985	2.25E-05
        0.986	2.21E-05
        0.987	2.16E-05
        0.988	2.12E-05
        0.989	2.08E-05
        0.99	2.03E-05
        0.991	1.99E-05
        0.992	1.95E-05
        0.993	1.91E-05
        0.994	1.87E-05
        0.995	1.84E-05
        0.996	1.78E-05
        0.997	1.70E-05
        0.998	1.63E-05
        0.999	1.56E-05
        1	1.53E-05];
    AlphaFit=fit(Alpha(:,1),Alpha(:,2),'linearinterp');
    varargout{1}=AlphaFit(wl);
end % This gives the absorption coefficient (alpha)
end
function DrawArrow(x,y,c,hs)

if nargin <4
    hs=0.1; % Fraction head size compared to line
end
if nargin <3
    c='k';
end

line(x,y,'color',c)
LVec=hs*[x(1)-x(2),y(1)-y(2),0]';
AHV1 = RM(3,30)*LVec;
AHV2 = RM(3,-30)*LVec;

line([x(2);x(2)+AHV1(1)],[y(2);y(2)+AHV1(2)],'color',c)
line([x(2);x(2)+AHV2(1)],[y(2);y(2)+AHV2(2)],'color',c)
end
function [DI,FI]=FindPeaks1D(Data,LevelVec,HalfWidth,Method,FW,xWin)

% Brute force somewhat kludgy method of finding peaks in a 1D vector.
% Gives position and FWHM in the DI (Row) dimension.  FWHM
% only works if Method is 1 or 2 (centroid or 2D Gaussian Fit)

if nargin<4
    Method=0; % 0 is Max, 1 Centroid, 2 is Gaussian Fit
end

if nargin<5
    FW=0.5;
end

DI=[];

for Lvl=LevelVec
    DIt=find(Data(:)>Lvl);
    if ~isempty(DIt)
        DI=[DI; DIt];
        I=1;
        while I<length(DI)
            Dist=sqrt((DI(I+1:end)-repelem(DI(I),length(DI)-I)').^2)>HalfWidth;
            DI=DI([true(I,1); Dist]);
            I=I+1;
        end
    end
end

if nargin==6
    CenIndex=DI>min(xWin)&DI<max(xWin);
    DI=DI(CenIndex);
end

FI=zeros(size(DI));
II=1:length(Data);

if Method
    for I=1:length(DI)
        IIndex=(-HalfWidth:HalfWidth)+DI(I);
        IIndex=IIndex(IIndex>0&IIndex<size(Data,1));
        ITemp=Data(IIndex);
        ITemp=(ITemp-min(ITemp))/(max(ITemp)-min(ITemp));
        switch Method
            case 1
                DI(I)=sum(ITemp(:)'.*II(IIndex))/sum(ITemp);
                Y=min(II(IIndex)):0.01:max(II(IIndex));
                YFit=fit(II(IIndex)',ITemp(:),'linear');
                PeakValue=max(YFit(Y));
                FI(I)=Y(find(YFit(Y)>FW*PeakValue,1,'last'))-Y(find(YFit(Y)>FW*PeakValue,1,'first'));
            case 2
                try
                    IFit=fit(II(IIndex)',ITemp,'gauss1');
                    DI(I)=IFit.b1;
                    FI(I)=2*IFit.c1*sqrt(-log(1/2));
                catch
                    DI(I)=sum(ITemp.*II(IIndex)')/sum(ITemp);
                    Y=min(II(IIndex)):0.01:max(II(IIndex));
                    YFit=fit(II(IIndex),ITemp,'linear');
                    PeakValue=max(YFit(Y));
                    FI(I)=Y(find(YFit(Y)>FW*PeakValue,1,'last'))-Y(find(YFit(Y)>FW*PeakValue,1,'first'));
                end
        end
    end
end
end
function RotMat=RM(V,t)

if length(V)>1
    RotMat=[cosd(t)+(1-cosd(t))*V(1)^2,V(1)*V(2)*(1-cosd(t))-V(3)*sind(t),V(1)*V(3)*(1-cosd(t))+V(2)*sind(t);
        V(1)*V(2)*(1-cosd(t))+V(3)*sind(t),cosd(t)+(1-cosd(t))*V(2)^2,V(2)*V(3)*(1-cosd(t))-V(1)*sind(t);
        V(1)*V(3)*(1-cosd(t))-V(2)*sind(t),V(2)*V(3)*(1-cosd(t))+V(1)*sind(t),cosd(t)+(1-cosd(t))*V(3)^2];
else
    switch V
        case 1
            RotMat=[1 0 0;
                0 cosd(t) -sind(t);
                0 sind(t) cosd(t)];
        case 2
            RotMat=[cosd(t) 0 sind(t);
                0 1 0;
                -sind(t) 0 cosd(t)];
        case 3
            RotMat=[cosd(t) -sind(t) 0;
                sind(t) cosd(t) 0;
                0 0 1];
        otherwise
            RotMat=eye(3);
    end
end
end
function FilteredData=LowPass(Data,pts,Method)

% This is similar to a box car routine except it multiplies in frequency
% space with a smoothed rectangle rather than a sinc function.  The variable
% pts is comparable to the width of the box in the boxcar function.

% Method 1 is smothed Boxcar in frequency space, Method 2 is Gaussian (both

if nargin<3
    Method=1;
end

faxis=linspace(-1/2,1/2,length(Data));

if Method==1
    % Smoothed box car, order 2 is a standard hanning window
    PW=PolyWindow(faxis(1),-1.1/(2*pts),-1/(2*pts),1/(2*pts),1.1/(2*pts),faxis(end),length(faxis),2)';
elseif Method==2
    W=-(pts)^2/log(0.5)*2;
    PW=exp(-faxis.^2*W);
end

FFTData=fftshift(fft(Data));
FilteredData=abs(ifft(FFTData(:).*PW(:)));

if false % Plot Frequency Space
    figure(107)
    clf
    plot(faxis,20*log10(abs(FFTData)),':k',faxis,20*log10(abs(FFTData.*PW)),'-b')
    axis('tight')
end

if false % Plot Real Space
    figure(108)
    clf
    plot(1:length(Data),Data,':k',1:length(Data),FilteredData,'-b')
    axis('tight')
end
end
function Dock
set(gcf,'WindowStyle','docked')
end