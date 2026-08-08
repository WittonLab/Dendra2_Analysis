function [IIxyzCrop,cropRegion,roiCrop] = dendra2_StackCrop(pat,fnPrefix,dayVec,saveMode,ROI)
% ---------------------------------------------------------------------- % 
% Function to align in Z and crop confocal image stacks for dendra2 
% imaging across days.
% 
% Inputs:
% pat      - parent directory of registered .tif image
% fn       - file name prefix (omitting the day number and '_SeriesReg.tif'
% dayVec   - (1,m) numerical vector of image days (i.e., [1, 1.1, 3, 7, 14])
% saveMode - mode for saving registered output to parent directory (0 = off; 1 = on) 
% ROI      - (optional) ROI information - structure output from dendra2_SeriesReg_dft.m 
%
% Outputs:
% IIxyzCrop    - (1,k) cell array containing crooped image series (in XYCZT format), across days (in k)
% cropRegion   - structure with info about crop region (position in XY and number of Z slices)
% roiCrop      - structure with info about user defined ROIs (from dendra2_SeriesReg.m) in aligned and croped images  
%
% by Jon Witton, 24 July 2026
% ---------------------------------------------------------------------- %

% --- Format inputs ---------------------------------------------------- %
if ~strcmp(pat(end),'\')                                                    % if parent directory name doesn't end with a backslash
    pat = [pat '\'];                                                        % add one
end

if nargin < 4, saveMode = 0; end                                            % default save mode is off

if nargin < 3, dayVec = [1, 1.1, 3, 7, 14]; end                             % default imaging days

nFiles = length(dayVec);                                                    % number of files to algin


% --- Define some constants -------------------------------------------- %
% image channels
tPMTchl = 2;                                                                % transmission PMT image channel 
GFPchl = 1;                                                                 % GFP image channel
RFPchl = 3;                                                                 % RFP image channel

% crop amounts
XYcropSize = 20;                                                            % size of area to crop along X and Y dimensions (in microns)

nCropSlicesMax = 10;                                                        % maximum number of Z-slices in cropped array
 
% --- Region Of Interest (ROI) data from previous dentra2_SliceReg ----- %
if ~exist('ROI','var')                                                      % if ROI data hasn't been input into function
    load([pat 'ROI.mat'],'ROI')                                             % load from parent folder 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- READ DATA -------------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Read each image file into a common cell array using Bioformats --- %
II = cell(1,nFiles);                                                        % pre-allocate cell array for all image files

xPixScale = nan(1,nFiles);                                                  % pre-allocate array of X-pixel scale 
yPixScale = nan(1,nFiles);                                                  % pre-allocate array of Y-pixel scale
zPixScale = nan(1,nFiles);                                                  % pre-allocate array of Z-pixel scale

for k = 1:nFiles
    disp(['Reading day ' num2str(dayVec(k)) '...'])                         % user update
    
    % --- Series registered .tif --- %
    rdr = BioformatsImage([pat fnPrefix num2str(dayVec(k)) '_SeriesReg.tif']);
    
    % image properties
    rowPix = rdr.height;                                                    % number of row (Y) pixels
    colPix = rdr.height;                                                    % number of column (X) pixels
    nChls = rdr.sizeC;                                                      % number of channels
    nSlices = rdr.sizeZ;                                                    % number of slices (in Z)
    nFrames = rdr.sizeT;                                                    % number of frames (in time)

    if ~(nFrames==1)
        error('Image contains more than one frame, i.e. timepoint')         % CATCH
    end
    
    % generate image array
    I = zeros(rowPix,colPix,nChls,nSlices,'uint16');                        % pre-allocate 16-bit 4D image array

    for cc = 1:nChls                                                        
        for ss = 1:nSlices
            I(:,:,cc,ss) = getPlane(rdr,ss,cc,1);                           % read image and add to preallocated 4D array (YXZCT format)
        end
    end
    
    % add 4D image array to parrent cell array
    II{k} = I;                                                              


    % --- metadata from original image file --- %
    rdrMeta = bfGetReader([pat fnPrefix num2str(dayVec(k)) '.ome.tif']);    % metadata reader

    meta = rdrMeta.getMetadataStore();                                      % get metadata
    
    px = meta.getPixelsPhysicalSizeX(0);                                    % X-pixel scale metadata
    xPixScale(k) = px.value().doubleValue();                                % X-pixel scale (in microns)

    py = meta.getPixelsPhysicalSizeY(0);                                    % Y-pixel scale metadata
    yPixScale(k) = py.value().doubleValue();                                % Y-pixel scale (in microns)

    pz = meta.getPixelsPhysicalSizeZ(0);                                    % Z-pixel scale metadata
    zPixScale(k) = pz.value().doubleValue();                                % Z-pixel scale (in microns)    
end


% --- CATCH: Check X,Y,Z pixel scales are consistent across all images - %
% --- (at least to an accuracy of 0.01 microns) ------------------------ %
disp('Checking pixel scale is consistent across files...')                  % user update

xPixScale2dp = (round(xPixScale*100))/100;                                  % round X-pixel scale to to decimal places (i.e., 0.01 of a micron)
yPixScale2dp = (round(yPixScale*100))/100;                                  % round Y-pixel scale to to decimal places (i.e., 0.01 of a micron)
zPixScale2dp = (round(zPixScale*100))/100;                                  % round Z-pixel scale to to decimal places (i.e., 0.01 of a micron)

if length(unique(xPixScale2dp))~=1                                          % CATCH in X
    error('X pixel scale is not consistent across images')
elseif length(unique(yPixScale2dp))~=1                                      % CATCH in Y
    error('Y pixel scale is not consistent across images')
elseif length(unique(zPixScale2dp))~=1                                      % CATCH in Z
    error('Z pixel scale is not consistent across images')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- CROP IMAGES ------------------------------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Cropping images to align X,Y, and Z axis positions across days...')   % user update

% --- Get start and end indices to align image Z-slices across days ---- %
Zdiff = ROI.ZprofilesPkDiffMn;                                              % Z-axis differences across days

startZslices = Zdiff + abs(min(Zdiff)) + 1;                                 % Z start slices for cropped image series

nCropSlices = nSlices - max(startZslices) + 1;                              % number of slices in cropped image
if nCropSlices > nCropSlicesMax                                             % if number of slices exceeds user defined limit
    nCropSlices = nCropSlicesMax;                                           % set number of slices in cropped image to the user defined limit
end

stopZslices = startZslices + nCropSlices - 1;                               % Z stop slices for cropped image series


% Centre cropped Z stack as close as possible to centre of the parent stack
halfwidthFull = ceil(nSlices/2);                                            % half-width of parent stack
halfwidthCrop = ceil(nCropSlices/2);                                        % half-width of cropped stack
midZslice = startZslices + halfwidthCrop;                                   % Z slices in middle position of cropped stack
midZdiff = midZslice-halfwidthFull;                                         % difference between middle of cropped stack and middle of parent stack

if midZdiff(1)>0 && min(startZslices)>1 && max(stopZslices)<nSlices         % if middle Z slice in crop is ABOVE mid Z slice in parent stack, and start and end points are not at the extreme edges of parent stack 
    
    while midZdiff(1)>0 && min(startZslices)>1 && max(stopZslices)<nSlices  % iteratively...
        startZslices = startZslices - 1;                                    % subtract 1 from start Z slice in crop...
        stopZslices = stopZslices - 1;                                      % ...and substract 1 from stop Z slice in crop 
        midZslice = startZslices + halfwidthCrop;                           % Recalculate mid Z slice in crop
        midZdiff = midZslice - halfwidthFull;                               % Recalculate difference between mid Z slice in crop and mid Z slice in parent stack 
    end

elseif midZdiff(1)<0 && min(startZslices)>1 && max(stopZslices)<nSlices     % if middle Z slice in crop is BELOW mid Z slice in parent stack, and start and end points are not at the extreme edges of parent stack

    while midZdiff(1)<0 && min(startZslices)>1 && max(stopZslices)<nSlices  % iteratively...
        startZslices = startZslices + 1;                                    % add 1 to start Z slice in crop...
        stopZslices = stopZslices + 1;                                      % ...and add 1 to stop Z slice in crop 
        midZslice = startZslices + halfwidthCrop;                           % Recalculate mid Z slice in crop
        midZdiff = midZslice - halfwidthFull;                               % Recalculate difference between mid Z slice in crop and mid Z slice in parent stack
    end

end


% --- Gets start and end indices to crop X and Y axes across days ------ %
% X-axis 
Xpad = round(XYcropSize*unique(xPixScale2dp));                              % number of X-pixels to be cropped at either end of the array 
nCropXpix = round(colPix-(2*Xpad));                                         % number of X pixels in the cropped image

startX = Xpad+1;                                                            % start index for cropping X-axis pixels    
stopX = startX + nCropXpix - 1;                                             % stop index for cropping X-axis pixels

% Y-axis
Ypad = round(XYcropSize*unique(yPixScale2dp));                              % number of Y-pixels to be cropped at either end of the array 
nCropYpix = round(rowPix-(2*Ypad));                                         % number of Y pixels in the cropped image

startY = Ypad+1;                                                            % start index for cropping Y-axis pixels
stopY = startY + nCropYpix - 1;                                             % stop index for cropping Y-axis pixels


% --- Crop each image to standardise X,Y,Z positions across days ------- %
IIxyzCrop = cell(1,nFiles);                                                 % pre-allocate cropped image array                                      

for k = 1:nFiles
    IIxyzCrop{k} = II{k}(startX:stopX, startY:stopY, :, startZslices(k):stopZslices(k)); % crop image
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------- MAKE SUMMARY FIGURES TO CHECK OUTPUT ----------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- t-PMT max projection with overlaid XY crop region  --------------- % 
H1 = figure('Name','t-PMT and GFP max projections with XY crop region','WindowState','maximized');

cropRegionXvec = [startX,startX,stopX,stopX,startX];                        % vector of X-axis line coordinates for cropped region of interest  
cropRegionYvec = [startY,stopY,stopY,startY,startY];                        % vector of Y-axis line coordinates for cropped region of interest

for k = 1:nFiles
    
    % transmission PMT channel
    Ipmt = squeeze(II{k}(:,:,tPMTchl,:));                                   % t-PMT channel for kth image 
    IpmtMx = uint16(max(single(Ipmt),[],3));                                % Max projection of t-PMT channel for kth image 
    
    axPMT = subplot(2,nFiles,k);                                            % make subplot for t-PMT image
    imagesc(IpmtMx)                                                         % display t-PMT image        
    axis image
    colormap(axPMT,gray)                                                    % grayscale colormap                        
    title(['Day ' num2str(dayVec(k)) ' - tPMT, Max Projection'])
    hold on
    plot(cropRegionXvec,cropRegionYvec,'--g','linewidth',1)                 % overlay line showing cropped region

    % GFP channel
    Igfp = squeeze(II{k}(:,:,GFPchl,:));                                    % GFP channel for kth image 
    IgfpMx = uint16(max(single(Igfp),[],3));                                % Max projection of GFP channel for kth image
    
    axGFP = subplot(2,nFiles,nFiles+k);                                     % make subplot for GFP image
    imagesc(IgfpMx)                                                         % display GFP image 
    axis image
    colormap(axGFP,parula)                                                  % parula colormap
    title(['Day ' num2str(dayVec(k)) ' - GFP, Max Projection'])
    hold on
    plot(cropRegionXvec,cropRegionYvec,'--g','linewidth',1)                 % overlay line showing cropped region

    % Overlay user selected ROIs
    for j = 1:numel(ROI.roiCoords)                                          % cycle through ROIs

        roi_Xvec = [ROI.roiCoords{j}(1,1), ROI.roiCoords{j}(1,1), ROI.roiCoords{j}(2,1), ROI.roiCoords{j}(2,1), ROI.roiCoords{j}(1,1)]; % vector of X-axis line coordinates for jth user ROI 
        roi_Yvec = [ROI.roiCoords{j}(1,2), ROI.roiCoords{j}(2,2), ROI.roiCoords{j}(2,2), ROI.roiCoords{j}(1,2), ROI.roiCoords{j}(1,2)]; % vector of Y-axis line coordinates for jth user ROI 

        plot(axPMT,roi_Xvec,roi_Yvec,'m')                                   % overlay line showing jth user ROI on t-PMT image
        plot(axGFP,roi_Xvec,roi_Yvec,'m')                                   % overlay line showing jth user ROI on GFP image
    end

end


% --- For each user ROI, plot X,Y,Z profile after crop of parent image - %

% First need to adjust ROI coordinates for cropped image series 
nROIs = numel(ROI.roiCoords);                                               % number of user defined ROIs (loaded from dendra2_seriesReg.m)

roiCoordsCrop = cell(size(ROI.roiCoords));                                  % pre-allocate array for user defined ROI coordinates in cropped images 
roiCoordsCrop_Xvec = cell(size(ROI.roiCoords));                             % pre-allocate array for X-axis user ROI line coordinates in cropped images 
roiCoordsCrop_Yvec = cell(size(ROI.roiCoords));                             % pre-allocate array for Y-axis user ROI line coordinates in cropped images  

for j = 1:nROIs                                                             % cycle through each ROI

    coords = ROI.roiCoords{j};                                              % ROI coordinates, based on original image scale
    cropCoords(:,1) = coords(:,1) - Xpad;                                   % X-axis ROI coordinates, based on cropped image scale
    cropCoords(:,2) = coords(:,2) - Ypad;                                   % Y-axis ROI coordinates, based on cropped image scale

    roiCoordsCrop{j} = cropCoords;                                          % add ROI coordinates for cropped image scale to output array
    
    roiCoordsCrop_Xvec{j} = [cropCoords(1,1), cropCoords(1,1), cropCoords(2,1), cropCoords(2,1), cropCoords(1,1)]; % vector of X-axis line coordinates for user ROI, based on cropped image scale
    roiCoordsCrop_Yvec{j} = [cropCoords(1,2), cropCoords(2,2), cropCoords(2,2), cropCoords(1,2), cropCoords(1,2)]; % vector of Y-axis line coordinates for user ROI, based on cropped image scale

end

% --- Cropped images with ROI overlay --- %
H2 = figure('Name','Cropped images with user ROI overlay','WindowState','maximized'); % make figure

% Custom green and red colormaps
bits = 2^16;                                                                % image bitdepth - 16-bit (uint16)
cmapGreen = [zeros(bits, 1), linspace(0, 1, bits)', zeros(bits, 1)];        % custom green colormap (for GFP channel)
cmapRed =   [linspace(0, 1, bits)', zeros(bits, 1), zeros(bits, 1)];        % custom red colormap (for RFP channel)

for k = 1:nFiles                                                            % cylce through files in image series
    
    % GFP channel
    IgfpCrop = squeeze(IIxyzCrop{k}(:,:,GFPchl,:));                         % cropped GFP channel for kth image in series
    IgfpCropMn = uint16(mean(single(IgfpCrop),3));                          % mean projection of cropped GFP channel

    axGFP = subplot(2,nFiles,k);                                            % make axis for GFP image
    imagesc(IgfpCropMn)                                                     % display mean projection of cropped GFP channel
    axis image
    colormap(axGFP,cmapGreen)                                               % apply green colormap
    colorbar                                                                % add a colorbar (needed as images of different days have different signal, but displayed images are scaled for visualisation
    title(['Day ' num2str(dayVec(k)) ' - GFP, Mean Projection'])
    hold on

    % RFP channel
    IrfpCrop = squeeze(IIxyzCrop{k}(:,:,RFPchl,:));                         % cropped RFP channel for kth image in series
    IrfpCropSum = uint16(mean(single(IrfpCrop),3));                         % mean projection of cropped RFP channel

    axRFP = subplot(2,nFiles,nFiles+k);                                     % make axis for RFP image

    imagesc(IrfpCropSum)                                                    % display mean projection of cropped RFP channel
    axis image
    colormap(axRFP,cmapRed)                                                 % apply green colormap
    colorbar                                                                % add a colorbar (needed as images of different days have different signal, but displayed images are scaled for visualisation
    title(['Day ' num2str(dayVec(k)) ' - RFP, Mean Projection'])
    hold on

    % Overlay ROIs
    for j = 1:nROIs                                                         % cyle through user defined ROIs
        plot(axGFP,roiCoordsCrop_Xvec{j},roiCoordsCrop_Yvec{j},'m')         % overlay line showing jth user ROI on cropped GFP image
        plot(axRFP,roiCoordsCrop_Xvec{j},roiCoordsCrop_Yvec{j},'m')         % overlay line showing jth user ROI on cropped RFP image
    end

end


% --- Display user ROIs with Z-axis profile post crop of parent image -- %
IIroiCrop = cell(nROIs,nFiles);                                             % pre-allocate cell array for user defined ROIs (in rows) in cropped images, across days (in columns) 
IIroiCropMnXY = cell(nROIs,nFiles);                                         % pre-allocate cell array for mean XY projections of user defined ROIs (in rows) in cropped images, across days (in columns)
IIroiCropMnXZ = cell(nROIs,nFiles);                                         % pre-allocate cell array for mean XZ projections of user defined ROIs (in rows) in cropped images, across days (in columns)
IIroiCropMnYZ = cell(nROIs,nFiles);                                         % pre-allocate cell array for mean YZ projections of user defined ROIs (in rows) in cropped images, across days (in columns)

H3 = cell(nROIs,1);                                                         % pre-allocate figure handle array

for j = 1:nROIs                                                             % cycle through user ROIs
    
    H3{j} = figure('Name',['ROI' num2str(j) '_post crop'],'WindowState','maximized'); % make figure

    for k = 1:nFiles                                                        % cycle through image files

        IgfpCrop = squeeze(IIxyzCrop{k}(:,:,GFPchl,:));                     % GFP channel in cropped image series for kth day 
        
        % --- get ROI --- %
        IIroiCrop{j,k} = IgfpCrop(roiCoordsCrop{j}(1,2):roiCoordsCrop{j}(2,2), roiCoordsCrop{j}(1,1):roiCoordsCrop{j}(2,1), :); % jth user ROI in cropped image 
        
        % XY projection
        IIroiCropMnXY{j,k} = uint16(mean(single(IIroiCrop{j,k}),3));        % mean XY projection of jth ROI 
        
        % XZ projection
        tmpXZproj = mean(single(IIroiCrop{j,k}),1);                         % mean XZ projection of jth ROI
        tmpXZproj_RotFlip = flipud(rot90(squeeze(tmpXZproj)));              % need to squeeze the array to convert from 3D to 2D, then rotate and flip to put it in the correct orientation
        IIroiCropMnXZ{j,k} = uint16(tmpXZproj_RotFlip);                     % finally, convert the array from 32-bit back to 16-bit

        % YZ projection
        tmpYZproj = mean(single(IIroiCrop{j,k}),2);                         % mean XZ projection of jth ROI
        tmpYZproj_RotFlip = flipud(rot90(squeeze(tmpYZproj)));              % need to squeeze the array to convert from 3D to 2D, then rotate and flip to put it in the correct orientation
        IIroiCropMnYZ{j,k} = uint16(tmpYZproj_RotFlip);                     % finally, convert the array from 32-bit back to 16-bit

        % --- display ROI --- %
        % XY projection
        subplot(3,nFiles,k);                                                % XY projection displayed on top row
        imagesc(IIroiCropMnXY{j,k})                                         % plot XY projection image of jth ROI
        axis image
        title(['Day ' num2str(dayVec(k)) ' - XY Mean Projection'])
        ylabel('Y pixels')
        xlabel('X pixels')

        % XZ projection
        subplot(3,nFiles,nFiles+k);                                         % XZ projection displayed on middle row
        imagesc(IIroiCropMnXZ{j,k})                                         % plot XZ projection image of jth ROI
        axis image
        title(['Day ' num2str(dayVec(k)) ' - XZ Mean Projection'])
        ylabel('Z slices')
        xlabel('X pixels')

        % YZ projection
        subplot(3,nFiles,2*nFiles+k);                                       % YZ projection displayed on bottom row
        imagesc(IIroiCropMnYZ{j,k})                                         % plot YZ projection image of jth ROI
        axis image
        title(['Day ' num2str(dayVec(k)) ' - YZ Mean Projection'])
        ylabel('Z slices')
        xlabel('Y pixels')
    end
end


% --- User ROI Z-axis profiles in cropped images ----------------------- %
H4 = figure('Name','User ROI Z-axis profiles for cropped image series','WindowState','maximized'); % make figure

IIroiCrop_Zprofile = cell(nROIs,nFiles);                                    % pre-allocate ROI Z-axis profiles cell array [(j,k) format, where ROI number (j) is in the rows and file number (k) is in the columns] 
IIroiCrop_ZprofileNorm = cell(nROIs,nFiles);                                % pre-allocate normalised ROI Z-axis profiles cell array [(j,k) format, where ROI number (j) is in the rows and file number (k) is in the columns] 
IIroiCrop_ZprofilePk = nan(nROIs,nFiles);                                   % pre-allocate array for number of the ROI Z-axis slice with peak signal

for k = 1:nFiles                                                            % cycle through image files
    for j = 1:nROIs                                                         % cycle through user ROIs

        % --- calculate Z-axis profile --- %
        profile = nan(nCropSlices,1);                                       % pre-allocate temporary Z-axis profile array     
        
        for ss = 1:nCropSlices                                              % cycle through slices  
            tmpIroiCrop = IIroiCrop{j,k}(:,:,ss);                           % current (ssth) frame, of jth ROI, for kth image file
            profile(ss) = mean(double(tmpIroiCrop(:)));                     % find mean intensity of current frame and add to Z axis profile    
        end
        
        smProfile = smooth(profile,3);                                      % smooth profile with 3 point moving average
        
        minSubProfile = smProfile-min(smProfile);                           % subtract minimum from profile
        pkNormProfile = minSubProfile/max(minSubProfile);                   % divide profile by maximum to normalise intensity between 0 and 1
    
        IIroiCrop_Zprofile{j,k} = smProfile;                                % Z-profile output                               
        IIroiCrop_ZprofileNorm{j,k} = pkNormProfile;                        % Normalised Z-profile output  
    
        [~,pkSliceInd] = max(pkNormProfile);                                % index of Z-slice corresponding to peak ROI intensity
        IIroiCrop_ZprofilePk(j,k) = pkSliceInd;                             % peak intensity Z-profile slice output  
    
        % --- display Z-axis profile --- %
        subplot(1,nROIs,j)
        hold on
        plot(IIroiCrop_ZprofileNorm{j,k},1:nCropSlices)                     % plot Z-axis profile
        set(gca,'ydir','reverse')                                           % set Y-axis direction as 'reverse', so slice number increases from top to bottom of the figure
        set(gca,'ylim',[1 nCropSlices])                                     % set Y-axis limits
        title(['ROI ' num2str(j)])
        xlabel('Normalised intensity')
        ylabel('Z-slice number')
        legend('location','northeast')   
    end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------- MAKE AND SAVE OUTPUTS ------------------------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Organise (some of the) output variables -------------------------- %
% Crop region coordinates
cropRegion.Xcoords = [min(cropRegionXvec) max(cropRegionXvec)];
cropRegion.Xvec = cropRegionXvec;
cropRegion.Ycoords = [min(cropRegionYvec) max(cropRegionYvec)];
cropRegion.Yvec = cropRegionYvec;
cropRegion.nSlices = nCropSlices; 

% User ROIs in cropped image series
roiCrop.Coords           = roiCoordsCrop;

roiCrop.II               = IIroiCrop;

roiCrop.IIMnProj.XY      = IIroiCropMnXY; 
roiCrop.IIMnProj.XZ      = IIroiCropMnXZ;
roiCrop.IIMnProj.YZ      = IIroiCropMnYZ;

roiCrop.Zprofile.Raw     = IIroiCrop_Zprofile;
roiCrop.Zprofile.Norm    = IIroiCrop_ZprofileNorm;
roiCrop.Zprofile.PkSlice = IIroiCrop_ZprofilePk;


% --- Save output variables and figures (to parent directory) ---------- %
if saveMode == 1
    
    % save cropped images 
    disp('Saving cropped images to parent directory...');                   % user update
    
    for k = 1:nFiles                                                        % cycle through recordings
        fnSave = [fnPrefix num2str(dayVec(k)) '_SeriesReg_Crop.tif'];       % save file name
        bfsave(IIxyzCrop{k}, [pat fnSave], 'dimensionOrder', 'XYCZT');      % save as bioFormats TIFF    

        disp(['Saved day ' num2str(dayVec(k)) '...'])                       % user update
    end  
    
    % save ROI data
    disp('Saving crop data to parent directory...');                        % user update
    
    save([pat 'CropData.mat'], 'IIxyzCrop', 'cropRegion', 'roiCrop')        % save output variables to parent directory                                      
        
    % save output figures 
    disp('Saving figures to parent directory...');                          % user update
    
    savefig(H1, [pat 't-PMT_GFP_CropRegionOverlay_Days.fig'])
    savefig(H2, [pat 'CroppedImages_ROIoverlay_Days.fig'])
    for j = 1:nROIs
        savefig(H3{j}, [pat 'ROI' num2str(j) ' postCrop_Days.fig'])
    end
    savefig(H4, [pat 'ROI_ZaxisProfiles_Days.fig'])

end

disp('Done!');                                                              % user update

end                                                                         % end function
