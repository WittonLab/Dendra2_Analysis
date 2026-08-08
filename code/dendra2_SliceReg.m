function ROI = dendra2_SliceReg(pat,fnPrefix,dayVec,regChl,nROIs,saveMode)
% ---------------------------------------------------------------------- % 
% Function to align Z slices across confocal image stacks for dendra2 
% imaging across days.
% 
% Inputs:
% pat      - parent directory of registered .tif image
% fn       - file name prefix (omitting the day number and '_StackReg.tif'
% dayVec   - (1,m) numerical vector of image days (i.e., [1, 1.1, 3, 7, 14])
% regChl   - channel used for slice registration (default is chl 1 - GFP)
% nROIs    - number of manually drawn ROIs used to align Z-slices (default = 3)
% saveMode - mode for saving registered output to parent directory (0 = off; 1 = on) 
%
% Outputs:
% ROI      - structure containing region of interest (ROI) variables, namely:
%            roiCoords         - ROI pixel coordinates                                                        
%            IIroi             - ROI image stacks
%            MnProjection      - ROI mean projection images
%            Zprofiles         - ROI Z-axis profiles
%            NormZprofiles     - ROI normalised Z-axis profiles
%            ZprofilesPk       - Z-slice correspondeing to peak ROI intensity 
%            ZprofilesPkDiff   - difference in peak Z-intensities between first file and subsequent files 
%            ZprofilesPkDiffMn - mean ROI difference in peak Z-intensities between first file and subsequent files
%
% by Jon Witton, 24 July 2026
% ---------------------------------------------------------------------- %

% --- Format inputs ---------------------------------------------------- %
if nargin < 6
    saveMode = 0;                                                           % default save mode is off
end

if nargin < 5
    nROIs = 3;                                                              % default number of ROIs used to align Z-slices
end

if nargin < 4
    regChl = 1;                                                             % default registration channel
end

if nargin < 3
    dayVec = [1, 1.1, 3, 7, 14];                                            % default imaging days
end

nFiles = length(dayVec);                                                    % number of files to algin

if ~strcmp(pat(end),'\')                                                    % if parent directory name doesn't end with a backslash
    pat = [pat '\'];                                                        % add one
end


% --- Read each image file into a common cell array using Bioformats --- %
II = cell(1,nFiles);                                                        % array for all image files

for k = 1:nFiles
    disp(['Reading day ' num2str(dayVec(k)) '...'])                         % user update

    rdr = BioformatsImage([pat fnPrefix num2str(dayVec(k)) '_SeriesReg.tif']);
    
    % image properties
    rowPix = rdr.height;                                                    % number of row (Y) pixels
    colPix = rdr.height;                                                    % number of column (X) pixels
    nChls = rdr.sizeC;                                                      % number of channels
    nSlices = rdr.sizeZ;                                                    % number of slices (in Z)
    nFrames = rdr.sizeT;                                                    % number of frames (in time)

    if ~(nFrames==1)
        error('Image contains more than one frame, i.e. timepoint')         % CATCH
    else 
        % generate image array
        I = zeros(rowPix,colPix,nChls,nSlices,'uint16');                    % pre-allocate 4D image array
    
        for cc = 1:nChls
            for ss = 1:nSlices
                I(:,:,cc,ss) = getPlane(rdr,ss,cc,1);                       % read image (YXZCT format)
            end
        end
        
        % Add image 4D array to parrent cell array
        II{k} = I;   
    end
end


% --- Make Z-projection of registration channel ------------------------ %
disp('Calculating projection images...')                                    % user update

IImn = zeros(rowPix,colPix,nFiles,'uint16');                                % pre-allocate mean projection image array
IImx = zeros(rowPix,colPix,nFiles,'uint16');                                % pre-allocate mean projection image array

for k = 1:nFiles

    I = squeeze(II{k}(:,:,regChl,:));                                       % image series for registration channel on kth day 
    
    IImn(:,:,k) = uint16(mean(single(I),3));                                % mean projection (preserving 16-bit depth)
    IImx(:,:,k) = uint16(max(single(I),[],3));                              % max projection (preserving 16-bit depth)
end


% --- Make reference images and draw ROIs ------------------------------ %
txt = [];                                                                   % make text array for user confirmation of ROI selection

while ~strcmp(txt,'Y') && ~strcmp(txt,'y')                                  % run the following code until user sets txt to 'Y' or 'y' (for "yes")

    close all                                                               % close any open figures    
     
    % Display projection images across days
    disp('Displaying projection images...')                                 % user update
    
    H1 = figure('Name','Registered','WindowState','maximized');
    
    for k = 1:nFiles
        
        % mean projection
        subplot(2,nFiles,k) 
        imagesc(IImn(:,:,k))                                                % display mean projection
        axis image
        title(['Day ' num2str(dayVec(k)) ' - Mean Projection'])
        hold on
        
        % max projection
        subplot(2,nFiles,nFiles+k)
        imagesc(IImn(:,:,k))                                                % display max projection
        axis image
        axis on
        title(['Day ' num2str(dayVec(k)) ' - Max Projection'])
        hold on
    end
    
    
    % Display projection images for Day 1 and draw ROIs on them
    H2 = figure('Name','DRAW ROIs ON PROJECTION IMAGES','WindowState','maximized');
    
    ax1 = subplot(1,2,1);
    imagesc(IImn(:,:,1));                                                   % display Day 1 mean projection
    axis image
    colorbar
    title(['Day ' num2str(dayVec(1)) ' - Mean Projection'])
    hold on
    
    ax2 = subplot(1,2,2);
    imagesc(IImx(:,:,1));                                                   % display Day 1 max projection
    axis image
    colorbar
    title(['Day ' num2str(dayVec(1)) ' - Max Projection'])
    hold on
    
    linkaxes([ax1 ax2],'xy')
    
    C = cell(nROIs,1);                                                      % pre-allocate user input coordinates array
    lineX = cell(nROIs,1);                                                  % pre-allocate array for X-axis ROI line coordinates
    lineY = cell(nROIs,1);                                                  % pre-allocate array for Y-axis ROI line coordinates
    
    for p = 1:nROIs
        
        C{p} = round(ginput(2));                                            % ROI coordinates - select top left and bottom right edges
        
        lineX{p} = [C{p}(1,1), C{p}(2,1), C{p}(2,1) C{p}(1,1), C{p}(1,1)];  % ROI X-axis coordinates
        lineY{p} = [C{p}(1,2), C{p}(1,2), C{p}(2,2) C{p}(2,2), C{p}(1,2)];  % ROI Y-axis coordinates
        
        plot(ax1,lineX{p},lineY{p},'color','m','linewidth',2)               % display ROI on mean projection
        plot(ax2,lineX{p},lineY{p},'color','m','linewidth',2)               % display ROI on max projection
    end
    
    % Update reference images with ROIs
    figure(H1)                                                              % make H1 reference image active
    
    for k = 1:nFiles*2
        
        % display registered images
        subplot(2,nFiles,k) 
        for p = 1:nROIs
            plot(lineX{p},lineY{p},'color','m','linewidth',1)               % overlay ROI
        end
    end
    
    
    % --- Check user satisfaction with ROI selection ----------------------- %
    prompt = 'Are you happy with your ROI selection? \n Press Y to continue / press any other key to re-do: ';
    txt = input(prompt,'s');
    
end                                                                         % end while loop  


% --- Display ROIs and calculate Z-axis profiles ----------------------- %
H3 = figure('Name','ROIs','WindowState','maximized');
H4 = figure('Name','ROI Z-axis profiles','WindowState','maximized');

IIroi = cell(nROIs,nFiles);                                                 % pre-allocate ROI cell array [(p,k) format, where ROI number (p) is in the rows and file number (k) is in the columns] 
IIroiMn = cell(nROIs,nFiles);                                               % pre-allocate ROI projection image cell array [(p,k) format, where ROI number (p) is in the rows and file number (k) is in the columns]
IIroi_Zprofile = cell(nROIs,nFiles);                                        % pre-allocate ROI Z-axis profiles cell array [(p,k) format, where ROI number (p) is in the rows and file number (k) is in the columns] 
IIroi_ZprofileNorm = cell(nROIs,nFiles);                                    % pre-allocate normalised ROI Z-axis profiles cell array [(p,k) format, where ROI number (p) is in the rows and file number (k) is in the columns] 
IIroi_ZprofilePk = nan(nROIs,nFiles);                                       % pre-allocate array for Z-profile slices containing peak intensity [(p,k) format, where ROI number (p) is in the rows and file number (k) is in the columns]

for k = 1:nFiles

    I = squeeze(II{k}(:,:,regChl,:));                                       % image series for registration channel on kth day 

    for p = 1:nROIs
        
        % --- get ROI --- %
        IIroi{p,k} = I(C{p}(1,2):C{p}(2,2), C{p}(1,1):C{p}(2,1), :);        % crop of pth ROI 
        
        IIroiMn{p,k} = uint16(mean(single(IIroi{p,k}),3));                  % mean projection of pth ROI 

        % --- display ROI image --- %
        figure(H3)                                                          % set H3 as the current figure 
        subplotPos = (p-1)*nFiles+k;                                        % subplot position on figure                         
        subplot(nROIs,nFiles,subplotPos);
        imagesc(IIroiMn{p,k})                                               % plot projection image of pth ROI
        axis image
        title(['ROI ' num2str(p) ' - Day ' num2str(dayVec(k))])

        % --- calculate Z-axis profile --- %
        profile = nan(nSlices,1);                                           % pre-allocate temporary array for Z-axis profile     
        
        for j = 1:nSlices      
            tmpIroi = IIroi{p,k}(:,:,j);                                    % jth frame of current (pth) ROI for current (kth) image file
            profile(j) = mean(double(tmpIroi(:)));                          % find mean intensity of current frame and add to Z axis profile    
        end
        
        smProfile = smooth(profile,3);                                      % smooth profile with 3 point average
        
        minSubProfile = smProfile-min(smProfile);                           % subtract minimum from profile
        pkNormProfile = minSubProfile/max(minSubProfile);                   % divide profile by maximum to normalise intensity between 0 and 1
    
        IIroi_Zprofile{p,k} = smProfile;                                    % find mean intensity of current frame and add to Z axis profile
        IIroi_ZprofileNorm{p,k} = pkNormProfile;
    
        [~,pkSliceInd] = max(pkNormProfile);                                % index of Z-slice corresponding to peak ROI intensity
        IIroi_ZprofilePk(p,k) = pkSliceInd;

        % --- display Z-axis profile --- %
        figure(H4)                                                          % set H4 as the current figure 
        subplot(1,nROIs,p), hold on;
        plot(IIroi_ZprofileNorm{p,k},1:nSlices)                             % plot Z-axis profile
        set(gca,'ydir','reverse')                                           % set Y-axis direction as 'reverse', so slice number increases from top to bottom of the figure
        set(gca,'ylim',[1 nSlices])                                         % set Y-axis limits
        title(['ROI ' num2str(p)])
        xlabel('Normalised intensity')
        ylabel('Z-slice number')
        legend('location','northeast')    
    end
end


% --- Find difference in location of peak ROI intensity across images -- %
IIroi_ZprofilePkDiff = IIroi_ZprofilePk - IIroi_ZprofilePk(:,1);            % all files minus day 1
IIroi_ZprofilePkDiffMn = round(mean(IIroi_ZprofilePkDiff,1));               % average across all ROIs


% --- Make output structure -------------------------------------------- %
ROI.roiCoords = C;                                                          % ROI pixel coordinates                                                        
ROI.oimg = IIroi;                                                           % ROI image stacks
ROI.MnProjection = IIroiMn;                                                 % ROI mean projection images
ROI.Zprofiles = IIroi_Zprofile;                                             % ROI Z-axis profiles
ROI.NormZprofiles = IIroi_ZprofileNorm;                                     % ROI normalised Z-axis profiles
ROI.ZprofilesPk = IIroi_ZprofilePk;                                         % Z-slice correspondeing to peak ROI intensity 
ROI.ZprofilesPkDiff = IIroi_ZprofilePkDiff;                                 % difference in peak Z-intensities between first file and subsequent files 
ROI.ZprofilesPkDiffMn = IIroi_ZprofilePkDiffMn;                             % mean ROI difference in peak Z-intensities between first file and subsequent files


% --- Save output structrue (to parent directory) ---------------------- %
if saveMode == 1
    disp('Saving output structure...');                                     % user update
    
    % save output structure to parent directory 
    save([pat 'ROIdata.mat'], 'ROI')                                           
    
    % save figures to parent directory 
    savefig(H1,[pat 'ROIsOnProjection_Days.fig'])
    savefig(H2,[pat 'ROIsOnProjection_Day1_Drawn.fig'])
    savefig(H3,[pat 'ROIcrops_Days.fig'])
    savefig(H4,[pat 'ROI_ZaxisProfiles_Days.fig'])
end

disp('Done!');                                                              % user update

end                                                                         % end function

