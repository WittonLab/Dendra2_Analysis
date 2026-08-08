function [II,IIr] = dendra2_SeriesReg_dft(pat,fnPrefix,dayVec,regChl,saveMode)
% ---------------------------------------------------------------------- % 
% Function to align in XY confocal image stacks for dendra2 imaging across
% days.
% 
% Inputs:
% pat      - parent directory of registered .tif image
% fn       - file name prefix (omitting the day number and '_StackReg.tif'
% regChl   - channel used for image registration (default is chl 1 - GFP)
% dayVec   - (1,m) numerical vector of image days (i.e., [1, 1.1, 3, 7, 14])
% saveMode - mode for saving registered output to parent directory (0 = off; 1 = on) 
%
% Outputs:
% II       - (1,m) cell array containing 4D image data arrays of format YXCZ
% IIr      - (1,m) cell array containing 4D image data arrays (YXCZ format) registered to Day 1 regChl 
%
% by Jon Witton, 24 July 2026
% ---------------------------------------------------------------------- %

% --- Format inputs ---------------------------------------------------- %
if nargin < 5
    saveMode = 0;                                                           % default save mode is off
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

    rdr = BioformatsImage([pat fnPrefix num2str(dayVec(k)) '_StackReg.tif']);
    
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

IIm = zeros(rowPix,colPix,nFiles,'uint16');                                 % pre-allocate projection image array

for k = 1:nFiles

    I = squeeze(II{k}(:,:,regChl,:));                                       % image series for registration channel on kth day 
    IIm(:,:,k) = uint16(mean(single(I),3));                                 % mean projection (preserving 16-bit depth)

end


% --- Register projection images (to first recording) ------------------ %
% --- (using Efficient Subpixel Registration) -------------------------- %

% Perform registration
disp('Registering projection images...')                                    % User update

usFac = 100;                                                                % upsampling factor for dft registration

IImr = zeros(rowPix,colPix,nFiles,'uint16');                                % pre-allocate registered image array
IImr(:,:,1) = IIm(:,:,1);                                                   % allocate day 1 projection image (not registered) to registered image array  

regOutput = cell(1,nFiles);                                                 % pre-allocate array for registration properties

target = double(IIm(:,:,1));                                                % target image (day 1 project image baseline)

for k = 2:nFiles                                                            % cycle through each recording

    source = double(IIm(:,:,k));                                            % frame to be registered

    [regOutput{k}, sourceReg_dft] = dftregistration(fft2(target),fft2(source),usFac); % perform dft registration 
    sourceReg = abs(ifft2(sourceReg_dft));                                  % registered source image

    IImr(:,:,k) = uint16(sourceReg);                                        % add registeered frame to registered frames array  

end


% --- Display reference image ------------------------------------------ %
disp('Displaying registered projection images...')                          % user update

day1Target = IImr(:,:,1);                                                   % registration target (Day 1 image)

figure('Name','Registered','WindowState', 'maximized')

for k = 1:nFiles
    
    % display registered images
    subplot(2,nFiles,k) 
    imagesc(IImr(:,:,k))                                                    % display
    axis image
    title(['Day ' num2str(dayVec(k))])
    
    % display registered imaged overlaid on day 1 target
    subplot(2,nFiles,nFiles+k)
    imshowpair(day1Target,IImr(:,:,k),'scaling','independent')              % registered image overlaid on Day 1 target image
    axis image
    axis on
    title(['Day ' num2str(dayVec(k)) ' vs. Day 1'])

end


% --- Check user satisfaction with registration ------------------------ %
prompt = 'Are you happy with this registration? \n Press any key to continue / press N to cancel: ';
txt = input(prompt,'s');
if strcmp(txt,'N') || strcmp(txt,'n')
    error('User break - registration has not worked')
end

close(gcf)                                                                  % close figure


% --- Apply transformations per each day to all images in series ------- %
disp('Applying transformation per day to all images in series...')          % User update

% 2D array for mapping dft registration output
Nr = ifftshift((-fix(rowPix/2):ceil(rowPix/2)-1));
Nc = ifftshift((-fix(colPix/2):ceil(colPix/2)-1));
[Nc,Nr] = meshgrid(Nc,Nr);

% Pre-allocate output and apply transformation
IIr = cell(1,nFiles);                                                       % pre-allocate array for all registered image files

IIr{1} = II{1};                                                             % add Day 1 (not registered) to output array                                                 

for k = 2:nFiles                                                            % cycle through recordings
    
    IIr{k} = zeros(rowPix,colPix,nChls,nSlices,'uint16');                   % pre-allocate array for kth day registered image series     

    IItmp = II{k};                                                          % image series to be registered

    % --- transformation parameters --- %
    phase = regOutput{k}(2);                                                % rotational phase
    deltar = regOutput{k}(3);                                               % row (Y) shift
    deltac = regOutput{k}(4);                                               % column (X) shift

    % --- cycle through images and apply transformation --- %
    for cc = 1:nChls
        for ss = 1:nSlices
            
            source = double(IItmp(:,:,cc,ss));                              % source image to be registered

            source_dft = fft2(source);                                      % dft of source image
            source_dft_t = source_dft .* exp(1i*2*pi*(-deltar*Nr/rowPix-deltac*Nc/colPix)); % apply translation
            source_dft_tr = source_dft_t * exp(1i*phase);                   % apply rotation
            sourceReg = abs(ifft2(source_dft_tr));                          % registered source image

            IIr{k}(:,:,cc,ss) = uint16(sourceReg);                          % convert registered image frame from 64-bit to 16-bit and add to output array   

        end
    end

end


% --- Save registered image series (to parent directory) --------------- %
if saveMode == 1
    disp('Saving registered images to parent directory...');                % user update
    
    for k = 1:nFiles                                                        % cycle through recordings
        fnSave = [fnPrefix num2str(dayVec(k)) '_SeriesReg.tif'];            % save file name
        bfsave(IIr{k}, [pat fnSave], 'dimensionOrder', 'XYCZT');            % save as bioFormats TIFF    
        
        disp(['Saved day ' num2str(dayVec(k)) '...'])                       % user update
    end   

end
   
disp('Done!');                                                              % user update

end                                                                         % end function


