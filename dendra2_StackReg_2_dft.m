function [I,Ir] = dendra2_StackReg_2_dft(pat,fn,tPMTchl,saveMode)
% ---------------------------------------------------------------------- % 
% Function to register confocal image stacks for dendra2 imaging.
% 
% Performs stack registration using a defined input image channel (default
% is the tranmission PMT channel).
% 
% Inputs:
% pat      - parent directory of OME TIFF image to be processed
% fn       - file name of OME TIFF image to be processed
% tPMTchl  - number of the transmission PMT image channel 
% saveMode - mode for saving registered output to parent directory (0 = off; 1 = on) 
%
% Output: 
% I        - image series
% Ir       - registered image series
%
% by Jon Witton, 24 July 2026
% ---------------------------------------------------------------------- %

% --- Format inputs ---------------------------------------------------- %
if nargin < 4
    saveMode = 0;                                                           % default save mode is off
end

if nargin < 3
    tPMTchl = 2;                                                            % default t-PMT channel (for image registration) is channel 2
end

if ~strcmp(pat(end),'\')                                                    % if parent directory name doesn't end with a backslash
    pat = [pat '\'];                                                        % add one
end


% --- Read image using Bioformats -------------------------------------- %
disp('Reading image...')                                                    % User update

rdr = BioformatsImage([pat '\' fn]);                                        % Bioformats image reader

if ~size(rdr.lut,2) == 2^16                                                 % CATCH - check the image is 16-bit
    error('Image is not 16-bit')
end

% image properties
rowPix = rdr.height;                                                        % number of row (Y) pixels
colPix = rdr.height;                                                        % number of column (X) pixels
nChls = rdr.sizeC;                                                          % number of channels
nSlices = rdr.sizeZ;                                                        % number of slices (in Z)
nFrames = rdr.sizeT;                                                        % number of frames (in time)

if ~(nFrames==1)
    error('Image contains more than one frame, i.e. timepoint')             % CATCH
else
    % generate image array 
    I = zeros(rowPix,colPix,nChls,nSlices,'uint16');                        % pre-allocate 4D image array
    
    for cc = 1:nChls
        for ss = 1:nSlices
            I(:,:,cc,ss) = getPlane(rdr,ss,cc,1);                           % read image (YXZCT format)
        end
    end
end


% --- Perform stack registration --------------------------------------- %
% --- (using Efficient Subpixel Registration) -------------------------- %

% Get registration coordinates
% (Register each frame in the stack to the preceeding one)
disp(['Getting registation coordinates using channel ' num2str(tPMTchl) '...']) % User update

usFac = 100;                                                                % upsampling factor for dft registration

regOutput = cell(nSlices,1);                                                % pre-allocate array for reistration properties

targetStack = zeros(rowPix,colPix,nSlices,'double');                        % pre-allocate array for stack of registered target frames 
targetStack(:,:,1) = double(I(:,:,tPMTchl,1));                              % first slice of target stack is first slice of raw image array 

for ss = 2:nSlices                                                          % cycle through Z-slices     

    source = double(I(:,:,tPMTchl,ss));                                     % current (source) frame (to be registered)
    target = targetStack(:,:,ss-1);                                         % target frame (preceeding frame in registered stack)

    [regOutput{ss}, sourceReg_dft] = dftregistration(fft2(target),fft2(source),usFac); % perform dft registration 
    sourceReg = abs(ifft2(sourceReg_dft));                                  % registered source image
  
    targetStack(:,:,ss) = sourceReg;                                        % add registered frame to target stack 

end


% --- Apply transformation for each Z-plane across all channels -------- %
disp('Applying transformation to all image channels...')                    % user update

% 2D array for mapping dft registration output
Nr = ifftshift((-fix(rowPix/2):ceil(rowPix/2)-1));
Nc = ifftshift((-fix(colPix/2):ceil(colPix/2)-1));
[Nc,Nr] = meshgrid(Nc,Nr);

% Pre-allocate output and apply transformation
Ir = zeros(rowPix,colPix,nChls,nSlices,'uint16');                           % pre-allocate registered output image array                         

for cc = 1:nChls                                                            % cycle through channels
    
    Ir(:,:,cc,1) = I(:,:,cc,1);                                             % add the first frame of the stack (not registered) to image array 
    
    for ss = 2:nSlices
        
        % --- transformation parameters --- %
        phase = regOutput{ss}(2);                                           % rotational phase
        deltar = regOutput{ss}(3);                                          % row (Y) shift
        deltac = regOutput{ss}(4);                                          % column (X) shift
    
        % --- apply ssth transformation to slice ss in channel cc --- % 
        source = single(I(:,:,cc,ss));                                      % current (source) frame (to be registered)

        source_dft = fft2(source);                                          % dft of source image
        source_dft_t = source_dft .* exp(1i*2*pi*(-deltar*Nr/rowPix-deltac*Nc/colPix)); % apply translation
        source_dft_tr = source_dft_t * exp(1i*phase);                       % apply rotation
        sourceReg = abs(ifft2(source_dft_tr));                              % registered source image

        Ir(:,:,cc,ss) = uint16(sourceReg);                                  % convert registered image frame from 64-bit to 16-bit and add to output array   
    
    end

end
     

% --- Save registered image series (to parent directory) --------------- %
if saveMode == 1
    disp('Saving registered image to parent directory...');                 % user update
    ind = strfind(fn,'.ome.tif');
    fnSave = [fn(1:ind-1) '_StackReg.tif'];                                 % save filename            
    bfsave(I, [pat fnSave], 'dimensionOrder', 'XYCZT');                     % save as OME TIF using bioFormats
end
   
disp('Done!');                                                              % user update

end                                                                         % end function         
