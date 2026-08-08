function dendra2_Master(fnBatch)
% dendra2_Master.m
% ---------------------------------------------------------------------- %
% Master function to run dendra2 image analysis for an image series across
% several days.
%
% The code expects:
% (1) Images to be in OME TIFF format (".ome.tif" suffix)
% (2) Image files to follow a consistent naming convension, with the
%     experiment day specifed after a space at the end of the filename, before 
%     ".ome.tif", e.g.,"Experiment 1_control_slice 2_plate 1_day 1.ome.tif'
% (3) All the image files for a single experiment (i.e., across multiple 
%     days) to be within the same folder
% (4) Data folder names to be specified in a batch file (provided as an 
%     Excel spreadsheet)
% 
% Input:
% fnBatch   : name of the Excel batch file (as a string)
% 
% By Jon Witton, 26 July 2026
% ---------------------------------------------------------------------- %

% --- Define some parameters ------------------------------------------- % 
% Image channels
tPMTchl = 2;                                                                % transmission PMT image channel
gfpChl = 1;                                                                 % GFP image channel
rfpChl = 3;                                                                 % RFP image channel

%dayVec = [1, 1.1, 3, 7, 14];                                               % vector of imaging days 

saveMode = 1;                                                               % save outputs after each step to parent directory? (0 = No, 1 = Yes)

% --- Get directory information from batch file ------------------------ %
T = readtable(fnBatch);                                                     % read batch file
nDirs = height(T);                                                          % number of directories (i.e., experiments) to process 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------- Process each directory in a for loop ------------------ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dd = 1:nDirs
    
    close all                                                               % close any open figures (e.g., from previous loop in dd)

    % --- Get directory and file names --- %
    pat = T.DirectoryNames{dd};                                             % get directory (i.e., path) name from batch file
    if ~strcmp(pat(end),'\')                                                % if path name does not end with a backslash
        pat = [pat '\'];                                                    % add one
    end

    disp('--------------------------------------------------------------')
    disp(['Processing directory: '  pat])
    disp('--------------------------------------------------------------')  % update user on which directory is being processed
    
    % --- Get file names --- %    
    filenames = dir([pat '*.ome.tif']);                                     % query the directory to get file names
    nFiles = numel(filenames);                                              % number of files in the directory to process

    % --- Make fnPrefix input  --- %
    ind = find(filenames(1).name==' ',1,'last');
    fnPrefix = filenames(1).name(1:ind);                                    % file name prefix, without the recording day

    % --- Make 'dayVec' variable if not already defined --- %
    if ~exist('dayVec','var')                                               % if dayVec variable does not exist
        dayVec = nan(1,nFiles);                                             % pre-allocate dayVec array

        for k = 1:nFiles
            ind1 = strfind(filenames(k).name,'day ')+length('day ')-1;      % file name index immediately before the day is stated
            ind2 = strfind(filenames(k).name,'.ome');                       % file name index immediately after the day is stated
            dayVec(k) = str2double(filenames(k).name(ind1+1:ind2-1));       % update day vector with recording day from file name
        end
        
        dayVec = sort(dayVec,'ascend');                                     % sort recording days in ascending order
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Register within stack (i.e., StackReg) to t-PMT channel ----- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['------------------------------------------------------\n'...
             '------ Registering within stack (i.e., StackReg) -----\n'...
             '------------------------------------------------------\n']);  % update user
    for k = 1:nFiles                                                        % cylce through image files       
        fn = [fnPrefix num2str(dayVec(k)) '.ome.tif'];                      % name of kth image file

        disp(['Processing file: ' fn]);                                     % update user on which file is being processed

        [I,Ir] = dendra2_StackReg_2_dft(pat,fn,tPMTchl,saveMode);           % Uses DFT-based subpixel registration (see Guizar-Sicairos et al. 2008. DOI: 10.1364/ol.33.000156)   
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------- Register (i.e, align in XY) each image in series --------- %
    % -------------------to GFP channel in Session 1 ------------------- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['------------------------------------------------------\n'...
             '-- Register (i.e, align in XY) each image in series --\n'...
             '------------ to GFP channel in Session 1 -------------\n'...
             '------------------------------------------------------\n']);  % update user

    [II,IIr] = dendra2_SeriesReg_dft(pat,fnPrefix,dayVec,gfpChl,saveMode);  % Uses FFT-based subpixel registration (see Guizar-Sicairos et al. 2008. DOI: 10.1364/ol.33.000156)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % - Get ROI information for GFP channel to align stacks across days -%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Allow user to define number of ROIs used for Z alignment
    prompt = 'Based on the registation, how many ROIs do you want to use for Z alignment (default = 3)?: ';
    str = input(prompt,'s');                                                % Ask user how many ROIs they want to use for Z alignment
    if isempty(str) || isnan(str2double(str))                               % If user does not provide a valid numeric input...
        nROIs = 3;                                                          % ...set default as 3 ROIs
    elseif isnumeric(str2double(str))                                       % If user does provide a valid numeric input...
        nROIs = str2double(str);                                            % ...set the input as the number of ROIs
    end

    fprintf(['------------------------------------------------------\n'...
             '-- Get ROI information to align stacks across days ---\n'...
             '------------------------------------------------------\n']);  % update user
    
    ROI = dendra2_SliceReg(pat,fnPrefix,dayVec,gfpChl,nROIs,saveMode);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Align image series in Z and crop to consistent size across days   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['------------------------------------------------------\n'...
             '- Align image series in Z across days and crop size --\n'...
             '------------------------------------------------------\n']);  % update user
    
    [IIxyzCrop,~,~] = dendra2_StackCrop(pat,fnPrefix,dayVec,saveMode,ROI);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -------- Calculate deltaF/F data for GFP and RFP channels -------- %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['------------------------------------------------------\n'...
             '-- Calculate deltaF/F data for GFP and RFP channels --\n'...
             '------------------------------------------------------\n']);  % update user
    [dFF0gfp,dFF0rfp] = dendra2_dFF0(pat,dayVec,saveMode,gfpChl,rfpChl,IIxyzCrop);


    % --- CATCH: Allow the user to terminate after each directory ------ %
    if nDirs>1 && dd<nDirs                                                  % if more than one directory is specified in the batch file, and the current directory is not the last one
        prompt2 = ['Are you happy to move onto the next data directory (i.e., experiment)?\n'...
                   'Press N to break, or any other key to continue: '];
        txt = input(prompt2,'s');                                           % ask the user if they are happy to progress to the next directory
        if strcmp(txt,'N') || strcmp(txt,'n')
            error(['Process terminated by user at directory: ' pat])        % If they input 'N' or 'n', terminate the process
        end
    end
end

end                                                                         % end function

