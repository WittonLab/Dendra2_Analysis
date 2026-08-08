function [dFF0gfp,dFF0rfp] = dendra2_dFF0(pat,dayVec,saveMode,gfpChl,rfpChl,IIxyzCrop)
% ---------------------------------------------------------------------- % 
% Function to calculate deltaF/F0 information for GFP and RFP image 
% channels for dendra2 image series across days.
% 
% Inputs:
% pat       - parent directory of cropped image series (output from dendra2_StackCrop.m)
% dayVec    - (1,m) numerical vector of image days (i.e., [1, 1.1, 3, 7, 14])
% saveMode  - Mode for saving output dFF0 data (0 = off; 1 = on) 
% GFPchl    - number of GFP image channel (default = 1)
% RFPchl    - number of RFP image channel (default = 3)
% IIxyzCrop - Output from dendra2_StackCrop.m (OPTIONAL: read from pat if not input)
%             Is a (1,k) cell array containing cropped images, across days (in k)
%             Each image array (within a cell) is a 4D array of XYCZT format
% 
% Outputs:
% dFF0gfp - (1,k) vector containing deltaF/F0 for GFP channel, across days (in k)
% dFF0rfp - (1,k) vector containing deltaF/F0 for GFP channel, across days (in k)
%
% by Jon Witton, 24 July 2026
% ---------------------------------------------------------------------- %

if nargin < 6
    tmp = load([pat 'CropData.mat'],'IIxyzCrop');                           % read IIxyzCrop from parent directory if not input
    IIxyzCrop = tmp.IIxyzCrop;
end

if nargin < 5, rfpChl = 3; end                                              % default RFP image channel

if nargin < 4, gfpChl = 1; end                                              % default GFP image channel

if nargin < 3, saveMode = 0; end                                            % default save mode is off

if nargin < 2, dayVec = [1, 1.1, 3, 7, 14]; end                             % default imaging days


% --- Define baselines ------------------------------------------------- %
gfpBaseline = 1;                                                            % recording/file/day number of baseline for GFP channel
rfpBaseline = 1.1;                                                          % recording/file/day number of baseline for RFP channel


% --- Get fluorescence intensity data ---------------------------------- %
nFiles = numel(IIxyzCrop);                                                  % number of files to process

IgfpSig = nan(1,nFiles);                                                    % preallocate array for GFP signal (across days)                                     
IrfpSig = nan(1,nFiles);                                                    % preallocate array for RFP signal (across days)    

for k = 1:nFiles                                                                                     

    I = IIxyzCrop{k};                                                       % image series (all channels) for kth file 
    
    Igfp = squeeze(I(:,:,gfpChl,:));                                        % gfp channel for kth file
    IgfpSig(k) = mean(Igfp(:));                                             % mean pixel intensity in GFP channel (across all Y-rows, X-columns, and Z-slices)

    Irfp = squeeze(I(:,:,rfpChl,:));                                        % gfp channel for kth file
    IrfpSig(k) = mean(Irfp(:));                                             % mean pixel intensity in RFP channel (across all Y-rows, X-columns, and Z-slices)

end


% --- Calculate deltaF/F0 ---------------------------------------------- %
disp('Calculating deltaF/F0...');                                           % user update

% GFP
gfpBaseInd = find(dayVec==gfpBaseline,1,'first');                           % index of recording for GFP channel baseline
F0gfp = IgfpSig(gfpBaseInd);                                                % F0 (baseline) for GFP channel
dFF0gfp = (IgfpSig - F0gfp) ./ F0gfp;                                       % calculate deltaF/F0 for GFP channel

% RFP
rfpBaseInd = find(dayVec==rfpBaseline,1,'first');                           % index of recording for RFP channel baseline
F0rfp = IrfpSig(rfpBaseInd);                                                % F0 (baseline) for RFP channel
dFF0rfp = (IrfpSig - F0rfp) ./ F0rfp;                                       % calculate deltaF/F0 for RFP channel


% --- Plot figures ----------------------------------------------------- %
% --- Raw signal figure --- %
H1 = figure('Name','Dendra2 Mean Raw Fluorescence');

% GFP
yyaxis('left');
plot(dayVec(gfpBaseInd:end),IgfpSig(gfpBaseInd:end),'-og')                  % display data                                                
ylabel('GFP Mean Raw Fluorescence (AU)')

%RFP
yyaxis('right');
plot(dayVec(rfpBaseInd:end),IrfpSig(rfpBaseInd:end),'-or')                  % display data                                                
ylabel('RFP Mean Raw Fluorescence (AU)')

xlabel('Day')


% --- dFF0 figure --- %
H2 = figure('Name','Dendra2 DeltaF/0');

% GFP
yyaxis('left');
plot(dayVec(gfpBaseInd:end),dFF0gfp(gfpBaseInd:end),'-og')                  % display data                                                
hold on
plot([min(xlim) max(xlim)],[0 0],'--g')                                     % plot referline line through y=0
ylabel('GFP \DeltaF/F0')

%RFP
yyaxis('right');
plot(dayVec(rfpBaseInd:end),dFF0rfp(rfpBaseInd:end),'-or')                  % display data                                                
hold on
plot([min(xlim) max(xlim)],[0 0],'--r')                                     % plot referline line through y=0
ylabel('RFP \DeltaF/F0')

xlabel('Day')


% --- Save output variables and figures (to parent directory) ---------- %
if saveMode == 1
    disp('Saving fluorescence signal data to parent directory...');         % user update

    % save data as .mat file
    save([pat 'signalData_Days.mat'], 'IgfpSig', 'dFF0gfp', 'IrfpSig', 'dFF0rfp'); % save deltaF/F variables  

    % save data as Excel file
    dataTab = table(dayVec', IgfpSig', dFF0gfp', IrfpSig', dFF0rfp', 'VariableNames', {'RecordingDay','GFP_RawSignal','GFP_dFF0','RFP_RawSignal','RFP_dFF0'}); % format data table
    writetable(dataTab,[pat 'signalData_Days.xlsx'])                        % write file
        
    % save figures     
    disp('Saving figures to parent directory...');                          % user update
    
    savefig(H1, [pat 'MeanRawFluorescence_Days.fig'])
    savefig(H2, [pat 'dFF0_Days.fig'])
end

disp('Done!');                                                              % user update

end                                                                         % end function         
