%% SOME CONSTANTS
clear;
clc; %clear any previous workspaces

windowLength = 0.03; %length of window is seconds (30ms)
overlap = 0.5; %percentage of window overlap
Po = 1.2; %kg/m^3
c = 343; %m/s
Zo = Po*c; % = 411.6 Pa.s/m
r = 0.075; %1/2 dist between ears

% SET THESE TO 'y' IF WANT TO CHECK, SET AS 'n' TO NOT TEST AND SAVE TIME
qualityCheck = 'n';
localisationCheck = 'y';
diffusionCheck = 'n';

allSSIM = zeros(5,8);

%% READING THROUGH FOLDERS
% BASED OFF "https://uk.mathworks.com/matlabcentral/answers/uploaded_files/30598/recurse_subfolders.m"

% Define a starting folder.
topLevelFolder = 'C:\Users\david\webMUSHRA\configs\resources\audio\Binaural_noEQ';

allSubFolders = genpath(topLevelFolder);
% Parse into a cell array.
remain = allSubFolders;
listOfFolderNames = {}; % create folder name array
while true
	[singleSubFolder, remain] = strtok(remain, ';');
	if isempty(singleSubFolder)
		break;
	end
	listOfFolderNames = [listOfFolderNames singleSubFolder];
end
numberOfFolders = length(listOfFolderNames);  

%% PROCESSING 

for Folder=2:numberOfFolders % first folder in list is top level folder, ignore
    cd (listOfFolderNames{Folder}) % changes folder for each loop
    
    mylist = dir('*.wav'); %make list of wav files in this folder
    listlength = length(mylist); %find the length of this list 
    mylist = struct2cell(mylist); % need cell, not struct to read eachfile with audioread uu
    
    fileAziAngle = zeros(listlength,1);
    fileSSIM = zeros(listlength,1);
    filediffusion = zeros(listlength,1);  
    
    %% Read In File
    for listIndex= 1:listlength %for the amount of files in this folder
        
        fn(1,1)= mylist(1,listIndex); 
        [x,fs]=audioread(fn{1,1}); %read current file, store vals in x, fs = sampling frequency
        sigLength = length(x); %store the length of this wav file
        
        origbin(1,1) = mylist(1,1);
        [orig, fs] = audioread(origbin{1,1});
        
        leftEar = x(:,1);
        rightEar = x(:,2);
        
        %% SET-UP CALCS FOR FILE
        % general
        runtime= round((length(x)-1)/fs); % find overal time of clip in secs
        numsamp = floor(windowLength*fs); %calculate #samples per window
        update = floor(fs*overlap); %find #samples for this overlap time 

        arraySize = floor(length(x)/update); %how many windows for this clip?
        maxXval = arraySize*update; % will truncate some of the signal to round value
        xpoints = (numsamp/2:update:maxXval)/fs; % times of center of each window for plotting

        % Localisation accuracy
        %ILD = zeros(numsamp, arraySize);
        aziAngle = zeros(arraySize, 1);
        
        % timbral quality
        SSIMcoeffs = zeros(arraySize,1);
        
        %diffusion
        IACCdiff = zeros(arraySize, 1);
        
        % for processing
        pos =1; %index value, update value
        window = 1; % index for windows
        counter = 0;
        
        %% PROCESS WINDOW BY WINDOW
        while (pos+numsamp) <= length(x) % to the end of the sampled data
             y=x(pos:pos+numsamp-1, :); % y has the x sampled data of this window
             origy = orig(pos:pos+numsamp-1, :); %sampled orig data
             t =(0:length(y)-1)/fs; % finds the time each sample is at 
                    
             %for check=1:length(y) %omni value may be negative, flip for calcs
             %   if y(check,1) < 0 
             %       y(check, 2) = y(check,2)*(-1);
             %   end
             %end
             
             %% LOCALISATION ACCURACY
             if localisationCheck == 'y'
                    leftSig= y(:,1); 
                    rightSig= y(:,2); 

                    leftSigLow = lowpass(leftSig, 700, fs);
                    leftSigHigh = highpass(leftSig, 700, fs);
                    rightSigLow = lowpass(rightSig, 700, fs);
                    rightSigHigh = highpass(rightSig, 700, fs);

                    %leftFFT = fft(leftSigHigh);
                    %rightFFT = fft(rightSigHigh);
                    %ILD(:,window) = leftFFT./rightFFT;

                    [IACC1, lagleft]= xcorr(leftSigLow,rightSigLow);
                    [IACCmaxleft,I] = max(abs(IACC1));
                    lagDiff = lagleft(I);
                    timeDiffLeft = lagDiff/fs;
                    
                    [IACC2, lagright]= xcorr(rightSigLow,leftSigLow);
                    [IACCmaxright,I] = max(abs(IACC2));
                    lagDiff = lagright(I);
                    timeDiffRight = lagDiff/fs;
                    
                    if timeDiffLeft >0
                        dist = timeDiffLeft*c;
                    else 
                        dist = timeDiffRight*c;
                    end
                    if dist > 2*r
                        dist = 0;
                    end
                    aziAngle(window,1) = (asin(dist/(2*r)))*(180/pi);
             end
                    
             %% TIMBRAL QUALITY
                    
             if qualityCheck == 'y'
                 leftTest = y(:,1);
                 rightTest = y(:,2);
                 leftOrig = origy(:,1);
                 rightOrig = origy(:,2);
                 
                 spgrambw(leftTest,fs,'Jwta');
                 saveas(gcf,'testleft.png');    
                 spgrambw(rightTest,fs,'Jwta');
                 saveas(gcf,'testright.png');

                 spgrambw(leftOrig,fs,'Jwta');
                 saveas(gcf,'refleft.png');
                 spgrambw(rightOrig,fs,'Jwta');
                 saveas(gcf,'refright.png');
                 
                 close all;
                        
                 compleftref = imread('refleft.png');
                 complefttest = imread('testleft.png');
                 ssimvalleft = ssim(compleftref,complefttest);
                    
                 comprightref = imread('refright.png');
                 comprighttest = imread('testright.png');
                 ssimvalright = ssim(comprightref,comprighttest);
                 
                 delete *.png;
                    
                 SSIMcoeffs(window,1) = (ssimvalleft + ssimvalright)/2;
             end
             
             if diffusionCheck == 'y'
                 leftSig= y(:,1); % y has the x sampled data up to input ms
                 RightSig= y(:,2); 

                 % bandpass filter the signals into octaves
                 leftOct3 = bandpass(leftSig, [500, 1000], fs);
                 leftOct4 = bandpass(leftSig, [1000, 2000], fs);
                 leftOct5 = bandpass(leftSig, [2000, 4000], fs);
                    
                 rightOct3 = bandpass(RightSig, [500, 1000], fs);
                 rightOct4 = bandpass(RightSig, [1000, 2000], fs);
                 rightOct5 = bandpass(RightSig, [2000, 4000], fs);
                   
                 IACR3 = xcorr(leftOct3, rightOct3, 'coeff'); % coeff gives normalised IACR
                 IACR4 = xcorr(leftOct4, rightOct4, 'coeff'); % coeff gives normalised IACR
                 IACR5 = xcorr(leftOct5, rightOct5, 'coeff'); % coeff gives normalised IACR

                 IACR3max = max(abs(IACR3));
                 IACR4max = max(abs(IACR4));
                 IACR5max = max(abs(IACR5));

                 IACCdiff(window,1) = (IACR3max + IACR4max +IACR5max)/3;
                 
                 % according to https://eprints.soton.ac.uk/398728/1/Michael_Cousins_ICA_2016_Final.pdf
             end
            
             % update
             pos=pos+update; % This gives an overlap window of Y%
             window = window +1; %update window counter
        end % end of window         
        
        fileAziAngle(listIndex,1) = sum(aziAngle)/length(aziAngle);
        fileSSIM(listIndex, 1) = sum(SSIMcoeffs)/length(SSIMcoeffs);
        for i=1:length(IACCdiff)
            if IACCdiff(i,1) > 0
               counter = counter+1;
            end
        end       
        filediffusion(listIndex,1) = sum(IACCdiff)/counter;
    end % end of file processing
    
    allAziAngles(:,Folder-1) = fileAziAngle(:,1);
    %allSSIM(:,Folder-1) = fileSSIM(:,1);
    %allDiffusion(:, Folder-1) = filediffusion(:,1);
    % send end results for that original file + compressed files into overall results array
    % matrix reads: every col = different file
    % row1: original, row2: bitrate1, row3: bitrate2, row4: bitrate3 etc.
    
end % end of folder processing