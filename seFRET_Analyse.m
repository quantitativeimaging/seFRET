% "seFRET_Analyse" - seFRET for Matlab
% Eric Rees, 22 Jan 2013
% Copyright 2013, University of Cambridge
%  License - GPL Version 3 or later, http://www.gnu.org/licenses/gpl.html
%  
%
% Tested with Matlab 2011b. 
% Requires image processing toolbox: for imerode()
% 
% doi:10.1098/rsif.2008.0381.focus (refer to Alan Elder paper on seFRET)
%
% Version history
% V_1-0 This is forked from "iFRET V 1-3"
%
% NOTES
% This script visualises normalised dFRET or aFRET images 
% from a user-selected set of raw seFRET intensity image data  
%
% How to use this script.
%
% First run the "seFRET_Calibrate" script to set AER, DER, Alpha, Beta
% Second, clean up the workspace if wanted - you only need those 4 numbers
%  - This script starts with a 'clearvars' command to do this.
% Third, run this script on a set of seFRET images for a real specimen
%
% Note that this script produces 8 figure by default, 
% AND automatically saves some images of the results
% - View the figure plotting commands below for information about these
% - User may (and should) delete or change Section 4 (Visualisation) 
% - as required.
%
% Useful workflows for the user
% 1. Change the size of "morphShape" to amend regions of interest
% 2. Record mean and std dev. of dFRET or aFRET in the chosen ROI
% 3. Save dFRET or aFRET figures (Masked? Colorbars? choose for yourself) 
%    -- For overlays, use the imwrite command (preserves image dimensions)
%    as well as saving a figure with the colorbar (you will need to paste
%    in this colorbar in a graphics program or something).
%

% SET INITIAL CONDITIONS
% AER = 0.7;  % Do not uncomment these inputs! Find them by calibration!
% DER = 0.5;  % These are test values suitable for the sample data ONLY
% Alpha = 0.3;
% Beta = 2.5;

clearvars -except AER DER Alpha Beta

ThresholdLow = 2000; % DmDx or AmAx signal must be sufficient to analyse
ThresholdHigh=30000; % High signals may suffer from nonlinear CCD response

flagDoNotErode = 1;  % 
flagScaleData    = 0;

% Dialog box to prompt user for inputs:
% (Comment out the following lines to use the values written above):
prompt = {'Set low threshold (darker pixels are background):',...
          'Set high threhold (to exclude bright, nonlinear pixels):',...
          'Mask erosion radius (excludes edge effects)'};
dlg_title = 'Please confirm parameters:';
num_lines = 1;
default   = {'2000', '30000', '15'};
answer    = inputdlg(prompt,dlg_title,num_lines,default);
ThresholdLow    = str2num( answer{1} );
ThresholdHigh   = str2num( answer{2} );
MaskErodeRadius = str2num( answer{3} );

% For heuristic masking, threshold and erode
morphShape = strel('disk',MaskErodeRadius); 


% 3. Process a real sample to determine its FRET image
% Find cFRET, using the values of AER, DER determined above, and images

% clear im*   % Clear cluttering image files

[FileNameAmAxSam,PathName] = uigetfile('*.tif',...
    'Open Sample File, Acceptor Emission, Acceptor Excitation'); 
imAmAx = double( imread([PathName FileNameAmAxSam]) );

[FileNameAmDxSam,PathName] = uigetfile('*.tif',...
    'Open Sample, Acceptor Emission, Donor Excitation'); 
imAmDx = double( imread([PathName FileNameAmDxSam]) );

[FileNameDmDxSam,PathName] = uigetfile('*.tif',...
    'Open Sample, Donor Emission, Donor Excitation'); 
imDmDx = double( imread([PathName FileNameDmDxSam]) );

if(flagScaleData)
    imAmAx = imresize(imAmAx,0.25);
    imAmDx = imresize(imAmDx,0.25);
    imDmDx = imresize(imDmDx,0.25);
end
    
% How do we want to threshold the sample data?
% Probably low values of imAmAx and imDmDx are desirable here
% Nonetheless, we will threshold in the same way as before
% Eroded mask for FRET estimation is within Low AND High Thresholds
mask = ( imAmAx > ThresholdLow & imAmAx < ThresholdHigh & ...
         imDmDx > ThresholdLow & imDmDx < ThresholdHigh & ...
         imAmDx > ThresholdLow & imAmDx < ThresholdHigh );
maskEroded = imerode(mask,morphShape);

% Mask for background determination is under Low Threshold ONLY
maskBGRaw = ( imAmAx < ThresholdLow & ...
           imDmDx < ThresholdLow & ...
           imAmDx < ThresholdLow );
maskBG = imerode(maskBGRaw,morphShape);

if( flagDoNotErode)
    maskEroded = mask;
    maskBG = maskBGRaw;
end

if( sum(maskEroded(:)) == 0 )
    error('Threshold excludes sample imAmAx or imDmDx entirely');
end


% Mess with images - try subtracting background
imAmAx = imAmAx - mean(imAmAx(maskBG));
imAmDx = imAmDx - mean(imAmDx(maskBG));
imDmDx = imDmDx - mean(imDmDx(maskBG));

% Calculated "cFRET", which refers to "corrected FRET"
cFRET = imAmDx - AER*imAmAx - DER*imDmDx;
cFRET(cFRET<0) = 0; % cFRET < 0 seems unphysical


% Use Alpha, Beta, and cFRET to determine dFRET and aFRET
dFRET = ( cFRET*(Alpha/DER) )./( imDmDx + cFRET*(Alpha/DER) );
aFRET = cFRET./( AER*imAmAx*Beta );


% Fix any "inf" and "nan" exceptions - set them to zero (and warn)
if( sum(isnan(aFRET(:))) || sum(isnan(aFRET(:))) || ...
    sum(isinf(dFRET(:))) || sum(isinf(dFRET(:)))  )
    warning('Some inf or nan values in dFRET or aFRET have been zeroed');
end
dFRET(isnan(dFRET)) = 0;
dFRET(isinf(dFRET)) = 0;
aFRET(isnan(aFRET)) = 0;
aFRET(isinf(aFRET)) = 0;

% Truncate dFRET and aFRET to the range [0 1] 
% - more extreme values are likely to be unphysical
if( sum(sum( (((dFRET<0) | (dFRET>1) | (aFRET<0) | (aFRET>1)).*maskEroded) )) > 0 )
    warning('Some aFRET or dFRET values outside 0-1 have been zeroed')
end
dFRET(dFRET < 0 ) = 0;
dFRET(dFRET > 1 ) = 1;
aFRET(aFRET < 0 ) = 0;
aFRET(aFRET > 1 ) = 1;


% 4. Visualisation
% Display some figures and save the image data (for overlays etc.):

colMax = 0.6; % Upper limit of colour scale for normalised FRET images

% 1) dFRET image, entire
figure(1)         
% imshow(dFRET,'border','tight')
imshow(dFRET)
colormap(jet)
caxis([0 colMax])
% set(gcf,'Position',[100,100,size(dFRET,2)+200,size(dFRET,2)]); % 
colorbar
title('dFRET image', 'FontSize',16)
set(gca,'FontSize',16,'fontweight','bold');

% 2) dFRET image, masked (non-good signal area shown as white)
%    This shows the region used to determine mean values
%    Display this figure with a colorbar, but save the image without one
figure(2)        
    % Upper limit of color axis range
myIm = dFRET;
myIm(not(maskEroded)) = 0;
imshow(myIm)
% set(gcf,'Position',[100,100,size(dFRET,2)+200,size(dFRET,2)]); % 
colormap(jet(256))
caxis([0 colMax])
cmap = colormap;
cmap(1,:) = [1,1,1]; % Set uncertain (zero dFRET) region to white
colormap(cmap);
colorbar
title('dFRET image, masked to show areas of FRET plotted on histogram'...
    , 'FontSize',14);
set(gca,'FontSize',16,'fontweight','bold');

% % Save png with white mask and colorbar
thisIm = getframe(gcf);
thisIm = thisIm.cdata;
% UNCOMMENT THE FOLLOWING TO SAVE THE GRAPH:
% imwrite(thisIm,'dFRET_im_masked_with_colorbar.png')

% ALTERNATIVE: uncomment the next line to save without colorbar
% UNCOMMENT THE FOLLOWING TO SAVE GRAPH:
% imwrite(myIm*256*(1/colMax),cmap,'dFRET_im_masked.png','png') 


% 3. Histogram of dFRET values (in the maskEroded region of "good" signal)
figure(3)
hist(dFRET(maskEroded), (0:0.01:1) )
xlim([0 1])
xlabel('dFRET', 'FontSize',16)
ylabel('Number of pixels', 'FontSize',16)
title('dFRET histogram', 'FontSize',16)
set(gca,'FontSize',16,'fontweight','bold');

% 4. aFRET image, entire
figure(4)
imshow(aFRET)
colormap(jet)
caxis([0 colMax])
colorbar
title('aFRET image', 'FontSize',16)
set(gca,'FontSize',16,'fontweight','bold');

% 5. aFRET image, masked (non-good signal area shown as white)
figure(5)        
myIm = aFRET;
myIm(not(maskEroded)) = 0;
imshow(myIm)
colormap(jet(256))
caxis([0 colMax])
cmap = colormap;
cmap(1,:) = [1,1,1]; % Set uncertain (zero dFRET) region to white
colormap(cmap);
colorbar
title('aFRET image, masked to show areas of FRET plotted on histogram'...
    , 'FontSize',14);
set(gca,'FontSize',16,'fontweight','bold');

% % Save png with white mask and a colorbar
thisIm = getframe(gcf);
thisIm = thisIm.cdata;
% UNCOMMENT THE FOLLOWING TO SAVE GRAPH:
% imwrite(thisIm,'aFRET_im_masked_with_colorbar.png')

% Save image without border or colorbar
% imwrite(myIm*256*(1/colMax),cmap,'aFRET_im_masked.png','png') 

% 6. Histogram of dFRET values (in the maskEroded region of "good" signal)
figure(6)
hist(aFRET(maskEroded), (0:0.01:1) )
xlim([0 1])
xlabel('aFRET', 'FontSize',16);
ylabel('Number of pixels', 'FontSize',16);
title('aFRET histogram', 'FontSize',16)
set(gca,'FontSize',16,'fontweight','bold');

% 7. Print mean and standard deviation of dFRET, aFRET signals
% Note that Standard deviation of pixel values may be a combination of:
% (a) Real variation across the sample area, AND
% (b) Random measurement errors
% CONSIDER INTENSITY-WEIGHTING dFRET and aFRET averages, for final results
calculated_dFRET_mean = mean(dFRET(maskEroded))
calculated_dFRET_std  =  std(dFRET(maskEroded))

calculated_aFRET_mean = mean(aFRET(maskEroded))
calculated_aFRET_std  =  std(aFRET(maskEroded))


% The following Figures 7 and 8 are not essential - commented out.
% % 8. Print imDmDx and imAmAx (after background subtraction)
% % Note that background subtraction has been applied - this may not longer
% % show that the raw data is in the linear 2000 - 30000 region!
% figure(7)
% imshow(imDmDx)
% colormap(jet)
% caxis([2000 30000]) 
% set(gcf,'Position') %,[100,100,size(dFRET,2)+200,size(dFRET,2)]); % 
% colorbar
% set(gca,'FontSize',16,'fontweight','bold');
% 
% figure(8)
% imshow(imAmAx)
% colormap(jet)
% caxis([2000 30000]) 
% set(gcf,'Position') % ,[100,100,size(dFRET,2)+200,size(dFRET,2)]); % 
% colorbar
% set(gca,'FontSize',16,'fontweight','bold');

% We have now determined the normalised FRET images, 
% within a mask which specifies only regions with "good" signal values
% And we have printed the mean aFRET, dFRET values to the console


% To Do:
% DONE 1. Background Subtraction - discuss 
% DONE 2. inf / nan errors
% 2b. Halting on error messages (done, changed to warnings)
% DONE 3. masking - for AER, DER selection - display - and for output
% DONE 4. Setting Alpha and Beta
% 5. Find best width for erosion filter (15 px is OK for example data)