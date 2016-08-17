% "seFRET_Calibrate" - seFRET example in Matlab
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
% 1.0 Forked from iFRET 1.3 (see comments and documentation in that file)
% 
% NOTES:
% This script will determine mean calibration values for seFRET, by
% processing user-selected reference image files. 
%
% Abbreviations:
% AOS    = Acceptor only sample
% DOS    = Donor only sample
% Linker = CFP-YFP linker sample, with known dFRET and aFRET
% AmAx   = Acceptor Emission, Acceptor Excitation spectral channel
% AmDx   = Acceptor Emission, Donor Excitation spectral channel
% DmDx   = Donor Emission, Donor Excitation spectral channel
% Im     = image (assumed 16 bit)
%
% ThresholdLow: pixels within a region darker than ThresholdLow are 
%               assumed to be background. Used for d.c. background 
%               estimation and subtraction
% ThresholdHigh: pixels within a region that is brighter than ThresholdLow
%                AND darker than ThresholdHigh are identified as signals.
%                The reason for ThresholdHigh is that some detectors may 
%                become nonlinear at high brightness, and such nonlinear
%                regions are unsuitable for this linear seFRET analysis.
%                
%     ThresholdLow and ThresholdHigh are set at the start of this script.
%     They may (and should) be optimised by the user for their own system.
% 
% KnownDFRET = the dFRET value of the Linker construct sample, 
%              which must be known (e.g. from TCSPC / lifetime measurement)
%              in order for alpha and beta to be fitted in this script.


% 0. SET INITIAL CONDITIONS
%    Users should check these parameters are OK.
%    In particular:
%    (a) Set KnownDFRET from lifetime calibration data, and
%    (b) Set reasonable Threshold levels for your images
%
clear all

ThresholdLow = 2000; % DmDx or AmAx signal must be sufficient to analyse
ThresholdHigh=30000; % High signals may suffer from nonlinear CCD response

KnownDFRET   = 0.30; % Specify the known dFRET of the linker images
                     % E.g. Use 0.30 for 17AA sample

MaskErodeRadius = 15;% Erode edges of masked regions by this much
        
flagDoNotErode = 1;

% Dialog box to prompt user for inputs:
% (Comment out the following lines to use the values written above):
prompt = {'Set low threshold (darker pixels are background):',...
          'Set high threhold (to exclude bright, nonlinear pixels):',...
          'dFRET value of Linker Construct:',...
          'Mask erosion radius (excludes edge effects)'};
dlg_title = 'Please confirm parameters:';
num_lines = 1;
default   = {'2000', '30000', '0.30', '15'};
answer    = inputdlg(prompt,dlg_title,num_lines,default);
ThresholdLow    = str2num( answer{1} );
ThresholdHigh   = str2num( answer{2} );
KnownDFRET      = str2num( answer{3} );
MaskErodeRadius = str2num( answer{4} );

% For heuristic masking, threshold and erode
morphShape = strel('disk',MaskErodeRadius); 


% 1a. AER
% Find AER from an acceptor-only sample (e.g. m-CHERRY only)
% Threshold the imAmAx signal
% Erode the threshold mask to guess a clean region for determining AER

% Locate folder containing AOS image set(s) to process
% - assuming the Leica "ch00, ch01..." etc name pattern. 
[FileNameAmAxAOS,PathName] = uigetfile('*.tif',...
    'Open Acceptor Only Sample, Acceptor Emission, Acceptor Excitation');      
imAmAxAOS = double( imread([PathName FileNameAmAxAOS]) );

[FileNameAmDxAOS,PathName] = uigetfile('*.tif',...
    'Open Acceptor Only Sample, Acceptor Emission, Donor Excitation'); 
imAmDxAOS = double( imread([PathName FileNameAmDxAOS]) );

% Eroded mask for AER determination is within Low AND High Thresholds
maskAOS = (imAmAxAOS > ThresholdLow & imAmAxAOS < ThresholdHigh & ...
           imAmDxAOS > ThresholdLow & imAmDxAOS < ThresholdHigh );
maskAOSeroded = imerode(maskAOS,morphShape);

% Mask for background determination starts under the Low Threshold ONLY
maskAOSbgRaw = (imAmAxAOS < ThresholdLow & imAmDxAOS < ThresholdLow);
maskAOSbg = imerode(maskAOSbgRaw,morphShape);

if (flagDoNotErode)
    maskAOSeroded = maskAOS;
    maskAOSbg = maskAOSbgRaw;
end

if(sum(maskAOSeroded(:)) == 0 )
    error('No data within cautious AOS mask')
  else if(sum(maskAOSbg(:)) == 0)
    error('No region could be found to define an AOS background')
       end
end

% Try to subtract the image background level
imAmAxAOS = imAmAxAOS - mean(imAmAxAOS(maskAOSbg));
imAmDxAOS = imAmDxAOS - mean(imAmDxAOS(maskAOSbg));

% Create an image matrix of the acceptor emission ratio
imRatioAOS = zeros(size(imAmAxAOS));
imRatioAOS(maskAOS) = imAmDxAOS(maskAOS)./imAmAxAOS(maskAOS);

% Find the mean AER (e.g. 0.5)
AER = mean(imRatioAOS(maskAOSeroded));
% Consider giving a warning if AER > 1 (unlikely), or AER < 0 (unphysical).

AER % Echo AER to MATLAB console


%%
% 1b. DER
% Find DER from a donor-only sample (e.g. GFP only)
% Threshold the imDmDx signal
% Erode the threshold mask to guess a clean region for determining DER

[FileNameAmDxDOS,PathName] = uigetfile('*.tif',...
    'Open Donor Only Sample, Acceptor Emission, Donor Excitation');  
imAmDxDOS = double( imread([PathName FileNameAmDxDOS]) );

[FileNameDmDxDOS,PathName] = uigetfile('*.tif',...
    'Open Donor Only Sample, Donor Emission, Donor Excitation');  
imDmDxDOS = double( imread([PathName FileNameDmDxDOS]) );

% Eroded mask for DER determination is within Low AND High Thresholds
maskDOS = ( imDmDxDOS > ThresholdLow & imDmDxDOS < ThresholdHigh & ...
            imAmDxDOS > ThresholdLow & imAmDxDOS < ThresholdHigh );
maskDOSeroded = imerode(maskDOS,morphShape);

% Mask for background determination starts under the Low Threshold ONLY
maskDOSbgRaw = ( imDmDxDOS < ThresholdLow & imAmDxDOS < ThresholdLow );
maskDOSbg = imerode(maskDOSbgRaw,morphShape);

if(flagDoNotErode)
    maskDOSeroded = maskDOS;
    maskDOSbg = maskDOSbgRaw;
end

if(sum(maskDOSeroded(:)) == 0 )
    error('No data within cautious DOS mask')
   else if(sum(maskDOSbg(:)) == 0 ) 
           error('No region could be found to define an DOS background')
        end
end

% Try subtracting background;
imAmDxDOS = imAmDxDOS-mean(imAmDxDOS(maskDOSbg));
imDmDxDOS = imDmDxDOS-mean(imDmDxDOS(maskDOSbg));

% Create an image matrix of the Donor Emission Ratio
imRatioDOS = zeros(size(imDmDxDOS));
imRatioDOS(maskDOS) = imAmDxDOS(maskDOS)./imDmDxDOS(maskDOS);

% Find the mean DER, e.g. DER = 0.75;
DER = mean(imRatioDOS(maskDOSeroded));

DER % Echo DER to MATLAB


%%
% 2. Process some reference images to determine calibration numbers
% "LIN" refers to an image taken using a Linker Construct - 
% But other references should work provided their dFRET and aFRET are equal 
% to the "KnownDFRET" parameter specified above.

[FileNameAmAxLIN,PathName] = uigetfile('*.tif',...
    'Open Linker Sample, Acceptor Emission, Acceptor Excitation');  
imAmAxLIN = double( imread([PathName FileNameAmAxLIN]) );

[FileNameAmDxLIN,PathName] = uigetfile('*.tif',...
    'Open Linker Sample, Acceptor Emission, Donor Excitation');  
imAmDxLIN = double( imread([PathName FileNameAmDxLIN]) );

[FileNameDmDxLIN,PathName] = uigetfile('*.tif',...
    'Open Linker Sample, Donor Emission, Donor Excitation'); 
imDmDxLIN = double( imread([PathName FileNameDmDxLIN]) );

% Mask and background subtraction
% Bug fixed here - check name consistency!
% Eroded mask for setting Alpha and Beta is within Low AND High Thresholds
maskLIN = ( imAmAxLIN > ThresholdLow & imAmAxLIN < ThresholdHigh & ...
            imDmDxLIN > ThresholdLow & imDmDxLIN < ThresholdHigh & ...
            imAmDxLIN > ThresholdLow & imAmDxLIN < ThresholdHigh );
maskLINeroded = imerode(maskLIN,morphShape);

% Mask for background determination starts under the Low Threshold ONLY
maskLINbgRaw = ( imAmAxLIN < ThresholdLow & ...
              imDmDxLIN < ThresholdLow & ...
              imAmDxLIN < ThresholdLow );
maskLINbg = imerode(maskLINbgRaw,morphShape);

if(flagDoNotErode)
    maskLINeroded = maskLIN;
    maskLINbg = maskLINbgRaw;
end


% Use bg mask to estimate and subtract background levels
imAmAxLIN = imAmAxLIN - mean(imAmAxLIN(maskLINbg));
imAmDxLIN = imAmDxLIN - mean(imAmDxLIN(maskLINbg));
imDmDxLIN = imDmDxLIN - mean(imDmDxLIN(maskLINbg));

% % Iterate FRET calculations to determine Alpha and Beta
% Alpha = 0.5; % This is an initial value, and will be iteratively adjusted
% Beta = 2;    % So is this

% Find "lcFRET" which refers to "Linker Corrected FRET"
lcFRET = imAmDxLIN - AER*imAmAxLIN - DER*imDmDxLIN;
lcFRET(lcFRET<0) = 0; % cFRET < 0 is unphysical

% Fit Alpha and Beta parameters
% The following step may go wrong in the presence of strong noise
% median(listOfAlphas) etc "should" be robust versus outliers...
% But this cannot be guaranteed. User should check plausible values found
% 

listOfAlphas = (DER*imDmDxLIN(maskLINeroded) ) ./ ...
               ( lcFRET((maskLINeroded))*(KnownDFRET^-1 -1) );

listOfBetas = lcFRET(maskLINeroded) ./ ...
              ( KnownDFRET*AER*imAmAxLIN(maskLINeroded) );

% Reject values of Alpha = inf, and Beta = 0, from lcFRET = 0...
listOfAlphas(listOfAlphas == inf) = [];
listOfBetas( listOfBetas  == 0)   = [];

% Hopefully this section has determined sensible Alpha and Beta values
% For the 1:1 linker, aFRET is the same as dFRET here - to determine Beta

Alpha = median(listOfAlphas(:) );
Beta  = median(listOfBetas(:) );

if( Alpha < 0 || Beta < 0 )
   error('Negative Alpha or Beta fitted to Linker image')
end

Alpha  % Echo alpha to MATLAB console
Beta   % Echo beta to MATLAB console


% End of calibration
% Now the stored values of AER, DER, Alpha and Beta 
% Can be used to process real specimen data
%
% 22 Jan 2013
% This script has been simplified - for ease of reading and basic use
% 1. It no longer processes batches of files to produce consistency stats
% 2. It no longer opens files using the Leica SP5 name pattern
%    Instead, users are promted for all files (tedious, but simple)
% 3. It solves for alpha, beta, instead of iterative fitting.
% 
% Users should check for reasonable values 
% For example, like AER 0.75, DER 0.5, Alpha 0.3, Beta 2.3
% 