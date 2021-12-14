function [wavelength,scatteringCrossSection,alpha,validity]  ...
    = dynamicThomsonCrossSection(wavelengthRange,probeWavelength,neTotal,Te,fract,Ti,A,Z,vI,J,scatAngle,order,superGaussian)
% [wavelength,scatteringCrossSection,alpha,validity]  ...
%    = dynamicThomsonCrossSection(wavelengthRange,probeWavelength,neTotal,Te,fract,Ti,A,Z,vI,J,scatAngle,o
if ~exist('superGaussian','var'); superGaussian =2.0;end

% These are the parameters driving the "dynamic" fiting algorithm
wavelengthResolutionThreshold = 1e-4;
angleThreshold =1e-3;
validity=1;
% Define inital resolution
res = 10;
%% Clean up the user input and format data for calculation
if abs(sum(fract)-1)>1e-2
    warning('Input fractions have been normalized');
    fract = fract./sum(fract);
end
Ti = double(Ti(:)');fract=double(fract(:)');A =double(A(:)');Z=double(Z(:)');vI = double(vI(:)');
neTotal = double(neTotal); Te = double(Te); J = double(J); scatAngle = double(scatAngle);
%% Determine populated populations, discard unpopulated inputs.
pop= fract > 0; fract=fract(pop);
Ti=Ti(pop);A=A(pop);Z=Z(pop);vI =vI(pop);

%% Calculate the "individual" population densities and associated partial electron densities
% Mean ionization state:
Zbar = sum(fract.*Z);
% Total ion density
niT = neTotal./Zbar;
% Population densities
ni =fract.*niT;
% Z weighted electron densities per population
ne = ni.*Z;

% Set the scattering information
probe = Photon(probeWavelength);
% Care for shoddy inputs
if length( wavelengthRange)>2
    wavelengthRange = minMax(wavelengthRange);
    warning('size of wavelength range is greater than 2. Using min and max as range. ')
end

% calculate the inital wavlength range
wavelength = wavelengthRange(1)-res:res:wavelengthRange(2)+res;

% calculate the inital cross-section. This is the starting point for
% optimization.

[scatteringCrossSection,~,alpha] = coreThomsonCrossSection(Photon(wavelength),probe,neTotal,ne,Te,Ti,A,Z,vI,J,scatAngle,order,superGaussian);
logScatteringCrossSection = real(log10(scatteringCrossSection));
% escape is initialized zero - this is used to escape the loop
escape =0;

%% Smoothing loop - This loop looks for discontinuities in the data which exceed the small angle threshold defined at the top of the page
% Iteration counter...
iterations=0;
while escape==0
    % Increment the iteration counter
    iterations=iterations+1;
    % Emergency exit...
    if iterations>100
        warning('Loop terminated after 100 iterations');
        escape = 1;
    end
    % Find the wavelength step-sizes
    deltaWavelength = diff(wavelength);
    % Find the difference in "angle" between adjacent parts of the log10(cross-section)curve
    angle = atan2(diff(logScatteringCrossSection),deltaWavelength);
    diffAngle = (diff(angle));
    % Find the wavelength that coresponds to each of these changes in "angle"
    wavelengthAngle =movmean(wavelength,2,'Endpoints','discard');
    % Find the wavelength step between adjacent angle changes
    deltaWavelengthDiffAngle =diff(wavelengthAngle);
    sepTest=(deltaWavelengthDiffAngle>wavelengthResolutionThreshold);
    % Find the angle changes (the angle is larger than the smallAngleThreshold seperation is greater than the wavlengthResolutionThreshold)
    changes =  uint32(find(abs(diffAngle)>angleThreshold  & sepTest ));
    %% Refine the fit around the detected areas where the angle change is large.
    if isempty(changes)
        % if no changes then exit
        escape = 1;
    else
        % Calculate new-wavelength values above and below the detected
        % position. The new points lie half way between old point pairs.
        
        % Find the unique change points. Changes+1 is needed as we need to
        % add points on both sides of the "angle" which is too steep.
        changes = unique([changes changes+1]);
        % Get a list of new wavelength points
        newWavelength = (wavelength(changes)+wavelength(changes+1))./2;
        % Calculate the coresponding wavelength values
        newScatteringCrossSection  = coreThomsonCrossSection(Photon(newWavelength),probe,neTotal,ne,Te,Ti,A,Z,vI,J,scatAngle,order,superGaussian);
        % Get the length of the current wavelength
        n= numel(wavelength);
        % Get the index where changes should be inserted
        changes= changes+1;
        % Find where the old wavelength value positions should be
        logical = false(1,n); logical(changes)=1;
        oldValuePositions = cumsum(logical)+(1:n);
        % Find where the new wavelength value positions should be
        newValuePositions = changes+uint32((0:length(changes)-1));
        % Indexing!
        wavelength([oldValuePositions newValuePositions]) = [wavelength newWavelength];
        scatteringCrossSection([oldValuePositions newValuePositions]) =[scatteringCrossSection newScatteringCrossSection];
        logScatteringCrossSection([oldValuePositions newValuePositions]) = [logScatteringCrossSection real(log10(newScatteringCrossSection))];
    end
end

%% Peak Finding loop - refines around peaks in the spectrum

% Find the peak locations:
% get the difference between adjacent peaks
Diff = diff(real(scatteringCrossSection));
% find Positive gradients
indexOver = Diff>0;
%find Negative gradients
indexUnder = Diff<0;
% Find peaks - discard any peaks appearing in negative values of the
% cross-section - these dont behave well...
indexPeak = [false indexOver] & [indexUnder false] & (scatteringCrossSection>0);
indexPeak =find(indexPeak);

% Loop over each peak
for ind  = 1: length(indexPeak)
    % Get the indices of the three points making up the "peak".
    peakIndexs = [-1:1]+indexPeak(ind);
    % Get the values
    threePointPeakWavelength = wavelength(peakIndexs);
    threePointPeakSpectrum = scatteringCrossSection(peakIndexs);
    % initialize the data store with the current center value of the peak.
    newStoreWavelength= threePointPeakWavelength(2);
    newStoreScatteringCrossSection = threePointPeakSpectrum(2) ;
    % initialize escape parameters
    escape = 0;
    tries = 0;
    warn =1;
    while escape == 0      
        % Save the current max value
        oldMax = threePointPeakSpectrum(2);
        % Calculate the new wavelength values with will be calculated
        newWavelength = [threePointPeakWavelength(1)+threePointPeakWavelength(2) threePointPeakWavelength(2)+threePointPeakWavelength(3)]./2;
        
        % And the coresponding TS values
        newScatteringCrossSection  = coreThomsonCrossSection(Photon(newWavelength),probe,neTotal,ne,Te,Ti,A,Z,vI,J,scatAngle,order,superGaussian);
       
        % Make arrays of the new 5 point peak
        fivePointPeakWavelength = [threePointPeakWavelength(1) newWavelength(1) threePointPeakWavelength(2) newWavelength(2) threePointPeakWavelength(3)];
        fivePointPeakSpectrum = [threePointPeakSpectrum(1) newScatteringCrossSection(1) threePointPeakSpectrum(2) newScatteringCrossSection(2) threePointPeakSpectrum(3)];       
        dW = diff(fivePointPeakWavelength);    
        if any(dW==0)
            % Escape if the deltas have colapsed to zero -this happens due
            % to limit of double precisison
            escape =1;
            if warn == 1
                warning('Peak may not be resolved. Peak cannot be found within Double precision.');
                validity=0;
            end
        else
            % Add the new wavelength values to the data store
            newStoreWavelength = [newStoreWavelength newWavelength];
            newStoreScatteringCrossSection = [newStoreScatteringCrossSection newScatteringCrossSection];
                      
            % Get the three-point sub peak within the 5 point peak
            [newMax, I] = max(real(fivePointPeakSpectrum));
            threePointPeakWavelength = fivePointPeakWavelength((-1:1)+I);
            threePointPeakSpectrum = fivePointPeakSpectrum((-1:1)+I);
            % Run tests for warning supression
            if abs((newMax-oldMax)./oldMax)<0.01
                tries = tries+1;
                % If the new max is more than 5% higher than the previous then reset tries             
            else % else increment
                tries = 0;
            end
            if tries>10 
                % If ten tests in a row are passed then the fit is assumed
                % to be good.
                warn =0;
            end
        end
        
    end
    % Sort the peak wavelength and spectum
    [newStoreWavelength,sortInd] = sort(newStoreWavelength);
    newStoreScatteringCrossSection = newStoreScatteringCrossSection(sortInd);
    % insert the data
    wavelength = [wavelength(1:peakIndexs(1)) newStoreWavelength wavelength(peakIndexs(3):end)];
    scatteringCrossSection = [scatteringCrossSection(1:peakIndexs(1)) newStoreScatteringCrossSection scatteringCrossSection(peakIndexs(3):end)];
    % Increase the index to take account of the points which have been
    % added.
    indexPeak = indexPeak+ length(newStoreWavelength)-1 ;
end
wavelength = wavelength(:);
scatteringCrossSection=scatteringCrossSection(:);

end