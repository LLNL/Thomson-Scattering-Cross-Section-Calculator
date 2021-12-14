% This is a demonstration script for running TS calculations on a
% self-refining wavelength scale. This is useful for properly captuing the
% amplitude of resonance in high alpha plasma.


clear;
close all;

% Set the wavlength range for the calculation
wavelengthRangeNM = [280  320];
% Set the insturment width
instrumentWidthNM = 0.291;
% Set the scattering information
probeWavelengthNM = 351;
scatteringAngle = 40;
% Make a list of electron density and temperature for this demo run
neList = [1.1  2.9 1.6 4].*1e26;
TeList = [1.4 1.3 0.83 0.5].*1e3;
% Here I run with two ion populations  (CH2 plasma)- this is extendable 
A = [12 1];
Z = [6 1];
fraction = [1 2];
fraction = fraction./sum(fraction);
Vi = [0 0];
J = 0;

% make a figure to display results
fig = figure('Color','white');
ax =axes;
lineColors =  ax.ColorOrder;
for ind = 1:length(neList)
    %% Run the calcuation
    % Parameters
    Te = TeList(ind);
    ne = neList(ind);
    Ti = [1 1].*Te;    
    
    % Calculation of cross section using mech refinement:
    [wavelength,spectrum,alpha,validity] = dynamicThomsonCrossSection(wavelengthRangeNM,probeWavelengthNM,ne,Te,fraction,Ti,A,Z,Vi,J,scatteringAngle,2);
    
    % The spectrum is calculated on a self-refining wavelength scale
    % This must be rebinned to a regular wavelength grid before applying
    % broadening
    minResNM = instrumentWidthNM/3;
    if minResNM == 0
        minResNM = 0.001;
    end
    [rebinedWavelength,rebinedSpectrum] = rebinData(wavelength,spectrum,wavelengthRangeNM(1), wavelengthRangeNM(2),minResNM);
    
    % Apply broadening
    rebinedSpectrum = gaussianBroadening((rebinedWavelength),rebinedSpectrum,instrumentWidthNM);
    
    % Plot spectra. Validity output is used to detect if peaks were
    % properly resolved. In very collective cases the peak height will not
    % be resolved and the feature will not have the proper amplitude.
    if validity ==1
        lineStyle = '-';
    else
        lineStyle = ':';
    end

    % Plot results:
    % Raw spectrum
    plot(wavelength,spectrum,'LineStyle',lineStyle,'Color',lineColors(ind,:));
    hold on;
    % After rebinning and broadening
    l(ind) = plot(rebinedWavelength,rebinedSpectrum,'LineStyle',lineStyle,'Color',lineColors(ind,:),'LineWidth',2);
    hold on;
    
    leg{ind} = ['ne = ' num2strScientific(ne,2) 'm^{-3} Te = ' num2strScientific(Te,3,3) ' eV \alpha = ' num2strScientific(alpha,2)];
end
ax.YLim = [0 1e-4];
xlabel('Wavelength [nm]');
ylabel('scatteringCrossSection [m^-1  nm^-1 sr^-1]');
legend(l,leg);
