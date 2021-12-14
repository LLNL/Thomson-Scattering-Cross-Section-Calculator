function [broadenedSpectrum] = gaussianBroadening(spectralRange,intensity,sigma)
% Calculates the particle inertial length to the spectrometer resolution.
%
%  [broadenedSpectrum] = gaussianBroadening(lambda,spectrum,width)
%  This function convolutes the spectrum with a gaussian width correpsonding
%  Symbol           Description                  
%
%  lambda                wavelength scale           [arbitary]
%  spectrum              Atomic Mass                [arbitary]
%  sigma          gaussian standard deviation       [same as lambda]
%  FWHM = 2*sqrt(2*ln(2))*sigma = 2.35*sigma
%
%  broadenedSpectrum     spectrum with broadening   [same as spectrum]
%
% Swadling Feb 2017

arguments
    spectralRange = [1:1000]
    intensity = [ones((size(spectralRange)./[1 2])-[0 1]) 500 ones(size(spectralRange)./[1 2])];
    sigma = 10;
end

s =size(spectralRange);
spectralRange = spectralRange(:);
intensity = intensity(:);


if isa(spectralRange,'Photon')
    lambda = spectralRange.wavelengthNM;
else
    lambda = spectralRange;
end


if sigma ==0
    broadenedSpectrum =intensity;
else
    wavelengthResolution = abs(lambda(2)-lambda(1));
    % Find range to calculate spectrometer function over
    if sigma >= wavelengthResolution
        padLength = round((3.*sigma)/wavelengthResolution);
        lambda_Spec = -padLength*wavelengthResolution:wavelengthResolution:wavelengthResolution*padLength;
        
        % Calulate spectrometer function (norm. gaussianGaussian)
        f = exp(-0.5.*(lambda_Spec./sigma).^2);
        f = f./sum(f);

        paddedIntensity = [ones(padLength,1).*intensity(1);intensity;ones(padLength,1).*intensity(end)];
        % Convulute the spectrum with the calculated gaussian  - retain only the
        % vaid range
        broadenedSpectrum = conv(paddedIntensity,f,'valid')./trapz(f);
       % broadenedSpectrum (1:length(f)) =nan;
        %broadenedSpectrum (end-length(f):end)=nan;
    else
        broadenedSpectrum = intensity;
    end
end
if nargout == 0 
    
    figure;
    plot(spectralRange,intensity);
    hold on;
    plot(spectralRange,broadenedSpectrum);
    wait = 0;
end

end