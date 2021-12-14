function nc = criticalDensity(lambda,cgsUnit)
%% Calculates the critical density for optical propagation in a plasma
%  nc = criticalDensity(lambda,unit)
%
%  Symbol           Description                    SI(default) CGS
%  lambda           wavelength                     [m]        [cm]
%  cgsUnit(optional)setting to 1 specifies CGS
%
%  nc               critialDensity                 [m^-3]     [cm^-3]  
%
% Written by Swadling Jan 2017

nc = pi./(lambda.^2*SI.re);
% Check for CGS unit choice
if nargin ==2 && cgsUnit ==1;   
    nc = nc.*1e-2; 
end
end
