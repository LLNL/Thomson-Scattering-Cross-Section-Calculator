function nu = refractiveIndex(ne,lambda,cgsUnit)
%% Calculates the refractive index for optical propagation in a plasma
%  nc = refractiveIndex(ne,lambda,unit)
%
%  Symbol           Description                    SI(default) CGS
%  ne               electrondensity                [m^-3]     [cm^-3]  
%  lambda           wavelength                     [m]        [cm]
%  cgsUnit(optional)setting to 1 specifies CGS
%
%  nu               refractive index               [m^-3]     [cm^-3]  
%
% Written by Swadling Jan 2017

if nargin<3
    cgsUnit=0;
end
    nu = sqrt(1-ne/criticalDensity(lambda,cgsUnit));
end