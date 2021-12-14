function lambda = debyeLength(ne,Te,cgsUnit)
% Calculates the electron debye length for the plasma
%  lambda = debyeLength(ne,Te,unit)
%
%  Symbol           Description                    SI(default) CGS
%  ne               electron density               [m-3]       [cm-3]  
%  Te                temperature                    [eV]        [eV]
%  
%  cgsUnit(optional)setting to 1 specifies CGS
%
%  nc               debye length                   [m]         [cm]  
%
% Written by Swadling Jan 2017
lambda = sqrt(SI.e0*SI.kbeV.*Te./(ne.*SI.e^2));
if nargin ==3 && cgsUnit ==1
    lambda = lambda.*1e-1; 
end
end
