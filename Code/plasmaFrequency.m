function w = plasmaFrequency(n,A,Z,cgsUnit)
% Calculates the particle plasma frequency - units in s-1
%  w = plasmaFrequency(ne,A,Z,cgsUnit)
%
%  Symbol           Description                    SI(default) CGS
%
%  n                electron/ion density           [m-3]       [cm-3]
%  A                Atomic Mass                    [AU set to 0 for electron]
%  Z                Ionization State               [-1 for electron]
%  cgsunit(optional)setting to 1 specifies CGS
%
%  w                plasma angular frequency       [s-1]       [s-1]
%
% Written by George Swadling Jan 2017

% Default to SI electron value
if nargin == 1
    A=0;
    Z=-1;
end
% Calculate particle mass
m =A.*SI.u;
m(m==0)=SI.me;
% Calculate angular frequency
w = sqrt(SI.e^2.*n.*Z^2./(m.*SI.e0));
if nargin ==4 && cgsUnit ==1
    w = w.*1e3;
end
end