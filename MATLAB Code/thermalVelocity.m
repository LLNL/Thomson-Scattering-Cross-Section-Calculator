function V = thermalVelocity(T,A,cgsUnit)
% Calculates the particle thermal velocity
%  V = thermalVelocity(T,A,unit)
%
%  Uses the formula: Vth = sqrt(2kT/m)
%  Most Probable speed in three dimentions
%
%  Symbol           Description                 SI(default) CGS
%  T                temperature                     [eV]        [eV]
%  A                Atomic mass of                  [Atomic Units]
%                   particle species                [A=0 for electron]
%  cgsUnit(optional)setting to 1 specifies CGS
%
%  V                electron thermal velocity       [ms-1]      [cms-1]
%
% Written by Swadling Jan 2017

% Convert atomic mass to SI mass. If A == 0 use electron mass;

% Default to SI electrons for single input
if nargin ==1
    A=0;
end

m = A.*SI.u;
m(A==0)=SI.me;
V = sqrt(2.*SI.kbeV.*T./m);
if nargin ==3 && cgsUnit ==1
    V = V.*1e2;
end
end
