function [Z, Zp] = plasmaDispersionFunction(x)
% Calculate the the plamsa dispersion function and its derivative.
% The plasma dispersion function is the thermal correction for the plasma
% dispersion relation. The outputs are complex numbers.
%
% [Z, Zp] = Plasma_Dispersion_Function(x)
%
% Input:-
% x - ratio of the phase velocity (omega/k) to the thermal velocity
%
% Outputs:-
% Z - plasma dispersion function
% Zp - derivative with x of the plasma disperison function
%
% Uses external Faddeeva routine to solve the imaginary error function -
% Dawson(z) = sqrt(pi)/2 * exp(-z^2) * erfi(z)
% erfi(z) is the imaginary error function
%
% By George 4/2013

% Calculate the dispersion function for vector x
try
    Z = 1i*sqrt(pi)*Faddeeva_w(x);
catch
    % Use Built in Function as a backup
    Z = -2.*(1-2.*x.*dawson(x))+1i.*(-2.*sqrt(pi).*x.*exp(-x.*x));
    warning('Precompiled Feddeva_w function is not working - Using built in Dawson Integral - Performance will be substantially reduced')
end
%Z = 1i*sqrt(pi).*exp(-x.*x).*erfc(-1i.*x);
% Calculate the derivative of the dispersion function
Zp = -2.*(x.*Z+1);
end