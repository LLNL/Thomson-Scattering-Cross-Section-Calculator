function [X,Y] = rebinData(x,y,minX,maxX,resolutionX)
%  [X,Y] = rebinData(x,y,minX,maxX,resolutionX)
%
% This function rebins data calculated on a non-linear scale onto a linear
% scale for ease of further processing.
%
% This function was developed for Thomson scattering analaysis 2/2018
%
% By George


x= x(:);y=y(:);

%% Linear interpolate data to exceed rebinning resolution:
inArea = trapz(x,y);
resolution = diff(x);
insertions = resolution>resolutionX;

while any(insertions)
    % Calculate new x values
    newx = (x([insertions;false])+ x([false;insertions]))./2;
    newy = interp1(x,y,newx);
    % Find where the new values should be inserted
    insertions = find(insertions)+1;
    n = numel(x);
    logical = false(n,1); logical(insertions)=1;
    oldValuePositions = cumsum(logical)+(1:n)';
    newValuePositions = insertions+((0:length(insertions)-1))';
    % Index the new and old values together in sequece
    x([oldValuePositions;newValuePositions]) = [x;newx];
    y([oldValuePositions;newValuePositions]) =[y;newy];
    % Recalculate for next loop:
    resolution = diff(x);
    insertions = resolution>resolutionX;
end

%% Rebin. This is an area preserving operation - this keeps data about the "brightness" of very narrow features.
X = [minX:resolutionX:maxX];
XEdges = (minX-resolutionX/2):resolutionX:(maxX+resolutionX/2);
% Calculate the cumulative trapezium integral of the data
cumulativeSum = cumtrapz(x,y);
% Interpolate the cumsum onto the bin edges, take difference and divide
% by the bin width to find a bin-average
IcumSum = interp1(x,cumulativeSum,XEdges,'linear');
Y =  diff(IcumSum)./resolutionX;

X= X(:);Y=Y(:);
end



