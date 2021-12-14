function [str,textString] = num2strScientific(num,sigFigs,sciNotThreshold)

if nargin <2
    sigFigs=3;
end
if nargin<3
    sciNotThreshold=1;
end
% [str] = num2strScientific(num,sigFigs)
%
% This function is a scientific notation version of num2str.
% num is a number (or array of numbers)
% sigFigs is the number of siginficant figures
% sciNotThreshold may be used to set the digit limit at which scientific notation is
% used.
% By George 2/2018




if numel(num) ==1
    num = round(num,sigFigs,'significant');
    exp = real(floor(log10(num)));
    decimal = num./(10.^exp);
    if (abs(exp) <= sciNotThreshold) || isinf(exp)
        str = num2str(num);
        
    else
        str = [num2str(decimal,['%.' num2str(sigFigs-1) 'f']) ' \times 10^{' num2str(exp) '}'];
    end
    textString = strrep(str,'\times','x');
    textString = strrep(textString,'}','');
    textString = strrep(textString,'{','');
else
    %{
    for ind = 1: length(num)
        if exp(ind) <= limit && exp(ind)>=0
            str = num2str(num);
        else
            str{ind} = [num2str(decimal(ind)) ' \times 10^{' num2str(exp(ind)) '}'];
        end
    end
    
    %}
    s=size(num);
    num = num(:);
    for ind = 1: length(num)
        
        [str{ind}, textString{ind}]= num2strScientific(num(ind),sigFigs,sciNotThreshold);
        
    end
    str = reshape(str,s);
    textString = reshape(textString,s);
end

