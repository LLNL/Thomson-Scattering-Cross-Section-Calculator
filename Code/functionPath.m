function [path] = functionPath()
% Returns the path to the currently running function / script
list = dbstack();
caller = list(2);
location = which(caller.name);


path =fileparts(location);
end

