
% Add all subdirectories of the parent directory of this
% script into the path

p0 = which('nmppath');
K = strfind(p0, filesep);
p1 = p0(1:K(end)-1);
mypath = genpath(p1);
path(path,mypath);

%% Get root dir for numipad
K = strfind(p1, filesep);
nmprootpath = p1(1:K(end)-1);


if( strcmp(p1(K(end-1)+1:K(end)-1),'numipad') )  % matlab interface for NUMIPAD has been installed
  nmprootdatapath = nmprootpath;  
else                                                  % matlab interface is beeing run locally
  K = strfind(nmprootpath, filesep);
  nmprootdatapath = strcat(nmprootpath(1:K(end)-1),'/tests');
  nmprootpath = nmprootpath(1:K(end)-1);
end

%  clear p0 K mypath p1
clear p0 K mypath p1 

