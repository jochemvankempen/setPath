function copyToKani(localDir,filename)
% copies file on local machine to remote server
%
% INPUT:
%
% localFile = full directory with filename
%   e.g. /home/vjochem/Monash052/jochem/'filename.mat'
% localDir  = directory on remote server
%   e.g. /export/kani/jochem/
%
% Jochem van Kempen 13-10-14

javaaddpath('/gpfs/M2Home/projects/Monash052/jochem/repositories/scripts/setPathFunctions/ganymed-ssh2-build250.jar')

kaniDir = ['/export/kani/jochem/'];

% if isdir('/gpfs/M2Home/projects/Monash052/jochem/')
%     if strncmp(localFile,'/gpfs/M2Home/projects/Monash052/jochem/',length('/gpfs/M2Home/projects/Monash052/jochem/'))
remoteDir = [kaniDir localDir(length('/gpfs/M2Home/projects/Monash052/jochem/')+1:end)];
%     end
% elseif isdir('/home/vjochem/M1/')%M1getDecodeSettings_TW_Iowa
% end

% scptomatlab(userName,hostName,password,localFolder,remotefilename)
scptomatlab('jochem','130.194.90.116','jochem123',...
    [localDir filesep filename],remoteDir)
