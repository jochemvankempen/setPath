function copyToKani_M2(localDir,filename)
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

javaaddpath('/gpfs/M1Scratch/Monash052/repositories/scripts/setPathFunctions/ganymed-ssh2-build250.jar')

kaniDir = ['/export/kani/jochem/'];

% if isdir('/gpfs/M2Home/projects/Monash052/jochem/')
%     if strncmp(localFile,'/gpfs/M2Home/projects/Monash052/jochem/',length('/gpfs/M2Home/projects/Monash052/jochem/'))
remoteDir = [kaniDir localDir(length('/gpfs/M2Home/projects/Monash052/jochem/')+1:end)];
%     end
% elseif isdir('/home/vjochem/M1/')%M1getDecodeSettings_TW_Iowa
% end

% scptomatlab(userName,hostName,password,localFolder,remotefilename)
scptomatlab('vjochem','m2.massive.org.au','1JPcoen2',...
    [localDir filesep filename],remoteDir)
