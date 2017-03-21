function [data, savefilename] = prepareSPGdata(dimension,experiment,Channels2check,iChan,iSession,subID,varargin)


switch experiment.name
    case 'CFS_localizer'
        [data, savefilename] = prepareSPG_FaceLocalizer(dimension,experiment,Channels2check,iChan,iSession,subID,varargin{1});
    case 'CFS'
        [data, savefilename] = prepareSPG_CFS(dimension,experiment,Channels2check,iChan,iSession,subID);
    case 'CFSRT'
        [data, savefilename] = prepareSPG_CFSRT(dimension,experiment,Channels2check,iChan,iSession,subID,varargin);
        
end


