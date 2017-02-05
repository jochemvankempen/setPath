function outNumString = removePeriodFromName(inNum)
% this function removes the period from names, useful for creating
% figures that cannot have periods in the filename. (eg. for texmaker).
% with the exception of the file extensions defined by fileformats
%
% INPUT
% inNum, string or number. With or without period. variable or full
% filename.
% 
% OUTPUT
% outNumString, string without period

fileformats = {'mat','png','eps','pdf','fig','jpeg','jpg'};

if isnumeric(inNum)
    outNumString = num2str(inNum);
else
    outNumString = inNum;
end

idx = strfind(outNumString,'.');
if ~isempty(idx)
    for iIdx = 1:length(idx)
        if length(outNumString)>=(idx(iIdx)+3)
            [~,fileFormatIdx1] = grep(fileformats,outNumString(idx(iIdx):idx(iIdx)+3));
        end
        if length(outNumString)>=(idx(iIdx)+4)
            [~,fileFormatIdx2] = grep(fileformats,outNumString(idx(iIdx):idx(iIdx)+4));
        end
        
        if (~exist('fileFormatIdx1','var') && ~exist('fileFormatIdx2','var')) || (isempty(find(fileFormatIdx1,1)) && isempty(find(fileFormatIdx2,1)))
            outNumString(idx(iIdx)-(iIdx-1))=[];
        end
    end
end

