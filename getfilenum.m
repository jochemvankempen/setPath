function filenum = getfilenum(session,nOrder)
if nargin < 2 
    nOrder = 4;
end
filenum = num2str(session);
while length(filenum)<nOrder
    filenum = ['0' filenum];
end;  %add zeros at the beginning
