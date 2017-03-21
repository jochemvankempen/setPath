function plotOnMap(subID,PLOT,savefilename,DIR,channelValues,colorScale)
%Plots any values for either uni- , bipolar or all channels on subject
%specific maps.
%INPUT:
%subID          = a cell with the subject name as a string
%PLOT           = Structure array with fields
%                   -visible ('on'/'off')
%                   -print (print figure in DIR.fig)
%                   -channels ('bipolar','unipolar','all')
%savefilename   = name of the figure you want to save, this script will
%   add ventral/temporal/frontal to this name for the different images
%DIR            = struct array with fields
%                   -fig (figure directory)
%channelValues  = vector with values for all channels for that subject
%colorScale     = color limits (by default minimum-maximum of
%   channelValues)
%
%Jochem van Kempen

subSpecsIowaFaceLocalizer

if ~isinteger(bCNT.img_temp)
    bCNT.img_temp = uint8(bCNT.img_temp);
    bCNT.img_vent = uint8(bCNT.img_vent);
end

switch subID
    case '153'
        bCNT.img_vent(430:512,380:444,:)=255;
        bCNT.img_temp(290:400,1:80,:)=255;
end

%
if ~isfield(PLOT,'alpha')
    PLOT.alpha = 1;
end


chanIdx = getChannelLocation_Iowa(bCNT,cntID);

if isempty(colorScale) 
    scaleDecode = [min(channelValues) max(channelValues)];
else
    scaleDecode = colorScale;
end


switch PLOT.channels
    case 'unipolar'
        ventChan2use    = chanIdx.ventral;
        tempChan2use    = chanIdx.temp;
        frontChan2use   = chanIdx.frontal;
    case 'bipolar'
        ventChan2use    = chanIdx.ventralBipolar;
        tempChan2use    = chanIdx.tempBipolar;
        frontChan2use   = chanIdx.frontalBipolar;
    case 'all'
        ventChan2use = [chanIdx.ventral ventChanBipIdx];
        tempChan2use = [chanIdx.temp tempChanBipIdx];
        frontChan2use = [chanIdx.frontal chanIdx.frontalBipolar];
end


%% ventral image
figureSettings
clear colormap cmap h imHandle
% colormap('hot')
r = 7;
theta = linspace(0,2*pi,16);
R = r*[cos(theta); sin(theta)];

colormap('default');

imHandle = imagesc(bCNT.img_vent);
if ~isempty(PLOT.alpha)
    alpha(imHandle,PLOT.alpha)
end
% imAdjust = imadjust(bCNT.img_vent,[0 0.5],[0 1]);
% imagesc(imAdjust)
axis image
axis off
hold on
for iChan = 1:length(ventChan2use)
    colval = (channelValues(ventChan2use(iChan))-scaleDecode(1))/diff(scaleDecode)*255;
    patch(bCNT.xcent(ventChan2use(iChan))+R(1,:),bCNT.ycent(ventChan2use(iChan))+R(2,:),colval)
end

% h = colorbar;
% set(get(h,'xlabel'),'String', 'A''','FontSize',20);
% set(h,'ylim',[0 255],'ytick',linspace(0,255,5),'yticklabel',num2str(linspace(scaleDecode(1),scaleDecode(2),5)',2),'FontSize',15);

h = colorbar('peer',gca);
set(h,'ylim',[0 255],'ytick',linspace(0,255,5),'yticklabel',num2str(linspace(scaleDecode(1),scaleDecode(2),5)',2),'FontSize',15,'colormap','Jet');
set(get(h,'xlabel'),'String', 'A''','FontSize',20);
colormap('Jet');

if PLOT.print
    tmpSavefilename = savefilename;
    savefilename = [savefilename '_ventral'];
    figureSave
    savefilename = tmpSavefilename;
else 
    pause
end
% set(c,'position',[0.84 0.05 0.1 0.9]);
% keyboard

%% temporal image
figureSettings
clear colormap cmap h imHandle
imHandle = imagesc(bCNT.img_temp);
if ~isempty(PLOT.alpha)
    alpha(imHandle,PLOT.alpha)
end
colormap('default');
axis image
axis off
hold on
for iChan = 1:length(tempChan2use)
    colval = (channelValues(tempChan2use(iChan))-scaleDecode(1))/diff(scaleDecode)*255;
    patch(bCNT.xcent(tempChan2use(iChan))+R(1,:),bCNT.ycent(tempChan2use(iChan))+R(2,:),colval)
    
end
h = colorbar;
set(get(h,'xlabel'),'String', 'A''','FontSize',20);
set(h,'ylim',[0 255],'ytick',linspace(0,255,5),'yticklabel',num2str(linspace(scaleDecode(1),scaleDecode(2),5)',2),'FontSize',15);
% set(h,'Location','South')
if PLOT.print
    tmpSavefilename = savefilename;
    savefilename = [savefilename '_temporal'];
    figureSave
    savefilename = tmpSavefilename;
else 
    pause
end
%% potential frontal image
if isfield(bCNT,'img_front')
    figureSettings
    clear colormap cmap h imHandle
    
    imHandle = imagesc(bCNT.img_front);
    alpha(imHandle,PLOT.alpha)
    colormap('default');
    axis image
    axis off
    hold on
    for iChan = 1:length(frontChan2use)
        colval = (channelValues(frontChan2use(iChan))-scaleDecode(1))/diff(scaleDecode)*255;
        patch(bCNT.xcent(frontChan2use(iChan))+R(1,:),bCNT.ycent(frontChan2use(iChan))+R(2,:),colval)
        
    end
    h = colorbar;
    set(get(h,'xlabel'),'String', 'A''','FontSize',20);
    set(h,'ylim',[0 255],'ytick',linspace(0,255,5),'yticklabel',num2str(linspace(scaleDecode(1),scaleDecode(2),5)',2),'FontSize',15);
    % set(h,'Location','South')
    if PLOT.print
        tmpSavefilename = savefilename;
        savefilename = [savefilename '_frontal'];
        figureSave
        savefilename = tmpSavefilename;
    end
    
end
















