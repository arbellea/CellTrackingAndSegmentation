function varargout = TrackingUtility(varargin)
% TRACKINGUTILITY MATLAB code for TrackingUtility.fig
%      TRACKINGUTILITY, by itself, creates a new TRACKINGUTILITY or raises the existing
%      singleton*.
%
%      H = TRACKINGUTILITY returns the handle to a new TRACKINGUTILITY or the handle to
%      the existing singleton*.
%
%      TRACKINGUTILITY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKINGUTILITY.M with the given input arguments.
%
%      TRACKINGUTILITY('Property','Value',...) creates a new TRACKINGUTILITY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrackingUtility_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrackingUtility_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrackingUtility

% Last Modified by GUIDE v2.5 22-Feb-2016 13:21:31

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TrackingUtility_OpeningFcn, ...
    'gui_OutputFcn',  @TrackingUtility_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TrackingUtility is made visible.
function TrackingUtility_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrackingUtility (see VARARGIN)

% Choose default command line output for TrackingUtility
handles.output = hObject;
handles.buttonDownFlag = false;
handles.buttonDownMovingFlag = false;
handles.loadingImageFlag = false;
handles.buttonDownStaticFlag = false;
handles.currentImageNumber = 0;
%set(hObject,'isRunning',false);
userData.isRunning = false;
set(hObject,'userdata',userData);
%set(gcf, 'WindowButtonDownFcn', @getMousePositionOnImage);
%set (gcf, 'WindowButtonMotionFcn', @getMousePositionOnImage);
%set(src, 'Pointer', 'crosshair'); % Optional
pan off % Panning will interfere with this code
% Update handles structure

guidata(hObject, handles);

% UIWAIT makes TrackingUtility wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TrackingUtility_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliderMain_Callback(hObject, eventdata, handles)
% hObject    handle to sliderMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderPos = get(hObject,'Value');
userData = get(handles.figure1,'userdata');
framePos = round(sliderPos*(userData.imageData.Frame_Num-1))+1;
set(handles.pmFrameNum,'Value',framePos);
handles.currentImageNumber = framePos;
guidata(hObject,handles);
trackingData = userData.cellData(:,:,get(handles.pmCellID,'value'));
zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
if ~get(handles.tglSegment,'value')
    loadNextImage(handles.axMain,userData.imageData,trackingData,framePos,zoomAx);
else
    plotSegImage(handles,userData,framePos,get(handles.pmCellID,'value'))
end
if isempty(userData.segData{get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value')})
    set(handles.pbDrawClear,'string','Draw');
else
    set(handles.pbDrawClear,'string','Clear');
end

set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');



% --- Executes during object creation, after setting all properties.
function sliderMain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnNew.
function btnNew_Callback(hObject, eventdata, handles)
% hObject    handle to btnNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName, cancleFlag]  = uigetimagefile();
if cancleFlag
    return
end
[path,name,ext] = fileparts(fileName);
ext = ['*',ext];
expr = inputdlg('Enter the filen name regular expresion:',...
             'Expression', 1,{name});
 if isempty(expr{1})
     return
 end
imageData = Load_Data(path,ext,expr{1});
cellData = nan([imageData.Frame_Num,2,1]);
segData = cell(imageData.Frame_Num,1);
contourData = cell(imageData.Frame_Num,1);
maxID = 1;
strID = {'1'};
strColorID = {'1'};
prc = [0.5,99.5];
[saveFileLocation,savePath] = uiputfile('*.mat','Save File Location','SaveData.mat');
if savePath==0
    return;
end
saveFileLocation = fullfile(savePath,saveFileLocation);

I = imread(imageData.Frame_name{1});
set(handles.axMain,'Units','pixels');
set(gcf,'Units','pixels');
mainAxPos = get(handles.axMain,'Position');
figPos = get(gcf,'position');
handles.axZeroPoint = [figPos(1:2),0 0]+mainAxPos+[0,mainAxPos(4),0 0];
set(handles.axMain,'Units','normalized');
set(gcf,'Units','normalized');
p = prctile(double(I(:)),prc);
him = imshow(I,p,'Parent',handles.axMain);

handles.currentImageNumber = 1;
userData.isRunning = false;
userData.imageData = imageData;
userData.cellData = cellData;
userData.segData = segData;
userData.contourData = contourData;
userData.maxID = maxID;
userData.strID = strID;
userData.strColorID = strColorID;
userData.Link = struct([]);
userData.prc = prc;
save(saveFileLocation,'-struct','userData');
set(gcf,'userdata',userData);
set(handles.pmFrameNum,'string',1:imageData.Frame_Num,'value',1)
set(handles.pmCellID,'string',strColorID,'value',1)
set(handles.sliderMain,'SliderStep',[1/(imageData.Frame_Num-1),1/(imageData.Frame_Num-1)],'Value',0,'Min',0,'Max',1);

xlim(handles.axMain,[1 userData.imageData.Width]);ylim(handles.axMain,[1 userData.imageData.Height]);
guidata(hObject,handles);

set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');
TrackCell;


% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fileName,pathName] = uigetfile('*.mat','Load File');
if fileName==0, return; end;
userData = load(fullfile(pathName,fileName));
I = imread(userData.imageData.Frame_name{1});
sI = size(I);
% set(handles.axMain,'Units','pixels');
% set(gcf,'Units','pixels');
% mainAxPos = get(handles.axMain,'Position');
% mainAxPos(3) = mainAxPos(4)*sI(1)./sI(2);
% set(handles.axMain,'Position',mainAxPos);
% figPos = get(gcf,'position');
% sliderPos = get(handles.sliderMain,'Position');
% %figPos(1:2) = 0;
% %figPos(3:4) = figPos(3:4)+size(I)-mainAxPos(3:4);
% %mainAxPos(3:4) = size(I);
% %sliderPos(3) = mainAxPos(3);
% handles.axZeroPoint = [figPos(1:2),0 0]+mainAxPos+[0,mainAxPos(4),0 0];
% set(handles.axMain,'Units','normalized');
% set(gcf,'Units','normalized');
% %set(handles.axMain,'position',mainAxPos);
% %set(handles.sliderMain,'position',sliderPos);
% %set(gcf,'position',figPos);
p = prctile(double(I(:)),userData.prc);
him = imshow(I,p,'Parent',handles.axMain);
handles.currentImageNumber = 1;
TrackCell;


%{
userData = get(gcf,'userdata');
userData.imageData = imageData;
userData.cellData = cellData;
userData.maxID = maxID;
userData.strID = strID;
userData.strColorID = strColorID;
userData.segData = segData;
userData.contourData = contourData;
userData.Link = Link;
%}
set(gcf,'userdata',userData);
set(handles.pmFrameNum,'string',1:userData.imageData.Frame_Num,'value',1)
set(handles.pmCellID,'string',userData.strColorID,'value',1)
set(handles.sliderMain,'SliderStep',[1/(userData.imageData.Frame_Num-1),1/(userData.imageData.Frame_Num-1)],'Value',0,'Min',0,'Max',1);
%set(handles.axMain,'Units','normalized');
xlim(handles.axMain,[1 userData.imageData.Width]);ylim(handles.axMain,[1 userData.imageData.Height]);
%set(gcf,'Units','normalized');
%set(handles.sliderMain,'Units','normalized');
guidata(hObject,handles);
%set(him, 'ButtonDownFcn',{@axMain_ButtonDownFcn,guidata(gcf)});
%set(him, 'hittest','off');
set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');
TrackCell;

% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[saveFileLocation,savePath] = uiputfile('*.mat','Save File Location','SaveData.mat');
if savePath==0
    return;
end
userData = get(gcf,'userData');
saveFileLocation = fullfile(savePath,saveFileLocation);
save(saveFileLocation,'-struct','userData');
set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');

% --- Executes on button press in btnPlayPause.
function btnPlayPause_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of btnPlayPause

set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');

% --- Executes on button press in btnNewCell.
function btnNewCell_Callback(hObject, eventdata, handles)
% hObject    handle to btnNewCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1,'userdata');
cellMax = userData.maxID;
userData.maxID = cellMax+1;
userData.strID{end+1} = sprintf('%d',cellMax +1);
userData.strColorID{end+1} = sprintf('%d',cellMax +1);
[userData.segData{:,end+1}] = deal([]);
[userData.contourData{:,end+1}] = deal([]);
set(handles.pmCellID,'String',userData.strColorID,'value',numel(userData.strID));
userData.cellData(:,:,size(userData.cellData,3)+1)=nan;
set(handles.figure1,'userdata',userData);
currentFrame = get(handles.pmFrameNum,'value');
if currentFrame==userData.imageData.Frame_Num
    set(handles.pmFrameNum,'value',1);
    set(handles.sliderMain,'value',0);
    if ~get(handles.tglSegment,'value')
        loadNextImage(handles.axMain,userData.imageData,userData.cellData(:,:,end),1,[]);
    else
        plotSegImage(handles,userData,1,get(handles.pmCellID,'value'))
    end
    
end
set(handles.pbDrawClear,'string','Draw');




set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');


% --- Executes on selection change in pmCellID.
function pmCellID_Callback(hObject, eventdata, handles)
% hObject    handle to pmCellID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmCellID contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmCellID
userData = get(gcf,'userData');
selectedCell = get(hObject,'value');
cellTrack = userData.cellData(:,:,selectedCell);
if ~get(handles.tglSegment,'value')
    trackEnd = find(~isnan(cellTrack(:,1)),1,'last');
    if trackEnd < size(cellTrack,1)
        trackEnd = trackEnd+1;
    end
    if isempty(trackEnd)
        trackEnd = 1;
    end
else
    trackEnd = get(handles.pmFrameNum,'value');
end

userData.strColorID{selectedCell} = userData.strID{selectedCell};
set(gcf,'userData',userData);
set(hObject,'string',userData.strColorID);
set(handles.pmFrameNum,'value',trackEnd);
set(handles.sliderMain,'value',(trackEnd-1)*1./(userData.imageData.Frame_Num-1));
if ~get(handles.tglSegment,'value')
    loadNextImage(handles.axMain,userData.imageData,userData.cellData(:,:,selectedCell),trackEnd,[]);
    TrackCell;
else
    plotSegImage(handles,userData,get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value'))
end


if isempty(userData.segData{get(handles.pmFrameNum,'value'),selectedCell})
    set(handles.pbDrawClear,'string','Draw');
else
    set(handles.pbDrawClear,'string','Clear');
end
set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');


% --- Executes during object creation, after setting all properties.
function pmCellID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmCellID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pmFrameNum.
function pmFrameNum_Callback(hObject, eventdata, handles)
% hObject    handle to pmFrameNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pmFrameNum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmFrameNum
framePos = get(hObject,'Value');
userData = get(handles.figure1,'userdata');
sliderPos = (framePos-1)*1./(userData.imageData.Frame_Num-1);
set(handles.sliderMain,'Value',sliderPos);
handles.currentImageNumber = framePos;
guidata(hObject,handles);
trackingData = userData.cellData(:,:,get(handles.pmCellID,'value'));
zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
if ~get(handles.tglSegment,'value')
    loadNextImage(handles.axMain,userData.imageData,trackingData,framePos,zoomAx);
else
    plotSegImage(handles,userData,framePos,get(handles.pmCellID,'value'))
end


if isempty(userData.segData{framePos,get(handles.pmCellID,'value')})
    set(handles.pbDrawClear,'string','Draw');
else
    set(handles.pbDrawClear,'string','Clear');
end

set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');
% --- Executes during object creation, after setting all properties.
function pmFrameNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmFrameNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnCreateCSV.
function btnCreateCSV_Callback(hObject, eventdata, handles)
% hObject    handle to btnCreateCSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(gcf,'userData');
[saveDir] = uigetdir('','Please Choose Save Directory');
mkdir(fullfile(saveDir,'ManualSeg'))
table = [];
LinkData = zeros(size(userData.cellData,3),1);
for l = 1:numel(userData.Link)
    LinkData(userData.Link(l).Children) = userData.Link(l).Mother; 
end
time = repmat(shiftdim(1:size(userData.cellData,1),-1),1,size(userData.cellData,3),1);
cellData = permute(userData.cellData,[2,3,1]);
cellData = [repmat([1:size(userData.cellData,3)],1,1,size(cellData,3));cellData;time;repmat(LinkData',1,1,size(cellData,3))];
cellData2 = reshape(cellData,5,[])';
cellData2 = cellData2(~isnan(cellData2(:,2)),:);
fid = fopen(fullfile(saveDir,'ManualSeg','TrackingData.txt'),'w');
fprintf(fid,'cellID\tcenter_col\tcenter_row\ttime\tmother\n')
fclose(fid);
%(fullfile(saveDir,'ManualSeg','TrackingData.txt'){'cellID','center_col','center_row','time','mother'})
dlmwrite(fullfile(saveDir,'ManualSeg','TrackingData.txt'),cellData2,'-append','delimiter','\t')
return;
for i = 1:size(userData.cellData,1)
    frameData = squeeze(userData.cellData(i,:,:));
    validCells = ~isnan(frameData(:,1));
    MotherCell
    frameTable = cat(2,find(validCells),frameData(validCells,:),LinkData(find(validCells)));
    table = cat(1,table,frameTable);
end





%{
function getMousePositionOnImage(src, event)
handles = guidata(src);
while handles.buttonDownStaticFlag
    handles = guidata(src);
end
handles.buttonDownMovingFlag = true;
guidata(src,handles);
if handles.buttonDownFlag&&~handles.loadingImageFlag
    handles.loadingImageFlag = true;
    handles.buttonDownMovingFlag = true;
    handles.currentImageNumber = handles.currentImageNumber+1;
    guidata(src,handles)
    cursorPoint = get(handles.axMain, 'CurrentPoint');
    curX = cursorPoint(1,1);
    curY = cursorPoint(1,2);
    
    xLimits = get(handles.axMain, 'xlim');
    yLimits = get(handles.axMain, 'ylim');
    if (curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits))
        disp(['Moving Cursor coordinates are (' num2str(curX) ', ' num2str(curY) ').']);
        guidata(src,handles)
        loadNextImage(src,handles)
        %handles = guidata(src);
    else
        handles.currentImageNumber = handles.currentImageNumber-1;
        disp('Cursor is outside bounds of image.');
    end
    
    pause(0.5);
    handles = guidata(src);
    handles.loadingImageFlag = false;
end
handles.buttonDownMovingFlag = false;
guidata(src,handles);
%}


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.buttonDownFlag = true;
%guidata(hObject,handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
return;
%{
handles.buttonDownFlag = false;
guidata(hObject,handles);
while handles.buttonDownMovingFlag||handles.loadingImageFlag||handles.buttonDownStaticFlag
    handles.buttonDownFlag = false;
    guidata(hObject,handles);
end
%}


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
return;

function loadNextImage(axMain,imageData,trackingData,currentImageNumber,zoomAx)
userData = get(gcf,'userData');
handles = guidata(gcf);
if ~isempty(imageData)&&currentImageNumber>0&&currentImageNumber<=imageData.Frame_Num
    
    I = imread(imageData.Frame_name{currentImageNumber});
    if size(I,3)==1
        I = repmat(I,1,1,3);
    end
    p = prctile(double(I(:)),userData.prc);
    imType = class(I);
    I = double((min(max(I,p(1)),p(2))-p(1)))./(p(2)-p(1))*double(intmax(imType));
    I = cast(I,imType);
    
    %% insert text ID
    
    text_loc = squeeze(userData.cellData(currentImageNumber,:,:));
    
    if size(userData.cellData,3)==1
        text_loc = text_loc';
    end
    [text_locUnique,uniqueCells] = unique(text_loc','rows');
    text_locUnique = text_locUnique';
    
    validCells = uniqueCells(~isnan(text_locUnique(1,:))) ;% ~isnan(text_loc(1,:));
    cellStringAll = userData.strID;
    cellNum = cellStringAll(validCells);
    text_loc = text_loc(:,validCells);
    text_loc = text_loc';
    
    %%

    %%
    if ~isempty(text_loc)&&get(handles.chkShowIDs,'value')
        text_loc(:,2) = text_loc(:,2)-5;
        I = insertText(I,int32(text_loc),cellNum,'FontSize',20,'TextColor','red','BoxOpacity',0);
    end
    lastIdx = find(~isnan(trackingData(:,1)),1,'last');
    fistIdx = find(~isnan(trackingData(:,1)),1,'first');
    if ~isempty(lastIdx)&&get(handles.chkShowIDs,'value')&&lastIdx<currentImageNumber&&currentImageNumber>=fistIdx
        I = insertText(I,int32(trackingData(lastIdx,:)-[0,5]),userData.strID(get(handles.pmCellID,'value')),'FontSize',20,'TextColor','green','BoxOpacity',0);
    end
    %%
    him = imshow(I,p,'parent',axMain);
    hold on;
    scatter(axMain,userData.cellData(currentImageNumber,1,:),userData.cellData(currentImageNumber,2,:),'*r','linewidth',2);
    if ~isempty(trackingData)&&~(all(isnan(trackingData(:))))&&(currentImageNumber>1)
        startid = find(~isnan(trackingData(:,1)),1,'first');
        if startid<=currentImageNumber
        plot(axMain,trackingData(startid:currentImageNumber-1,1),trackingData(startid:currentImageNumber-1,2),'-*','linewidth',2);
        
        
        
        if ~isnan(trackingData(currentImageNumber,1))
            plot(axMain,trackingData(currentImageNumber-1:currentImageNumber,1),trackingData(currentImageNumber-1:currentImageNumber,2),'-*r','linewidth',2);
        end
        scatter(trackingData(startid,1),trackingData(startid,2),'*g','linewidth',2);
        end
    end
    
    
    
    %%
    
    
    if ~isempty(zoomAx)
        set(axMain,'xlim',zoomAx.xlim,'ylim',zoomAx.ylim)
    end
    %set(him,'hittest','off');
    drawnow
    hold off;
end


function TrackCell(mouseStat,varargin)
% Michael Hanchak - Dayton, Ohio, USA
% simpiflied version of "view3d.m" found on Matlab web site.
% only has xy rotation and zoom
% left mouse - rotate
% right mouse - zoom
% double click left - center

% "rot3d off" to turn off
%handles = guidata(gcf);

if nargin<1
    set(gcf,'WindowButtonDownFcn',{@TrackCell,'down'});
    set(gcf,'WindowButtonUpFcn',''); %%set(gcf,'WindowButtonUpFcn',{@TrackCell,'up'});
    set(gcf,'WindowButtonMotionFcn','');
    return;
    %axis vis3d
    %   view all
elseif numel(varargin)==2
    mouseStat = varargin{2};
    handles = guidata(gcf);
    switch mouseStat
        case 'down'
            
            cursorPoint = get(handles.axMain, 'CurrentPoint');
            curX = cursorPoint(1,1);
            curY = cursorPoint(1,2);
            
            xLimits = get(handles.axMain, 'xlim');
            yLimits = get(handles.axMain, 'ylim');
            if ~((curX > min(xLimits) && curX < max(xLimits) && curY > min(yLimits) && curY < max(yLimits)))
                disp('Cursor is outside bounds of image.');
                return;
            end
            userData = get(gcf,'userdata');
            if strcmp(get(gcf,'SelectionType'),'normal')
                
                if userData.isRunning
                    userData.isRunning = false;
                    set(gcf,'userdata',userData);
                else
                    %currentGlobalPos = get(0,'PointerLocation');
                    currentImageNumber = get(handles.pmFrameNum,'value');
                    zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
                    %currentPos = abs(currentGlobalPos-handles.axZeroPoint(1:2)).*([zoomAx.xlim(2)-zoomAx.xlim(1)+1,zoomAx.ylim(2)-zoomAx.ylim(1)+1])./handles.axZeroPoint(3:4)+[round(zoomAx.xlim(1))-1,round(zoomAx.ylim(1))-1];
                    cellId = get(handles.pmCellID,'value');
                    currentPos = [curX,curY];
                    userData.cellData(currentImageNumber,:,cellId) = currentPos;
                    set(gcf,'userdata',userData)
                    if currentImageNumber<userData.imageData.Frame_Num
                        currentImageNumber = currentImageNumber+1;
                        
                        set(handles.pmFrameNum,'value',currentImageNumber);
                        set(handles.sliderMain,'value',(currentImageNumber-1)*1./(userData.imageData.Frame_Num-1));
                        loadNextImage(handles.axMain,userData.imageData,userData.cellData(:,:,cellId),currentImageNumber,zoomAx);
                        
                        if isempty(userData.segData{get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value')})
                            set(handles.pbDrawClear,'string','Draw');
                        else
                            set(handles.pbDrawClear,'string','Clear');
                        end
                    else
                        loadNextImage(handles.axMain,userData.imageData,userData.cellData(:,:,cellId),currentImageNumber,zoomAx);
                    end
                    
                    disp(currentPos);
                end
                
            end
            %{
            if strcmp(get(gcf,'SelectionType'),'open')
                
                if userData.isRunning
                    userData.isRunning = false;
                    set(gcf,'userdata',userData);
                else
                    userData.isRunning = true;
                    
                    set(gcf,'userdata',userData);
                    cellId = get(handles.pmCellID,'value');
                    currentTrack = userData.cellData(:,:,cellId);
                    while userData.isRunning
                        currentImageNumber = get(handles.pmFrameNum,'value');
                        
                        cursorPoint = get(handles.axMain, 'CurrentPoint');
                        curX = cursorPoint(1,1);
                        curY = cursorPoint(1,2);
            
                        zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
                        currentPos =[curX,curY];
                        currentTrack(currentImageNumber,:) = currentPos;
                        if currentImageNumber<userData.imageData.Frame_Num
                            currentImageNumber = currentImageNumber+1;
                            set(handles.pmFrameNum,'value',currentImageNumber);
                            set(handles.sliderMain,'value',(currentImageNumber-1)*1./(userData.imageData.Frame_Num-1));
                            loadNextImage(handles.axMain,userData.imageData,currentTrack,currentImageNumber,zoomAx);
                            
                            if isempty(userData.segData{get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value')})
                                set(handles.pbDrawClear,'string','Draw');
                            else
                                set(handles.pbDrawClear,'string','Clear');
                            end
                            disp(currentPos);
                            pause(0.8);
                        end
                        userData = get(gcf,'userdata');
                    end
                    userData.isRunning = false;
                    userData.cellData(:,:,cellId) = currentTrack;
                    set(gcf,'userdata',userData)
                    
                end
            end
            %}
        case 'up'
            return
            rdata.buttonDown = false;
            %fprintf('!!!!!!!!Button Up Detected!!!!!!!!!!!!!!')
            for i = 1:5
                set(handles.axMain,'userdata',rdata);
                pause(0.001); fprintf('!!!!!!!!Button Up Detected!!!!!!!!!!!!!!')
            end
        case 'off'
            set(gcf,'WindowButtonDownFcn','');
            set(gcf,'WindowButtonUpFcn','');
            set(gcf,'WindowButtonMotionFcn','');
    end
    
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over btnNewCell.
function btnNewCell_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btnNewCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'units','pixels');
set(handles.axMain,'units','pixels');
figPos = get(gcf,'position');
mainAxPos = get(handles.axMain,'position');
handles.axZeroPoint = [figPos(1:2),0 0]+mainAxPos+[0,mainAxPos(4),0 0];
set(gcf,'units','normalized');
set(handles.axMain,'units','normalized');
set(handles.sliderMain,'units','normalized');
guidata(hObject,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1,'userData');
xLim = xlim(handles.axMain);
yLim = ylim(handles.axMain);
if ~isfield(userData,'imageData')
    return;
end
Extr = [userData.imageData.Height,userData.imageData.Width];

if ~ismac&&strcmpi(eventdata.Modifier{1},'control');
    eventdata.Modifier{1} = 'command';
end

if ~isempty(eventdata.Modifier)&&numel(eventdata.Modifier)
    switch eventdata.Modifier{1}
        case 'command'
            switch eventdata.Key
                
                case 'c'
                    if get(handles.tglSegment,'value')&&strcmpi(get(handles.pbDrawClear,'string'),'Clear')
                        pbDrawClear_Callback(handles.pbDrawClear,[],handles)
                    elseif ~get(handles.tglSegment,'value')
                        
                        button = questdlg(sprintf('How much would you like to delete?'),'Delete Cell Track','All','Future','Cancle','Cancle');
                        switch button
                            case 'All'
                                DeleteCellTrack(get(handles.pmCellID,'value'),1,true);
                                %{
                                button2 = questdlg({sprintf('Preparing to delete full track of cell number %s'...
                                    ,userData.strID{get(handles.pmCellID,'value')}),...
                                    'Would you also like to delete tracks of any daughter cell?'},'Delete Daughter Cells?','Yes','No','Cancle','Cancle');
                                switch  button2
                                    case 'Yes'
                                        DeleteCellTrack(get(handles.pmCellID,'value'),1,true);
                                    case 'No'
                                        DeleteCellTrack(get(handles.pmCellID,'value'),1,false);
                                end
                                   %}  
                            case 'Future'
                                DeleteCellTrack(get(handles.pmCellID,'value'),get(handles.pmFrameNum,'value'),true);
                                %{
                                
                                button2 = questdlg({sprintf('Preparing to delete full track of cell number %s'...
                                    ,userData.strID{get(handles.pmCellID,'value')}),...
                                    'Would you also like to delete tracks of any daughter cell?'},'Delete Daughter Cells?','Yes','No','Cancle','Cancle');
                                switch  button2
                                    case 'Yes'
                                        DeleteCellTrack(get(handles.pmCellID,'value'),get(handles.pmFrameNum,'value'),true);
                                    case 'No'
                                        DeleteCellTrack(get(handles.pmCellID,'value'),get(handles.pmFrameNum,'value'),false);
                                end
                                %}
                        end
                               
                    end
                    
                case 'd'
                    if get(handles.tglSegment,'value')&&strcmpi(get(handles.pbDrawClear,'string'),'Draw')
                        pbDrawClear_Callback(handles.pbDrawClear,[],handles)
                    end
                    
                case 'm'
                    if ~get(handles.tglSegment,'value')
                        btnMitosis_Callback(handles.btnMistosis, [], handles);
                    end
                case 'n'
                    if ~get(handles.tglSegment,'value')
                        btnNewCell_Callback(handles.btnNewCell, [], handles);
                    end
                case 's'
                    set(handles.tglSegment,'value',~get(handles.tglSegment,'value'));
                    tglSegment_Callback(handles.tglSegment,[],handles);
                case 'v'
                    set(handles.chkShowIDs,'value',~get(handles.chkShowIDs,'value'));
                    chkShowIDs_Callback(handles.chkShowIDs,[],handles);
                case 'uparrow'
                    zoom(handles.axMain,1.1);
                case 'downarrow'
                    zoom(handles.axMain,1/1.1);
                    
            end
        case 'shift'
            switch eventdata.Key
                case 'rightarrow'
                    framePos = get(handles.pmFrameNum,'value');
                    if framePos<userData.imageData.Frame_Num
                        framePos = framePos + 1;
                        set(handles.pmFrameNum,'value',framePos);
                        step = get(handles.sliderMain,'sliderstep');
                        set(handles.sliderMain,'value',min(get(handles.sliderMain,'value')+step(1),get(handles.sliderMain,'max')))
                        if ~get(handles.tglSegment,'value')
                            trackingData = userData.cellData(:,:,get(handles.pmCellID,'value'));
                            zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
                            loadNextImage(handles.axMain,userData.imageData,trackingData,framePos,zoomAx);
                        else
                            plotSegImage(handles,userData,framePos,get(handles.pmCellID,'value'))
                        end
                        if isempty(userData.segData{get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value')})
                            set(handles.pbDrawClear,'string','Draw');
                        else
                            set(handles.pbDrawClear,'string','Clear');
                        end
                        
                    end
                    
                    
                case 'leftarrow'
                    framePos = get(handles.pmFrameNum,'value');
                    if framePos>1
                        framePos = framePos - 1;
                        set(handles.pmFrameNum,'value',framePos);
                        step = get(handles.sliderMain,'sliderstep');
                        set(handles.sliderMain,'value',max(get(handles.sliderMain,'value')-step(1),get(handles.sliderMain,'min')));
                        if ~get(handles.tglSegment,'value')
                            trackingData = userData.cellData(:,:,get(handles.pmCellID,'value'));
                            zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
                            loadNextImage(handles.axMain,userData.imageData,trackingData,framePos,zoomAx);
                        else
                            plotSegImage(handles,userData,framePos,get(handles.pmCellID,'value'))
                        end
                        if isempty(userData.segData{get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value')})
                            set(handles.pbDrawClear,'string','Draw');
                        else
                            set(handles.pbDrawClear,'string','Clear');
                        end
                    end
                case 'uparrow'
                    cellID = get(handles.pmCellID,'value');
                    cellIDStr = get(handles.pmCellID,'String');
                    if get(handles.tglSegment,'value')
                        currentFrame = get(handles.pmFrameNum,'value');
                        
                        validCellID = find(~isnan(userData.cellData(currentFrame,1,:)));
                        nextCellID  = validCellID(find(validCellID>cellID,1,'first'));
                        if isempty(nextCellID)
                            nextCellID = cellID;
                        end
                    else
                        nextCellID = cellID+1;
                    end
                    if nextCellID<=numel(cellIDStr)
                        set(handles.pmCellID,'value',nextCellID);
                        pmCellID_Callback(handles.pmCellID,[],handles)
                    end
                case 'downarrow'
                    cellID = get(handles.pmCellID,'value');
                    if get(handles.tglSegment,'value')
                        currentFrame = get(handles.pmFrameNum,'value');
                        
                        validCellID = find(~isnan(userData.cellData(currentFrame,1,:)));
                        prevCellID  = validCellID(find(validCellID<cellID,1,'last'));
                        if isempty(prevCellID)
                            prevCellID = cellID;
                        end
                    else
                        prevCellID = cellID-1;
                    end
                    if prevCellID>=1
                        set(handles.pmCellID,'value',prevCellID);
                        pmCellID_Callback(handles.pmCellID,[],handles)
                    end
            end
        otherwise
    end
    
%    ctrNew = strcmpi(eventdata.Modifier,'command')&&strcmpi(eventdata.Key,'n');
else
    switch eventdata.Key
        case 'uparrow'
            height = (yLim(2)-yLim(1));
            step = height*0.05;
            yLimNew = yLim - step;
            if yLimNew(1)<0
                yLimNew = [0 height];
            end
            ylim(handles.axMain,yLimNew);
            
            
        case 'downarrow'
            height = (yLim(2)-yLim(1));
            step = height*0.05;
            yLimNew = yLim + step;
            if yLimNew(2) > Extr(1)
                yLimNew = [Extr(1)-height,Extr(1)];
            end
            ylim(handles.axMain,yLimNew);
            
        case 'leftarrow'
            width = (xLim(2)-xLim(1));
            step = width*0.05;
            xLimNew = xLim - step;
            if xLimNew(1)<0
                xLimNew = [0 width];
            end
            xlim(handles.axMain,xLimNew);
        case 'rightarrow'
            width = (xLim(2)-xLim(1));
            step = width*0.05;
            xLimNew = xLim + step;
            if xLimNew(2)>Extr(2)
                xLimNew = [Extr(2)-width,Extr(2)];
            end
            xlim(handles.axMain,xLimNew);
    end
end





% --- Executes on button press in btnMitosis.
function btnMitosis_Callback(hObject, eventdata, handles)
% hObject    handle to btnMitosis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(gcf,'userData');
currentCell = get(handles.pmCellID,'Value');
cellStrings =  userData.strID;
cellStringsColor =  userData.strColorID;
if iscell(cellStrings)
    currentCellString = cellStrings{currentCell};
else
    currentCellString = cellStrings;
    cellStrings = {cellStrings};
end
[x,y] = ginput;
newString = cell(numel(x),1);
tokensAll = regexp(cellStrings,'(?<base>\d+).(?<gen>\w)(?<genid>\d+)','tokens');
valid = cellfun(@(t) ~isempty(t),tokensAll);
base = cellfun(@(t) t{1}{1},tokensAll(valid),'uniformoutput',false);
gen = cellfun(@(t) t{1}{2},tokensAll(valid),'uniformoutput',false);
genID = cellfun(@(t) t{1}{3},tokensAll(valid),'uniformoutput',false);
for i = 1:numel(x)
    if isempty(strfind(currentCellString,'.'))
        newString{i} = sprintf('%s.A%d',currentCellString,i);
    else
        t = regexp(currentCellString,'(?<base>\d+).(?<gen>\w)(?<genid>\d+)','tokens');
        t = t{1};
        
        if t{2}=='Z'
            error('Too many generations for cell #%s',t{1})
        end
        newGen = char(t{2}+1);
        prevGenID = str2double(cat(1,genID(strcmpi(t{1},base)&strcmpi(newGen,gen))));
        if isempty(prevGenID)
            newGenID = i;
        else
            newGenID = max(prevGenID)+i;
        end
        newString{i} = sprintf('%s.%s%d',t{1},newGen,newGenID);
            
    %newString{i} = sprintf('%s.%d',currentCellString,i);
    end
    userData.cellData(:,:,end+1) = nan;%userData.cellData(:,:,currentCell);
    [userData.segData{:,end+1}] = deal([]);
    [userData.contourData{:,end+1}] = deal([]);
    currentImageNumber = get(handles.pmFrameNum,'value');
    userData.cellData(currentImageNumber,:,end) = [x(i),y(i)];
    userData.Link(end+1).Mother = currentCell;
    
end
userData.Link(end).Children = (size(userData.cellData,3)-numel(x)+1):(size(userData.cellData,3));
userData.Link(end).Time = get(handles.pmFrameNum,'value');
cellStrings = [cellStrings;newString];
if numel(x)>1
    newString(2:end) = cellfun(@(ns) sprintf('<HTML><FONT COLOR="red">%s</HTML>',ns),newString(2:end),'uniformoutput',false);
    
end
cellStringsColor = [cellStringsColor;newString];
userData.strID = cellStrings;
userData.strColorID = cellStringsColor;
set(gcf,'userdata',userData);
cellId = (numel(cellStrings)-numel(x)+1);
set(handles.pmCellID,'String',cellStringsColor,'value',cellId);
zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
currentImageNumber = currentImageNumber+1;

if numel(x)>1
    if ~isfield(userData,'returnTo')
        userData.returnTo =numel(cellStrings)-numel(x)+2:numel(cellStrings);
    else
        userData.returnTo = [userData.returnTo,numel(cellStrings)-numel(x)+2:numel(cellStrings)];
    end
end
set(handles.pmFrameNum,'value',currentImageNumber);
set(handles.sliderMain,'value',(currentImageNumber-1)*1./(userData.imageData.Frame_Num-1));
loadNextImage(handles.axMain,userData.imageData,userData.cellData(:,:,cellId),currentImageNumber,zoomAx);

if isempty(userData.segData{get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value')})
    set(handles.pbDrawClear,'string','Draw');
else
    set(handles.pbDrawClear,'string','Clear');
end
set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');


% --- Executes on button press in tglSegment.
function tglSegment_Callback(hObject, eventdata, handles)
% hObject    handle to tglSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglSegment
set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');
tglState = get(handles.tglSegment,'value');
userData = get(gcf,'userdata');
currentImageNumber = get(handles.pmFrameNum,'value');
cellID = get(handles.pmCellID,'value');
if tglState
    TrackCell([],[],'off');
    set(handles.pbDrawClear,'visible','on');
    plotSegImage(handles,userData,currentImageNumber,cellID);
    if isempty(userData.segData{get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value')})
        set(handles.pbDrawClear,'string','Draw');
    else
        set(handles.pbDrawClear,'string','Clear');
    end
    %{
    currentImageNumber = get(handles.pmFrameNum,'value');
    I = imread(userData.imageData.Frame_name{currentImageNumber});
    p = prctile(double(I(:)),[0.1,99.9]);
    
    centAll = squeeze(userData.cellData(currentImageNumber,:,:));
    trackedCells = ~isnan(centAll(1,:));
    allSeg = userData.segData(currentImageNumber,:);
    emptySeg = cellfun(@(seg) isempty(seg),allSeg);
    toSeg = emptySeg&&trackedCells;
    shapeInserter = vision.ShapeInserter('shape','Polygons','BorderColor','custom','CustomBorderColor',[1,0,0]);
    
    while tglState && any(toSeg)
        cellID = find(toSeg,1,'first');
        Ishow = I;
        for i = 1:numel(userData.segData(currentImageNumber,:))
            Ishow = step(shapeInserter,Ishow,userData.segData{currentImageNumber,i});
        end
        him = imshow(Ishow,p,'parent',handles.axMain);
        hold on;
        scatter(handles.axMain,userData.cellData(currentImageNumber,1,:),userData.cellData(currentImageNumber,2,:),'*r','linewidth',2);
        scatter(handles.axMain,userData.cellData(currentImageNumber,1,cellID),userData.cellData(currentImageNumber,2,cellID),'*g','linewidth',2);
        
        %% Plot Available Contours
        hold off;
        %% Add Segmenation Part
        DrawCell(handles,currentImageNumber,cellID)
        %%
        toSeg(currentImageNumber) = false;
        tglState = get(hObject,'value');
    end
    if ~any(toSeg)&&tglState
        set(hObject,'value',0);
    end
    %}
else
    set(handles.pbDrawClear,'visible','off');
    trackingData = userData.cellData(:,:,cellID);
    zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
    loadNextImage(handles.axMain,userData.imageData,trackingData,currentImageNumber,zoomAx)
    
    if isempty(userData.segData{get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value')})
        set(handles.pbDrawClear,'string','Draw');
    else
        set(handles.pbDrawClear,'string','Clear');
    end
    TrackCell;
end

function DrawCell(handles,currentImageNumber,cellID)

userData = get(gcf,'userdata');
freezeGUI('off');
cellContour = imfreehand(handles.axMain);
freezeGUI('on');
Mask = cellContour.createMask;
Contour = bwperim(Mask);
userData.segData{currentImageNumber,cellID} = find(Mask(:));
userData.contourData{currentImageNumber,cellID} = find(Contour(:));
cellContour.delete;
set(gcf,'userdata',userData);

function DeleteCell(currentImageNumber,cellID)
userData = get(gcf,'userdata');
userData.segData(currentImageNumber,cellID) = {[]};
userData.contourData(currentImageNumber,cellID) = {[]};
set(gcf,'userdata',userData);




% --- Executes on button press in pbDrawClear.
function pbDrawClear_Callback(hObject, eventdata, handles)
% hObject    handle to pbDrawClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');
userData = get(gcf,'userdata');
currentImageNumber = get(handles.pmFrameNum,'value');
cellID = get(handles.pmCellID,'value');
switch get(hObject,'String')
    case 'Draw'
        DrawCell(handles,currentImageNumber,cellID)
        set(hObject,'String','Clear')
    case 'Clear'
        DeleteCell(currentImageNumber,cellID)
        set(hObject,'String','Draw')
end
userData = get(gcf,'userdata');
plotSegImage(handles,userData,currentImageNumber,cellID);

function plotSegImage(handles,userData,currentImageNumber,cellID)
zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
I = imread(userData.imageData.Frame_name{currentImageNumber});
p = prctile(double(I(:)),userData.prc);
imType = class(I);
I = double((min(max(I,p(1)),p(2))-p(1)))./(p(2)-p(1))*double(intmax(imType));
I = cast(I,imType);

if size(I,3)==1
    R = I;
    G = I;
    B = I;
else
    R = I(:,:,1);
    G = I(:,:,2);
    B = I(:,:,3);
end

ContoursIdx = cat(1,userData.contourData{currentImageNumber,:});
ContoursLastIdx = userData.contourData{currentImageNumber,cellID};
ContoursIdx = setdiff(ContoursIdx,ContoursLastIdx);
R(ContoursIdx) = intmax(class(R));
G(ContoursIdx) = 0;
B([ContoursIdx(:);ContoursLastIdx(:)]) = 0;
G(ContoursLastIdx) = intmax(class(R));
R(ContoursLastIdx) = 0;
IShow = cat(3,R,G,B);
text_loc = squeeze(userData.cellData(currentImageNumber,:,:));
    if size(userData.cellData,3)==1
        text_loc = text_loc';
    end
    [text_locUnique,uniqueCells] = unique(text_loc','rows');
    text_locUnique = text_locUnique';
    
validCells = uniqueCells(~isnan(text_locUnique(1,:))) ;% ~isnan(text_loc(1,:));
cellStringAll = userData.strID;

if any(validCells==cellID)
    showCell = true;
    validCells = setdiff(validCells,cellID);
    showCellLoc = text_loc(:,cellID);
else
    showCell = false;
    
end
cellNum = cellStringAll(validCells);
text_loc = text_loc(:,validCells);
text_loc = text_loc';


if ~isempty(text_loc)&&get(handles.chkShowIDs,'value')
    text_loc(:,2) = text_loc(:,2)-5;
    IShow = insertText(IShow,int32(text_loc),cellNum,'FontSize',20,'TextColor','red','BoxOpacity',0);
    if showCell
        IShow = insertText(IShow,int32(showCellLoc'-[0,5]),cellStringAll(cellID),'FontSize',20,'TextColor','green','BoxOpacity',0);
    end
end
imshow(IShow,'parent',handles.axMain);
set(handles.axMain,'xlim',zoomAx.xlim,'ylim',zoomAx.ylim)
hold on;
scatter(userData.cellData(currentImageNumber,1,:),userData.cellData(currentImageNumber,2,:),'*r','linewidth',2);
scatter(userData.cellData(currentImageNumber,1,cellID),userData.cellData(currentImageNumber,2,cellID),'*g','linewidth',2);
hold off;

function freezeGUI(status)

handles = guidata(gcf);
set(handles.btnCreateCSV,'enable',status);
set(handles.btnLoad,'enable',status);
set(handles.btnMitosis,'enable',status);
set(handles.btnNew,'enable',status);
set(handles.btnPlayPause,'enable',status);
set(handles.btnSave,'enable',status);
set(handles.pbDrawClear,'enable',status);
set(handles.pmCellID,'enable',status);
set(handles.pmFrameNum,'enable',status);
set(handles.sliderMain,'enable',status);
set(handles.btnPlayPause,'enable',status);
set(handles.tglSegment,'enable',status);
set(handles.btnNewCell,'enable',status);
handles.Frozen = status;
guidata(gcf,handles);

function DeleteCellTrack(cellID,framePos,Recursive)

userData = get(gcf,'userData');
userData.cellData(framePos:end,:,cellID) = nan;
if all(isnan(userData.cellData(:,1,cellID)))&&size(userData.cellData(1)>1)
[userData.segData{framePos:end,cellID}] = deal([]);
[userData.contourData{framePos:end,cellID}] = deal([]);
else
    %{
    userData.segData(framePos:end,cellID) = [];
    userData.contourData(framePos:end,cellID) = [];
    userData.cellData(framePos:end,:,cellID) = [];
    userData.strID{cellID} = [];
    userData.strColorID{cellID} = [];
    %}
end
set(gcf,'userdata',userData);
%if Recursive
%end
%set(handles.pmCellID,'string', userData.strColorID,'value',max(numel(userData.strColorID),cellID));


% --- Executes on button press in chkShowIDs.
function chkShowIDs_Callback(hObject, eventdata, handles)
% hObject    handle to chkShowIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkShowIDs
if isempty(get(gcf,'userData'))
    return;
end
if get(handles.tglSegment,'value');
    plotSegImage(handles,get(gcf,'userData'),get(handles.pmFrameNum,'value'),get(handles.pmCellID,'value'))
else
    userData = get(gcf,'userData');
    trackingData = userData.cellData(:,:,get(handles.pmCellID,'value'));
    zoomAx.xlim = get(handles.axMain,'xLim'); zoomAx.ylim = get(handles.axMain,'yLim');
    loadNextImage(handles.axMain,userData.imageData,trackingData,get(handles.pmFrameNum,'value'),zoomAx);
end
set(hObject,'enable','off');
drawnow update;
set(hObject,'enable','on');