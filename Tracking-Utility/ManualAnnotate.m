function varargout = ManualAnnotate(varargin)
% MANUALANNOTATE MATLAB code for ManualAnnotate.fig
%      MANUALANNOTATE, by itself, creates a new MANUALANNOTATE or raises the existing
%      singleton*.
%
%      H = MANUALANNOTATE returns the handle to a new MANUALANNOTATE or the handle to
%      the existing singleton*.
%
%      MANUALANNOTATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUALANNOTATE.M with the given input arguments.
%
%      MANUALANNOTATE('Property','Value',...) creates a new MANUALANNOTATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ManualAnnotate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ManualAnnotate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ManualAnnotate

% Last Modified by GUIDE v2.5 08-Jun-2016 15:45:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ManualAnnotate_OpeningFcn, ...
    'gui_OutputFcn',  @ManualAnnotate_OutputFcn, ...
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


% --- Executes just before ManualAnnotate is made visible.
function ManualAnnotate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ManualAnnotate (see VARARGIN)

% Choose default command line output for ManualAnnotate
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
linkaxes([handles.axes1,handles.axes2]);
% UIWAIT makes ManualAnnotate wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ManualAnnotate_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pbNewCell.
function pbNewCell_Callback(hObject, eventdata, handles)
% hObject    handle to pbNewCell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellString = get(handles.popCellList,'string');
if ischar(cellString)
    cellString = {1};
else
    newCell = str2double(cellString{end})+1;
    cellString{end+1} = newCell;
end
set(handles.popCellList,'string',cellString);
set(handles.popCellList,'value',numel(cellString));
DrawCell(1,handles);
DrawCell(2,handles);







% --- Executes on button press in pbClear2.
function pbClear2_Callback(hObject, eventdata, handles)
% hObject    handle to pbClear2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedCell = get(handles.popCellList,'value');
userData = get(handles.figure1,'userdata');
userData.cellsContours{selectedCell,2} = [];
userData.cellsCents{selectedCell,2} = [];
set(handles.figure1,'userdata',userData);
createContourMask(handles.figure1);


% --- Executes on button press in pbClear1.
function pbClear1_Callback(hObject, eventdata, handles)
% hObject    handle to pbClear1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selectedCell = get(handles.popCellList,'value');
userData = get(handles.figure1,'userdata');
userData.cellsContours{selectedCell,1} = [];
userData.cellsCents{selectedCell,1} = [];
set(handles.figure1,'userdata',userData);
createContourMask(handles.figure1);

% --- Executes on selection change in popCellList.
function popCellList_Callback(hObject, eventdata, handles)
% hObject    handle to popCellList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popCellList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popCellList


% --- Executes during object creation, after setting all properties.
function popCellList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popCellList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbLoad.
function pbLoad_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button = questdlg('Would you like to continue an existing project or start a new project?','Continue?','New','Continue','Cancel');
switch button
    case 'New'
        [fileName1,pathName1] = uigetfile('*.*');
        if fileName1==0
            return;
        end
        [fileName2,pathName2] = uigetfile(fullfile(pathName1,'*.*'));
        if fileName2==0
            return;
        end
        set(handles.pbNewCell,'enable','on')
        set(handles.popCellList,'String','');
        I1 = imread(fullfile(pathName1,fileName1));
        I2 = imread(fullfile(pathName2,fileName2));
        imType = class(I1);
        p = prctile(double(cat(1,I1(:),I2(:))),[0.1,99.9]);
        I1 = double((min(max(I1,p(1)),p(2))-p(1)))./(p(2)-p(1))*double(intmax(imType));
        I2 = double((min(max(I2,p(1)),p(2))-p(1)))./(p(2)-p(1))*double(intmax(imType));
        I1 = cast(I1,imType);
        I2 = cast(I2,imType);
        handles.im1 = imshow(I1,[],'parent',handles.axes1);
        %set(handles.im1,'buttonDownFcn',@DrawCell,'hittest','on','tag','image1');
        set(handles.im1,'hittest','on','tag','image1');
        handles.im2 = imshow(I2,[],'parent',handles.axes2);
        %set(handles.im2,'buttonDownFcn',@DrawCell,'hittest','on','tag','image2');
        set(handles.im2,'hittest','on','tag','image2');
        %userData = get(handles.figure1,'UserData');
        userData = [];
        userData.LoadData.pathName1 = pathName1;
        userData.LoadData.pathName2 = pathName2;
        userData.LoadData.fileName1 = fileName1;
        userData.LoadData.fileName2 = fileName2;
        userData.Images.Raw{1} = repmat(I1,1,1,3);
        userData.Images.Raw{2} = repmat(I2,1,1,3);
        userData.Images.Show{1} = repmat(I1,1,1,3);
        userData.Images.Show{2} = repmat(I2,1,1,3);
        userData.Images.contourMasks{1} = false(size(I1));
        userData.Images.contourMasks{2} = false(size(I2));
        userData.Cells = {};
        userData.cellsContours  = {};
        userData.shapeInserter = vision.ShapeInserter('Shape','Polygons','BorderColor','Custom','CustomBorderColor',[intmax(imType),0,0]);
        set(handles.figure1,'UserData',userData)
        guidata(handles.figure1,handles);
    case 'Continue'
        [fileName,pathName] = uigetfile('*.mat');
        try
            load(fullfile(pathName,fileName));
        catch err
            error('File is corrupt');
        end
        if ~exist('userData','var')
            error('File is corrupt');
        end
        set(handles.figure1,'userData',userData); %#ok<NODEF>
        set(handles.pbNewCell,'enable','on');
        
        cellString = num2cell(1:max(1,size(userData.cellsContours,1)));
        set(handles.popCellList,'string',cellString);
        handles.im1 = imshow(userData.Images.Show{1},[],'parent',handles.axes1);
        set(handles.im1,'hittest','on','tag','image1');
        handles.im2 = imshow(userData.Images.Show{2},[],'parent',handles.axes2);
        
    case 'Cancel'
end


% --- Executes on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
% hObject    handle to pbSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1,'UserData');
[FileName,PathName] = uiputfile(fullfile(userData.LoadData.pathName1,'*.mat'),'Save Manual Segmentation','ManualSeg.mat');
save(fullfile(PathName,FileName),'userData');


% --- Executes on mouse press over axes background.
function DrawCell(frameNum,handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~exist('handles','var')||isempty(handles)
    handles = guidata(hObject);
end
userData = get(handles.figure1,'userdata');
cellString = get(handles.popCellList,'string');
selectedCell = get(handles.popCellList,'value');
axes = [handles.axes1,handles.axes2];

if isempty(cellString)
    cellString = '1';
    set(handles.popCellList,'string',cellString);
end



switchFreezeGUI(handles.figure1	,false)
cellContour = imfreehand(axes(frameNum));
switchFreezeGUI(handles.figure1,true)
userData.Cells{selectedCell,frameNum} = cellContour.getPosition;
Mask = bwperim(cellContour.createMask);
userData.cellsContours{selectedCell,frameNum} = find(Mask(:));
userData.cellsCents{selectedCell,frameNum} = mean(cellContour.getPosition);
set(handles.figure1,'userdata',userData);
%mask = bwperim(cellContour.createMask);
%userData.Images.contourMasks{frameNum} = userData.Images.contourMasks{frameNum}&mask;
createContourMask(handles.figure1)
cellContour.delete;

function createContourMask(hObject)

handles = guidata(hObject);
userData = get(handles.figure1,'userdata');
xLim = xlim(handles.axes1);
yLim = ylim(handles.axes1);
for frame = 1:size(userData.Cells,2)
    ContoursIdx = cat(1,userData.cellsContours{:,frame});
    R = userData.Images.Raw{frame}(:,:,1);
    G = userData.Images.Raw{frame}(:,:,2);
    B = userData.Images.Raw{frame}(:,:,3);
    R(ContoursIdx) = intmax(class(R));
    G(ContoursIdx) = 0;
    B(ContoursIdx) = 0;
    userData.Images.Show{frame} = cat(3,R,G,B);
    cellNum = 1:size(userData.cellsContours,1);
    text_loc = userData.cellsCents(:,frame);
    validCells = cellfun(@(cell) ~isempty(cell),userData.cellsCents(:,frame));
    cellNum = cellstr(num2str(cellNum(validCells)'));
    text_loc = text_loc(validCells,:);
    text_loc = vertcat(text_loc{:});
    if ~isempty(text_loc)
        text_loc(:,2) = text_loc(:,2)-5;
        
        userData.Images.Show{frame} = insertText(userData.Images.Show{frame},int32(text_loc),cellNum,'FontSize',20,'TextColor','red','BoxOpacity',0);
    end
    %polygons = cellfun(@(cell) reshape(cell',1,[]),userData.Cells(:,frame),'uniformoutput',false);
    %polygons = uint16(vertcat(polygons{:}));
    %userData.Images.Show{frame} = step(userData.shapeInserter,userData.Images.Raw{frame},polygons);
    
    if frame==1
        
        handles.im1 = imshow(userData.Images.Show{frame},[],'parent',handles.axes1);
        %set(handles.im1,'buttonDownFcn',@DrawCell,'hittest','on','tag','image1');
        set(handles.im1,'hittest','on','tag','image1');
        
        
        
        xlim(handles.axes1,xLim);
        ylim(handles.axes1,yLim);
    else
        handles.im2 = imshow(userData.Images.Show{frame},[],'parent',handles.axes2);
        %set(handles.im2,'buttonDownFcn',@DrawCell,'hittest','on','tag','image2');
        set(handles.im2,'hittest','on','tag','image2');
        xlim(handles.axes2,xLim);
        ylim(handles.axes2,yLim);
    end
end
set(handles.figure1,'userdata',userData);
guidata(handles.figure1,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1,'userData');
xLim = xlim(gca);
yLim = ylim(gca);
Extr = size(userData.Images.Raw{1});
cntrlButton = 'control';
if ismac
   cntrlButton= 'command';
end
    if ~isempty(eventdata.Modifier)&numel(eventdata.Modifier)
        switch lower(eventdata.Modifier{1})
            case cntrlButton
                switch eventdata.Key
                    case 'n'
                        if ~(strcmpi(get(handles.pbNewCell,'Enable'),'off'))
                        pbNewCell_Callback(hObject, eventdata, handles)
                        end
                    case 'c'
                        if ~(strcmpi(get(handles.pbNewCell,'Enable'),'off'))
                        if isequal(gca,handles.axes1)
                            pbClear1_Callback(hObject, eventdata, handles);
                        elseif isequal(gca,handles.axes2)
                            pbClear2_Callback(hObject, eventdata, handles);
                        end
                        end
                    case 'd'
                        if ~(strcmpi(get(handles.pbNewCell,'Enable'),'off'))
                        if isequal(gca,handles.axes1)
                            pbDraw1_Callback(hObject, eventdata, handles);
                        elseif isequal(gca,handles.axes2)
                            pbDraw2_Callback(hObject, eventdata, handles);
                        end
                        end
                    case 'uparrow'
                        zoom(gca,1.1);
                    case 'downarrow'
                        zoom(gca,1/1.1);
                end
                
            otherwise
        end
        
        ctrNew = strcmpi(eventdata.Modifier,'command')&&strcmpi(eventdata.Key,'n');
    else
        switch eventdata.Key
            case 'uparrow'
                height = (yLim(2)-yLim(1));
                step = height*0.05;
                yLimNew = yLim - step;
                if yLimNew(1)<0
                    yLimNew = [0 height];
                end
                ylim(gca,yLimNew);
                
                
            case 'downarrow'
                height = (yLim(2)-yLim(1));
                step = height*0.05;
                yLimNew = yLim + step;
                if yLimNew(2) > Extr(1)
                    yLimNew = [Extr(1)-height,Extr(1)];
                end
                ylim(gca,yLimNew);
                
            case 'leftarrow'
                width = (xLim(2)-xLim(1));
                step = width*0.05;
                xLimNew = xLim - step;
                if xLimNew(1)<0
                    xLimNew = [0 width];
                end
                xlim(gca,xLimNew);
            case 'rightarrow'
                width = (xLim(2)-xLim(1));
                step = width*0.05;
                xLimNew = xLim + step;
                if xLimNew(2)>Extr(2)
                    xLimNew = [Extr(2)-width,Extr(2)];
                end
                xlim(gca,xLimNew);
        end
    end
    

%{
if ctrNew&&strcmpi(get(handles.pbNewCell,'enable'),'on')
    zoom(gca,'off'); pan(gca,'off');
    pbNewCell_Callback(hObject, eventdata, handles)
end
%}


% --- Executes on button press in pbDraw1.
function pbDraw1_Callback(hObject, eventdata, handles)
% hObject    handle to pbDraw1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%selectedCell = get(handles.popCellList,'value');
DrawCell(1,handles);



% --- Executes on button press in pbDraw2.
function pbDraw2_Callback(hObject, eventdata, handles)
% hObject    handle to pbDraw2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%selectedCell = get(handles.popCellList,'value');
DrawCell(2,handles);


% --- Executes on button press in pbDone.
function pbDone_Callback(hObject, eventdata, handles)
% hObject    handle to pbDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1,'userData');
saveDir = uigetdir(userData.LoadData.pathName1);
valid = cellfun(@(cell) ~isempty(cell), userData.cellsContours);
if any(~valid(:))
    warndlg('Missing Segmentations')
    return;
end

for frame = 1:2
    Seg = zeros(size(userData.Images.Raw{1},1),size(userData.Images.Raw{1},2));
    for cell = 1:size(userData.cellsContours(:,frame))
        
        Mask = poly2mask(userData.Cells{cell,frame}(:,1),userData.Cells{cell,frame}(:,2),size(userData.Images.Raw{1},1),size(userData.Images.Raw{1},2));
        Seg(Mask) = cell;
    end
    fileName = sprintf('ManualSeg_%d.tif',frame);
    imwrite(uint16(Seg),fullfile(saveDir,fileName))
end


% --- Executes on button press in pbCenter.
function pbCenter_Callback(hObject, eventdata, handles)
% hObject    handle to pbCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1,'userData');
M = poly2mask(userData.Cells{1,1}(:,1),userData.Cells{1,1}(:,2),512,641);


function switchFreezeGUI(hObject,flag)
handles = guidata(hObject);
if ~flag
    set(handles.pbClear1,'Enable','off');
    set(handles.pbClear2,'Enable','off');
    set(handles.pbCenter,'Enable','off');
    set(handles.pbDraw1,'Enable','off');
    set(handles.pbDraw2,'Enable','off');
    set(handles.pbLoad,'Enable','off');
    set(handles.pbNewCell,'Enable','off');
    set(handles.pbSave,'Enable','off');
    set(handles.popCellList,'Enable','off');
else
    set(handles.pbClear1,'Enable','on');
    set(handles.pbClear2,'Enable','on');
    set(handles.pbCenter,'Enable','on');
    set(handles.pbDraw1,'Enable','on');
    set(handles.pbDraw2,'Enable','on');
    set(handles.pbLoad,'Enable','on');
    set(handles.pbNewCell,'Enable','on');
    set(handles.pbSave,'Enable','on');
    set(handles.popCellList,'Enable','on');
end


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
if eventdata.VerticalScrollCount<0
    zoom(gca,1.1.^abs(eventdata.VerticalScrollCount))
else
    zoom(gca,1./(1.1).^abs(eventdata.VerticalScrollCount))
end


           


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if any(strcmpi(get(gco,'Type'),{'axes','image'}))
    %pan(gca,'on');
end


% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if any(strcmpi(get(gco,'Type'),{'axes','image'}))
  %  pan(gca,'off');
end