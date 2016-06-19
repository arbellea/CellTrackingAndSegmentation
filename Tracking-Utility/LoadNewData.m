function varargout = LoadNewData(varargin)
% LOADNEWDATA MATLAB code for LoadNewData.fig
%      LOADNEWDATA, by itself, creates a new LOADNEWDATA or raises the existing
%      singleton*.
%
%      H = LOADNEWDATA returns the handle to a new LOADNEWDATA or the handle to
%      the existing singleton*.
%
%      LOADNEWDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADNEWDATA.M with the given input arguments.
%
%      LOADNEWDATA('Property','Value',...) creates a new LOADNEWDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LoadNewData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LoadNewData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LoadNewData

% Last Modified by GUIDE v2.5 17-Jan-2016 14:05:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoadNewData_OpeningFcn, ...
                   'gui_OutputFcn',  @LoadNewData_OutputFcn, ...
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


% --- Executes just before LoadNewData is made visible.
function LoadNewData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LoadNewData (see VARARGIN)

% Choose default command line output for LoadNewData
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LoadNewData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LoadNewData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function txtRawDir_Callback(hObject, eventdata, handles)
% hObject    handle to txtRawDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtRawDir as text
%        str2double(get(hObject,'String')) returns contents of txtRawDir as a double


% --- Executes during object creation, after setting all properties.
function txtRawDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtRawDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pbBrowseRaw.
function pbBrowseRaw_Callback(hObject, eventdata, handles)
% hObject    handle to pbBrowseRaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1,'userData');
userData.rawDir = uigetdir;
set(handles.figure1,'userData',userData);


function txtRawExp_Callback(hObject, eventdata, handles)
% hObject    handle to txtRawExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtRawExp as text
%        str2double(get(hObject,'String')) returns contents of txtRawExp as a double


% --- Executes during object creation, after setting all properties.
function txtRawExp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtRawExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function guessRegExp(hObject);

userdata = get(hObject,'userdata');
fileInfo ;
