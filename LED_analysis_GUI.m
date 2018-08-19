function varargout = LED_analysis_GUI(varargin)
% LED_ANALYSIS_GUI MATLAB code for LED_analysis_GUI.fig
%      LED_ANALYSIS_GUI, by itself, creates a new LED_ANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = LED_ANALYSIS_GUI returns the handle to a new LED_ANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      LED_ANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LED_ANALYSIS_GUI.M with the given input arguments.
%
%      LED_ANALYSIS_GUI('Property','Value',...) creates a new LED_ANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LED_analysis_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LED_analysis_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LED_analysis_GUI

% Last Modified by GUIDE v2.5 25-Nov-2017 10:53:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LED_analysis_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LED_analysis_GUI_OutputFcn, ...
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


% --- Executes just before LED_analysis_GUI is made visible.
function LED_analysis_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LED_analysis_GUI (see VARARGIN)

% Choose default command line output for LED_analysis_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LED_analysis_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LED_analysis_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in start_analysis.
function start_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to start_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
WL_limit(1,1)=str2double(get(handles.WL_low1,'String'));
WL_limit(1,2)=str2double(get(handles.WL_high1,'String'));
Vth_def= str2double(get(handles.Vth_def1,'String'));
data_folder=get(handles.Address1,'String');
[all_analysis]=LED_all_analysis8(data_folder,WL_limit,Vth_def);



function Address1_Callback(hObject, eventdata, handles)
% hObject    handle to Address1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Address1 as text
%        str2double(get(hObject,'String')) returns contents of Address1 as a double
data_folder=get(handles.Address1,'String');

% --- Executes during object creation, after setting all properties.
function Address1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Address1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vth_def1_Callback(hObject, eventdata, handles)
% hObject    handle to Vth_def1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vth_def1 as text
%        str2double(get(hObject,'String')) returns contents of Vth_def1 as a double
Vth_def= str2double(get(handles.Vth_def1,'String'));

% --- Executes during object creation, after setting all properties.
function Vth_def1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vth_def1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on start_analysis and none of its controls.
function start_analysis_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to start_analysis (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)




function WL_low1_Callback(hObject, eventdata, handles)
% hObject    handle to WL_low1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WL_low1 as text
%        str2double(get(hObject,'String')) returns contents of WL_low1 as a double
WL_limit(1,2)=str2double(get(handles.WL_low1,'String'));

% --- Executes during object creation, after setting all properties.
function WL_low1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WL_low1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WL_high1_Callback(hObject, eventdata, handles)
% hObject    handle to WL_high1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WL_high1 as text
%        str2double(get(hObject,'String')) returns contents of WL_high1 as a double
WL_limit(1,1)=str2double(get(handles.WL_high1,'String'));

% --- Executes during object creation, after setting all properties.
function WL_high1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WL_high1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
