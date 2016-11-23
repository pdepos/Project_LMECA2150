function varargout = interfaceTurbine(varargin)
% INTERFACETURBINE MATLAB code for interfaceTurbine.fig
%      INTERFACETURBINE, by itself, creates a new INTERFACETURBINE or raises the existing
%      singleton*.
%
%      H = INTERFACETURBINE returns the handle to a new INTERFACETURBINE or the handle to
%      the existing singleton*.
%
%      INTERFACETURBINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFACETURBINE.M with the given input arguments.
%
%      INTERFACETURBINE('Property','Value',...) creates a new INTERFACETURBINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before interfaceTurbine_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to interfaceTurbine_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help interfaceTurbine

% Last Modified by GUIDE v2.5 20-Nov-2016 15:13:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @interfaceTurbine_OpeningFcn, ...
                   'gui_OutputFcn',  @interfaceTurbine_OutputFcn, ...
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


% --- Executes just before interfaceTurbine is made visible.
function interfaceTurbine_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to interfaceTurbine (see VARARGIN)

% Choose default command line output for interfaceTurbine
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
initialize_gui(hObject, handles);


% UIWAIT makes interfaceTurbine wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = interfaceTurbine_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in recuperator.
function recuperator_Callback(hObject, eventdata, handles)
% hObject    handle to recuperator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recuperator
if get(hObject,'Value')==1
    % change cycle figure
    axes(handles.cycle);
    img = imread('TurbineGasRecup.PNG');
    image(img);
    set(gca,'Visible','off');
    
    
else
    initialize_gui(hObject, handles);
end


function T1v_Callback(hObject, eventdata, handles)
% hObject    handle to T1v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T1v as text
%        str2double(get(hObject,'String')) returns contents of T1v as a double
T1 = str2double(get(hObject, 'String'));
handles.metricdata.T1 = T1;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function T1v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T1v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rv_Callback(hObject, eventdata, handles)
% hObject    handle to rv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rv as text
%        str2double(get(hObject,'String')) returns contents of rv as a double
r = str2double(get(hObject, 'String'));
handles.metricdata.r = r;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function rv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pev_Callback(hObject, eventdata, handles)
% hObject    handle to Pev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pev as text
%        str2double(get(hObject,'String')) returns contents of Pev as a double
Pe = str2double(get(hObject, 'String'));
handles.metricdata.Pe = Pe;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Pev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kccv_Callback(hObject, eventdata, handles)
% hObject    handle to kccv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kccv as text
%        str2double(get(hObject,'String')) returns contents of kccv as a double
kcc = str2double(get(hObject, 'String'));
handles.metricdata.kcc = kcc;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function kccv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kccv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function T3v_Callback(hObject, eventdata, handles)
% hObject    handle to T3v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T3v as text
%        str2double(get(hObject,'String')) returns contents of T3v as a double
T3 = str2double(get(hObject, 'String'));
handles.metricdata.T3 = T3;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function T3v_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T3v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xv_Callback(hObject, eventdata, handles)
% hObject    handle to xv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xv as text
%        str2double(get(hObject,'String')) returns contents of xv as a double
x = str2double(get(hObject, 'String'));
handles.metricdata.x = x;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function xv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yv_Callback(hObject, eventdata, handles)
% hObject    handle to yv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yv as text
%        str2double(get(hObject,'String')) returns contents of yv as a double
y = str2double(get(hObject, 'String'));
handles.metricdata.y = y;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function yv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zv_Callback(hObject, eventdata, handles)
% hObject    handle to zv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zv as text
%        str2double(get(hObject,'String')) returns contents of zv as a double
z = str2double(get(hObject, 'String'));
handles.metricdata.z = z;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function zv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in simulate.
function simulate_Callback(hObject, eventdata, handles)
% hObject    handle to simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% with or without recuperator
if get(handles.recuperator,'Value')==0
   
    %%%% Parameters %%%%
    T1 = handles.metricdata.T1; % [°C]
    r = handles.metricdata.r;
    Pe = handles.metricdata.Pe; % [MW]
    kcc = handles.metricdata.kcc;
    T3 = handles.metricdata.T3; % [°C]
    x = handles.metricdata.x;
    y = handles.metricdata.y;
    z = handles.metricdata.z;
    
    %%%% Fixed data %%%%
    etaPiC = 0.9;
    etaPiT = 0.9;
    
    %%%% Simulation %%%%
    [state,Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten,Exergy_losses,labels_Ex,etaRotex,etaCyclex,etaCombex,etaTotex] = mainTurbineGaz(T1,r,etaPiC,kcc,T3,etaPiT,Pe,x,y,z);
    
else
    
end

function initialize_gui(fig_handle, handles)

% Cycle figure
axes(handles.cycle);
img = imread('TurbineGas.PNG');
image(img);
set(gca,'Visible','off');

% default parameters
handles.metricdata.T1 = 15;
handles.metricdata.r = 18;
handles.metricdata.Pe = 230;
handles.metricdata.kcc = 0.95;
handles.metricdata.T3 = 1400;
handles.metricdata.x = 1;
handles.metricdata.y = 4;
handles.metricdata.z = 0;

set(handles.T1v,'String',handles.metricdata.T1);
set(handles.rv,'String',18);
set(handles.Pev,'String',230);
set(handles.kccv,'String',0.95);
set(handles.T3v,'String',1400);
set(handles.xv,'String',1);
set(handles.yv,'String',4);
set(handles.zv,'String',0);
guidata(handles.figure1, handles);