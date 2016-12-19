function varargout = SCGUIV3(varargin)
% SCGUIV3 MATLAB code for SCGUIV3.fig
%      SCGUIV3, by itself, creates a new SCGUIV3 or raises the existing
%      singleton*.
%
%      H = SCGUIV3 returns the handle to a new SCGUIV3 or the handle to
%      the existing singleton*.
%
%      SCGUIV3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCGUIV3.M with the given input arguments.
%
%      SCGUIV3('Property','Value',...) creates a new SCGUIV3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SCGUIV3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SCGUIV3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SCGUIV3

% Last Modified by GUIDE v2.5 19-Dec-2016 00:58:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SCGUIV3_OpeningFcn, ...
                   'gui_OutputFcn',  @SCGUIV3_OutputFcn, ...
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


% --- Executes just before SCGUIV3 is made visible.
function SCGUIV3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SCGUIV3 (see VARARGIN)

handles.metricdata.pel  = 35;
handles.metricdata.ta   = 15;
handles.metricdata.tmax = 520;
handles.metricdata.texh = 120;
handles.metricdata.pmax = 40;
handles.metricdata.fh   = 0;
handles.metricdata.lambda = 1.05;

axes(handles.axes6);
cla reset;
matlabImage = imread('rankine.PNG');
image(matlabImage)
axis off
axis image

% Choose default command line output for SCGUIV3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SCGUIV3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SCGUIV3_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in reheating.
function reheating_Callback(hObject, eventdata, handles)
% hObject    handle to reheating (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reheating


% --- Executes on button press in pushbuttonsimulate.
function pushbuttonsimulate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonsimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


P_el = 1000*handles.metricdata.pel;
Ta   = handles.metricdata.ta;
Tmax = handles.metricdata.tmax;
Texh = handles.metricdata.texh;
Pmax = handles.metricdata.pmax;
FH   = handles.metricdata.fh;
lambda = handles.metricdata.lambda;

% Fixed Params
Triv = 15;
Tmin = Triv + 18;

if get(handles.reheating, 'value') == 0
    RH = 'off';
    if FH == 0
        axes(handles.axes6);
        cla reset;
        matlabImage = imread('rankine.PNG');
        image(matlabImage)
        axis off
        axis image
    elseif FH > 0 
        axes(handles.axes6);
        cla reset;
        matlabImage = imread('feedheater.PNG');
        image(matlabImage)
        axis off
        axis image 
    end
elseif get(handles.reheating, 'value') == 1
    RH = 'on';
    if FH == 0
        axes(handles.axes6);
        cla reset;
        matlabImage = imread('reheating.PNG');
        image(matlabImage)
        axis off
        axis image
    elseif FH > 0 
        axes(handles.axes6);
        cla reset;
        matlabImage = imread('fhandrh.PNG');
        image(matlabImage)
        axis off
        axis image 
    end
end

%% States of the cycle
Steam_Cycle = SCStates (Tmax, Tmin, Pmax, FH, RH);
handles.steamcycle = Steam_Cycle;

if FH == 0
    A = length(Steam_Cycle);
else 
    A = length(Steam_Cycle) - 1;
end

skipped = 0;
States = cell(A,7);
for i = 1:A(1)
    if Steam_Cycle{i}.p == 0
        skipped = 1;
    end
    if skipped == 0 
        States{i,1} = Steam_Cycle{i}.States;
        States{i,2} = Steam_Cycle{i}.t;
        States{i,3} = 100*Steam_Cycle{i}.p;
        States{i,4} = Steam_Cycle{i}.x;
        States{i,5} = Steam_Cycle{i}.h;
        States{i,6} = Steam_Cycle{i}.s;
        States{i,7} = Steam_Cycle{i}.e;  
    else
        States{i,1} = Steam_Cycle{i+1}.States;
        States{i,2} = Steam_Cycle{i+1}.t;
        States{i,3} = 100*Steam_Cycle{i+1}.p;
        States{i,4} = Steam_Cycle{i+1}.x;
        States{i,5} = Steam_Cycle{i+1}.h;
        States{i,6} = Steam_Cycle{i+1}.s;
        States{i,7} = Steam_Cycle{i+1}.e;      
    end
end

set(handles.states,'Data',States);

%% Analysis
fuel = get(handles.fueltype,'SelectedObject');
if fuel == handles.ch4
    x = 1;
    y = 4;
    z = 0;
    comb = 'CH4';
elseif fuel == handles.c12h23
    x = 12;
    y = 23;
    z = 0;
    comb = 'C12H23';
end

[ FlowRates, eta_en, EnergyDistribution, eta_ex, ExergyDistribution ] = ...
    SCAnalysis (Steam_Cycle, P_el, FH, RH, x, y ,z, Texh, Ta, lambda);

set(handles.flowrates,'Data',FlowRates');

for i = 1:length(EnergyDistribution)
    EnergyDistribution{2,i} = EnergyDistribution{2,i}/1000;
end
for i = 1:length(ExergyDistribution)
    ExergyDistribution{2,i} = ExergyDistribution{2,i}/1000;
end

handles.energydistribution = EnergyDistribution;
handles.exergydistribution = ExergyDistribution;

mecen  = num2str(eta_en{2,1});
gen    = num2str(eta_en{2,2});
cyclen = num2str(eta_en{2,3});
toten  = num2str(eta_en{2,4});
set(handles.texttoten, 'String', toten);
set(handles.textmecen, 'String', mecen);
set(handles.textgen,   'String', gen);
set(handles.textcyclen,'String', cyclen);
%mecex   = num2str(eta_ex{2,1});
totex   = num2str(eta_ex{2,2});
gex     = num2str(eta_ex{2,3});
combex  = num2str(eta_ex{2,4});
chimnex = num2str(eta_ex{2,5});
transex = num2str(eta_ex{2,6});
rotex   = num2str(eta_ex{2,7});
cyclex  = num2str(eta_ex{2,8});
set(handles.texttotex,  'String', totex);
%set(handles.textmecex,  'String', mecex);
set(handles.textgex,    'String', gex);
set(handles.textcyclex, 'String', cyclex);
set(handles.textchimnex,'String', chimnex);
set(handles.textcombex, 'String', combex);
set(handles.texttransex,'String', transex);
set(handles.textrotex,  'String', rotex);

tsdiagramplot(handles);
pieenplot(handles);

guidata(hObject, handles);



function editpower_Callback(hObject, eventdata, handles)
% hObject    handle to editpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editpower as text
%        str2double(get(hObject,'String')) returns contents of editpower as a double
pel = str2double(get(hObject, 'String'));
handles.metricdata.pel = pel;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function editpower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editpmax_Callback(hObject, eventdata, handles)
% hObject    handle to editpmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editpmax as text
%        str2double(get(hObject,'String')) returns contents of editpmax as a double
pmax = str2double(get(hObject, 'String'));
handles.metricdata.pmax= pmax;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function editpmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editpmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editfh_Callback(hObject, eventdata, handles)
% hObject    handle to editfh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editfh as text
%        str2double(get(hObject,'String')) returns contents of editfh as a double
fh = str2double(get(hObject, 'String'));
handles.metricdata.fh = fh;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function editfh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editfh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editta_Callback(hObject, eventdata, handles)
% hObject    handle to editta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editta as text
%        str2double(get(hObject,'String')) returns contents of editta as a double
ta = str2double(get(hObject, 'String'));
handles.metricdata.ta = ta;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function editta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edittexh_Callback(hObject, eventdata, handles)
% hObject    handle to edittexh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittexh as text
%        str2double(get(hObject,'String')) returns contents of edittexh as a double
texh = str2double(get(hObject, 'String'));
handles.metricdata.texh = texh;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edittexh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittexh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edittmax_Callback(hObject, eventdata, handles)
% hObject    handle to edittmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittmax as text
%        str2double(get(hObject,'String')) returns contents of edittmax as a double
tmax = str2double(get(hObject, 'String'));
handles.metricdata.tmax = tmax;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edittmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editlambda_Callback(hObject, eventdata, handles)
% hObject    handle to editlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlambda as text
%        str2double(get(hObject,'String')) returns contents of editlambda as a double
lambda = str2double(get(hObject, 'String'));
handles.metricdata.lambda = lambda;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function editlambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editlambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textcyclex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textcyclex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textgex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textgex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textmecex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textmecex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function texttotex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to texttotex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function textrotex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textrotex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textchimnex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textchimnex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function textcombex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textcombex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function texttransex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to texttransex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function texttoten_CreateFcn(hObject, eventdata, handles)
% hObject    handle to texttoten (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textmecen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textmecen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textgen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function textcyclen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textcyclen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonts.
function pushbuttonts_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tsdiagramplot(handles)

function tsdiagramplot(handles)
axes(handles.axes4)
cla reset 

if get(handles.reheating, 'value') == 0
    RH = 'off';
elseif get(handles.reheating, 'value') == 1
    RH = 'on';
end

%% States of the cycle
FH   = handles.metricdata.fh;
state = handles.steamcycle;

ind1    = 1;
ind2    = 2;
step2 = 500;
eta_SiTHP = 0.90;
eta_SiTLP = 0.88;

switch RH
    case 'off' 
        if FH == 0
           alpha = 6;
           ind3  = 5;
           ind4  = 6;
        elseif FH > 0
           alpha = 10;
           beta  = 4;
           ind3   = 5;
           ind4   = 6;
           ind5   = 7;
           ind6   = 8;
           n = alpha + beta * FH;
        else
            warning('Negative number of feedheaters not allowed')
        end
    case 'on'
        if FH == 0
            alpha = 8;
            ind3HP  = 5;
            ind4HP  = 6;
            ind3LP  = 7;
            ind4LP  = 8;
        elseif FH > 0
            alpha = 12;
            beta  = 4;
            ind3HP  = 5;
            ind4HP  = 6;
            ind3LP  = 7;
            ind4LP  = 8;
            ind5    = 9;
            ind6    = 10;
            n = alpha + beta * FH + 4; 
        else
            warning('Negative number of feedheaters not allowed')
        end
    otherwise
        warning('Unexpected reheating entry.')
end


% Determining the Bache Position
if FH > 4
    index_bleed_bache = alpha + 1;
    position_bache = 1;
    while state{index_bleed_bache}.p < 4.6
        position_bache = position_bache + 1;
        index_bleed_bache = index_bleed_bache + beta;
    end
end

%% Saturated States Curve
S = linspace(0.00139,9.1531,2*step2);
Tsat = arrayfun(@(s) XSteam('Tsat_s',s), S);

%% Main Commom Points
% Curve 1 to 2
s12 = linspace(state{ind1}.s,state{ind2}.s,step2);
p12 = linspace(state{ind1}.p,state{ind2}.p,step2);
t12 = arrayfun(@(p,s) XSteam('T_ps',p,s), p12, s12);

%% RH OFF, FH = 0
if alpha == 6
    % Curve 2 to 3
    s23 = linspace(state{ind2}.s,state{ind3}.s,step2);
    p23 = linspace(state{ind2}.p,state{ind3}.p,step2);
    t23 = arrayfun(@(p,s) XSteam('T_ps',p,s), p23, s23);
    % Curve 3 to 4
    p34   = linspace(state{ind3}.p,state{ind4}.p,step2);
    s34_s = linspace(state{ind3}.s,state{ind3}.s,step2);
    h34_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34,s34_s);
    h34   = zeros(1,step2);
    for i = 1:step2
        h34(i) = state{ind3}.h + eta_SiTLP*(h34_s(i) - state{ind3}.h);
    end
    t34   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34,h34);
    s34   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34,h34);

    % Curve Inside Condensor 
    % 4 to 4L (vapor to sat liquid)
    s4L  = XSteam('sL_p',state{ind4}.p);
    s44L = linspace(state{ind4}.s,s4L,step2);
    p44L = linspace(state{ind4}.p,state{ind4}.p,step2);
    t44L = arrayfun(@(p,s) XSteam('T_ps',p,s), p44L, s44L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L1 = linspace(s4L,state{ind1}.s,step2);
    p4L1 = linspace(state{ind4}.p,state{ind1}.p,step2);
    t4L1 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4L1, s4L1);
    
    s41 = [s44L, s4L1];
    t41 = [t44L, t4L1];
    
Sif = [s23, s34, s41];
Tif = [t23, t34, t41];

%% RH ON, FH = 0
elseif alpha == 8 
    % Main Points
    % Curve 2 to 3HP
    s23 = linspace(state{ind2}.s,state{ind3HP}.s,step2);
    p23 = linspace(state{ind2}.p,state{ind3HP}.p,step2);
    t23 = arrayfun(@(p,s) XSteam('T_ps',p,s), p23, s23);
    
    % Curve 3HP to 4HP
    p34HP   = linspace(state{ind3HP}.p,state{ind4HP}.p,step2);
    s34HP_s = linspace(state{ind3HP}.s,state{ind3HP}.s,step2);
    h34HP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34HP,s34HP_s);
    h34HP   = zeros(1,step2);
    for i = 1:step2
        h34HP(i) = state{ind3HP}.h + eta_SiTHP*(h34HP_s(i) - state{ind3HP}.h);
    end
    t34HP   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34HP,h34HP);
    s34HP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34HP,h34HP);
    % Curve 4HP to 3LP
    s4HP3LP = linspace(state{ind4HP}.s,state{ind3LP}.s,step2);
    p4HP3LP = linspace(state{ind4HP}.p,state{ind3LP}.p,step2);
    t4HP3LP = arrayfun(@(p,s) XSteam('T_ps',p,s), p4HP3LP, s4HP3LP);  
    % Curve 3BP to 4BP
    p34LP   = linspace(state{ind3LP}.p,state{ind4LP}.p,step2);
    s34LP_s = linspace(state{ind3LP}.s,state{ind3LP}.s,step2);
    h34LP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34LP,s34LP_s);
    h34LP   = zeros(1,step2);
    for i = 1:step2
        h34LP(i) = state{ind3LP}.h + eta_SiTLP*(h34LP_s(i) - state{ind3LP}.h);
    end
    t34LP   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34LP,h34LP);
    s34LP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34LP,h34LP);
    
    s34 = [s34HP, s4HP3LP, s34LP];
    t34 = [t34HP, t4HP3LP, t34LP];

    % Curve Inside Condensor 
    % 4LP to 4L (vapor to sat liquid)
    s4LPL  = XSteam('sL_p',state{ind4LP}.p);
    s4LP4L = linspace(state{ind4LP}.s,s4LPL,step2);
    p4LP4L = linspace(state{ind4LP}.p,state{ind4LP}.p,step2);
    t4LP4L = arrayfun(@(p,s) XSteam('T_ps',p,s), p4LP4L, s4LP4L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L1 = linspace(s4LPL,state{ind1}.s,step2);
    p4L1 = linspace(state{ind4LP}.p,state{ind1}.p,step2);
    t4L1 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4L1, s4L1);
    
    s4LP1 = [s4LP4L, s4L1];
    t4LP1 = [t4LP4L, t4L1];
    
    
    Sif = [s23, s34, s4LP1];
    Tif = [t23, t34, t4LP1];
    
%% RH OFF, FH # 0
elseif alpha == 10
    % Main Points
    % Curve 2 to 3
    s23 = linspace(state{ind2}.s,state{ind3}.s,step2);
    p23 = linspace(state{ind2}.p,state{ind3}.p,step2);
    t23 = arrayfun(@(p,s) XSteam('T_ps',p,s), p23, s23);
    % Curve 3 to 4
    p34   = linspace(state{ind3}.p,state{ind4}.p,step2);
    s34_s = linspace(state{ind3}.s,state{ind3}.s,step2);
    h34_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34,s34_s);
    h34   = zeros(1,step2);
    for i = 1:step2
        h34(i) = state{ind3}.h + eta_SiTLP*(h34_s(i) - state{ind3}.h);
    end
    t34   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34,h34);
    s34   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34,h34);
    
    % Curve Inside Condensor 
    % 4 to 4L (vapor to sat liquid)
    s4L = XSteam('sL_p',state{ind4}.p);
    s44L = linspace(state{ind4}.s,s4L,step2);
    p44L = linspace(state{ind4}.p,state{ind4}.p,step2);
    t44L = arrayfun(@(p,s) XSteam('T_ps',p,s), p44L, s44L);
    % 4L to 5 sat liq to undercooled liquid at exit of condensor
    s4L5 = linspace(s4L,state{ind5}.s,step2);
    p4L5 = linspace(state{ind4}.p,state{ind5}.p,step2);
    t4L5 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4L5, s4L5);
    
    s45 = [s44L, s4L5];
    t45 = [t44L, t4L5];
    
    % Condensor Pump 5 to 6
    s56 = linspace(state{ind5}.s,state{ind6}.s,step2);
    p56 = linspace(state{ind5}.p,state{ind6}.p,step2);
    t56 = arrayfun(@(p,s) XSteam('T_ps',p,s), p56, s56);
    % Curves from 6 to 1
    if FH > 4
        index_bleed_bache = alpha + 1;
        position_bache = 1;
        while state{index_bleed_bache}.p < 4.6
            position_bache = position_bache + 1;
            index_bleed_bache = index_bleed_bache + beta;
        end
        % from 6 to bache
        s6bache = linspace(state{ind6}.s,state{index_bleed_bache + 2}.s,step2);
        p6bache = linspace(state{ind6}.p,state{index_bleed_bache + 2}.p,step2);
        t6bache = arrayfun(@(p,s) XSteam('T_ps',p,s), p6bache, s6bache);
        
        % bache pump
        ind6BP = alpha + beta*(position_bache + 1) - 2;
        sBP = linspace(state{index_bleed_bache + 2}.s,state{ind6BP}.s,step2);
        pBP = linspace(state{index_bleed_bache + 2}.p,state{ind6BP}.p,step2);
        tBP = arrayfun(@(p,s) XSteam('T_ps',p,s), pBP, sBP);
        
        % from 6BP to 1
        s6BP1 = linspace(state{ind6BP}.s,state{ind1}.s,step2);
        p6BP1 = linspace(state{ind6BP}.p,state{ind1}.p,step2);
        t6BP1 = arrayfun(@(p,s) XSteam('T_ps',p,s), p6BP1, s6BP1);
        
        s61 = [s6bache, sBP, s6BP1];
        t61 = [t6bache, tBP, t6BP1];
    else
        s61 = linspace(state{ind6}.s,state{ind1}.s,step2);
        p61 = linspace(state{ind6}.p,state{ind1}.p,step2);
        t61 = arrayfun(@(p,s) XSteam('T_ps',p,s), p61, s61); 
    end

    Sif = [s23, s34, s45, s56, s61];
    Tif = [t23, t34, t45, t56, t61];
    
    % Additional Points
    % Curves 4i to 7i 
    s47i = zeros(FH,step2);
    p47i = zeros(FH,step2);
    t47i = zeros(FH,step2);
    for i = 1:FH
        ind4i = alpha + beta*(i-1) + 1;
        ind7i = alpha + beta*i - 1;
        if (FH > 4) && (i == position_bache)
            p_bleedbache = state{ind7i}.p; % Pressure behind the isenthalpic valve 
            h_bleedbache = state{ind4i}.h; 
            s_bleedbache = XSteam('s_ph',p_bleedbache,h_bleedbache);
            s47i(i,:) = linspace(s_bleedbache,state{ind7i}.s,step2);
            p47i(i,:) = linspace(p_bleedbache,state{ind7i}.p,step2);
            t47i(i,:) = arrayfun(@(p,s) XSteam('T_ps',p,s), p47i(i,:), s47i(i,:)); 
            
            s4valve = linspace(state{ind4i}.s,s_bleedbache,step2);
            p4valve = linspace(state{ind4i}.p,p_bleedbache,step2);
            t4valve = arrayfun(@(p,s) XSteam('T_ps',p,s), p4valve, s4valve);
        else
        s47i(i,:) = linspace(state{ind4i}.s,state{ind7i}.s,step2);
        p47i(i,:) = linspace(state{ind4i}.p,state{ind7i}.p,step2);
        t47i(i,:) = arrayfun(@(p,s) XSteam('T_ps',p,s), p47i(i,:), s47i(i,:)); 
        end
    end
    % Additional curves were intentionally left out because they would not
    % bring much interesting information and would render the calculations
    % heavier.
    
    
%% RH ON, FH # 0
elseif alpha == 12
    % Main Points
    % Curve 2 to 3HP
    s23 = linspace(state{ind2}.s,state{ind3HP}.s,step2);
    p23 = linspace(state{ind2}.p,state{ind3HP}.p,step2);
    t23 = arrayfun(@(p,s) XSteam('T_ps',p,s), p23, s23);
    
    % Curve 3HP to 4HP
    p34HP   = linspace(state{ind3HP}.p,state{ind4HP}.p,step2);
    s34HP_s = linspace(state{ind3HP}.s,state{ind3HP}.s,step2);
    h34HP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34HP,s34HP_s);
    h34HP   = zeros(1,step2);
    for i = 1:step2
        h34HP(i) = state{ind3HP}.h + eta_SiTHP*(h34HP_s(i) - state{ind3HP}.h);
    end
    t34HP   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34HP,h34HP);
    s34HP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34HP,h34HP);
    % Curve 4HP to 3LP
    s4HP3LP = linspace(state{ind4HP}.s,state{ind3LP}.s,step2);
    p4HP3LP = linspace(state{ind4HP}.p,state{ind3LP}.p,step2);
    t4HP3LP = arrayfun(@(p,s) XSteam('T_ps',p,s), p4HP3LP, s4HP3LP);  
    % Curve 3BP to 4BP
    p34LP   = linspace(state{ind3LP}.p,state{ind4LP}.p,step2);
    s34LP_s = linspace(state{ind3LP}.s,state{ind3LP}.s,step2);
    h34LP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34LP,s34LP_s);
    h34LP   = zeros(1,step2);
    for i = 1:step2
        h34LP(i) = state{ind3LP}.h + eta_SiTLP*(h34LP_s(i) - state{ind3LP}.h);
    end
    t34LP   = arrayfun(@(p,h) XSteam('T_ph',p,h),p34LP,h34LP);
    s34LP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34LP,h34LP);
        
    s34 = [s34HP, s4HP3LP, s34LP];
    t34 = [t34HP, t4HP3LP, t34LP];
    
    % Curve Inside Condensor 
    % 4LP to 4L (vapor to sat liquid)
    s4LPL  = XSteam('sL_p',state{ind4LP}.p);
    s4LP4L = linspace(state{ind4LP}.s,s4LPL,step2);
    p4LP4L = linspace(state{ind4LP}.p,state{ind4LP}.p,step2);
    t4LP4L = arrayfun(@(p,s) XSteam('T_ps',p,s), p4LP4L, s4LP4L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L5 = linspace(s4LPL,state{ind5}.s,step2);
    p4L5 = linspace(state{ind4LP}.p,state{ind5}.p,step2);
    t4L5 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4L5, s4L5);
    
    s45 = [s4LP4L, s4L5];
    t45 = [t4LP4L, t4L5];
    
    % Condensor Pump 1
    s56 = linspace(state{ind5}.s,state{ind6}.s,step2);
    p56 = linspace(state{ind5}.p,state{ind6}.p,step2);
    t56 = arrayfun(@(p,s) XSteam('T_ps',p,s), p56, s56);
    % Curves from 6 to 1
    if FH > 4
        % from 6 to bache
        s6bache = linspace(state{ind6}.s,state{index_bleed_bache + 2}.s,step2);
        p6bache = linspace(state{ind6}.p,state{index_bleed_bache + 2}.p,step2);
        t6bache = arrayfun(@(p,s) XSteam('T_ps',p,s), p6bache, s6bache);
        
        % bache pump
        ind6BP = alpha + beta*(position_bache + 1) - 2;
        sBP = linspace(state{index_bleed_bache + 2}.s,state{ind6BP}.s,step2);
        pBP = linspace(state{index_bleed_bache + 2}.p,state{ind6BP}.p,step2);
        tBP = arrayfun(@(p,s) XSteam('T_ps',p,s), pBP, sBP);
        
        % from 6BP to 1
        s6BP1 = linspace(state{ind6BP}.s,state{ind1}.s,step2);
        p6BP1 = linspace(state{ind6BP}.p,state{ind1}.p,step2);
        t6BP1 = arrayfun(@(p,s) XSteam('T_ps',p,s), p6BP1, s6BP1);
        
        s61 = [s6bache, sBP, s6BP1];
        t61 = [t6bache, tBP, t6BP1];
    else
        s61 = linspace(state{ind6}.s,state{ind1}.s,step2);
        p61 = linspace(state{ind6}.p,state{ind1}.p,step2);
        t61 = arrayfun(@(p,s) XSteam('T_ps',p,s), p61, s61); 
    end
    
    Sif = [s23, s34, s45, s56, s61];
    Tif = [t23, t34, t45, t56, t61];
    
    % Additional Points
    % Curve 4HP to 7HP
    s4HP7 = linspace(state{ind4HP}.s,state{n - beta + 3}.s,step2);
    p4HP7 = linspace(state{ind4HP}.p,state{n - beta + 3}.p,step2);
    t4HP7 = arrayfun(@(p,s) XSteam('T_ps',p,s), p4HP7, s4HP7);

    % Curves 4i to 7i 
    s47i = zeros(FH + 1,step2);
    p47i = zeros(FH + 1,step2);
    t47i = zeros(FH + 1,step2);
    for i = 1:FH
        ind4i = alpha + beta*(i-1) + 1;
        ind7i = alpha + beta*i - 1;
        if (FH > 4) && (i == position_bache)
            p_bleedbache = state{ind7i}.p; % Pressure behind the isenthalpic valve 
            h_bleedbache = state{ind4i}.h; 
            s_bleedbache = XSteam('s_ph',p_bleedbache,h_bleedbache);
            s47i(i,:) = linspace(s_bleedbache,state{ind7i}.s,step2);
            p47i(i,:) = linspace(p_bleedbache,state{ind7i}.p,step2);
            t47i(i,:) = arrayfun(@(p,s) XSteam('T_ps',p,s), p47i(i,:), s47i(i,:)); 
            
            s4valve = linspace(state{ind4i}.s,s_bleedbache,step2);
            p4valve = linspace(state{ind4i}.p,p_bleedbache,step2);
            t4valve = arrayfun(@(p,s) XSteam('T_ps',p,s), p4valve, s4valve);
        else
        s47i(i,:) = linspace(state{ind4i}.s,state{ind7i}.s,step2);
        p47i(i,:) = linspace(state{ind4i}.p,state{ind7i}.p,step2);
        t47i(i,:) = arrayfun(@(p,s) XSteam('T_ps',p,s), p47i(i,:), s47i(i,:)); 
        end
    end
    s47i(FH+1,:) = s4HP7;
    t47i(FH+1,:) = t4HP7;
    % Additional curves were intentionally left out because they would not
    % bring much interesting information and would render the calculations
    % heavier.
    
end


%% Individual Cycle Points
A = size(state);
if FH > 0
    T_Cycle = ones(1,A(1) - 1);
    S_Cycle = ones(1,A(1) - 1);
    for i=1:A(1) - 1
        if i < alpha + 2
            T_Cycle(i) = state{i}.t;
            S_Cycle(i) = state{i}.s;
        elseif i >= alpha + 2
            T_Cycle(i) = state{i + 1}.t;
            S_Cycle(i) = state{i + 1}.s;
        end
    end
else
    T_Cycle = ones(1,A(1));
    S_Cycle = ones(1,A(1));
    for i = 1:A(1)
        T_Cycle(i) = state{i}.t;
        S_Cycle(i) = state{i}.s;
    end
end



Smain = [s12, Sif];
Tmain = [t12, Tif];


axes(handles.axes4)
plot(S_Cycle,T_Cycle,'*r',S,Tsat,'-b',Smain,Tmain,'-r'); hold on
if FH > 0
    size_s47i = size(s47i);
    for i = 1:size_s47i(1)
        plot(s47i(i,:),t47i(i,:),'-.g'); hold on
    end
    if FH > 4
        plot(s4valve,t4valve,'-.r'); hold on
    end
end
xlabel('s [kJ/kgK]'); hold on;
ylabel('t [°C]'); hold on;
grid on; 

% --- Executes on button press in pushbuttonhs.
function pushbuttonhs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonhs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hsdiagramplot(handles)

function hsdiagramplot(handles)
axes(handles.axes4);
cla reset;

if get(handles.reheating, 'value') == 0
    RH = 'off';
elseif get(handles.reheating, 'value') == 1
    RH = 'on';
end

%% States of the cycle
FH   = handles.metricdata.fh;
state = handles.steamcycle;

ind1    = 1;
ind2    = 2;
step2 = 500;
eta_SiTHP = 0.90;
eta_SiTLP = 0.88;

switch RH
    case 'off' 
        if FH == 0
           alpha = 6;
           ind3  = 5;
           ind4  = 6;
        elseif FH > 0
           alpha = 10;
           beta  = 4;
           ind3   = 5;
           ind4   = 6;
           ind5   = 7;
           ind6   = 8;
           n = alpha + beta * FH;
        else
            warning('Negative number of feedheaters not allowed')
        end
    case 'on'
        if FH == 0
            alpha = 8;
            ind3HP  = 5;
            ind4HP  = 6;
            ind3LP  = 7;
            ind4LP  = 8;
        elseif FH > 0
            alpha = 12;
            beta  = 4;
            ind3HP  = 5;
            ind4HP  = 6;
            ind3LP  = 7;
            ind4LP  = 8;
            ind5    = 9;
            ind6    = 10;
            n = alpha + beta * FH + 4; 
        else
            warning('Negative number of feedheaters not allowed')
        end
    otherwise
        warning('Unexpected reheating entry.')
end


% Determining the Bache Position
if FH > 4
    index_bleed_bache = alpha + 1;
    position_bache = 1;
    while state{index_bleed_bache}.p < 4.6
        position_bache = position_bache + 1;
        index_bleed_bache = index_bleed_bache + beta;
    end
end

%% Saturated States Curve
S = linspace(0.00139,9.1531,2*step2);
Psat = arrayfun(@(s) XSteam('psat_s',s), S); 
Hsat = arrayfun(@(p,s) XSteam('h_ps',p,s), Psat, S); % entalpie saturee pour chaque entropie

%% Main Commom Points
% Curve 1 to 2
s12 = linspace(state{ind1}.s,state{ind2}.s,step2);
p12 = linspace(state{ind1}.p,state{ind2}.p,step2);
h12 = arrayfun(@(p,s) XSteam('h_ps',p,s), p12, s12);

%% RH OFF, FH = 0
if alpha == 6
    % Curve 2 to 3
    s23 = linspace(state{ind2}.s,state{ind3}.s,step2);
    p23 = linspace(state{ind2}.p,state{ind3}.p,step2);
    h23 = arrayfun(@(p,s) XSteam('h_ps',p,s), p23, s23);
    % Curve 3 to 4
    p34   = linspace(state{ind3}.p,state{ind4}.p,step2);
    s34_s = linspace(state{ind3}.s,state{ind3}.s,step2);
    h34_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34,s34_s);
    h34   = zeros(1,step2);
    for i = 1:step2
        h34(i) = state{ind3}.h + eta_SiTLP*(h34_s(i) - state{ind3}.h);
    end
    s34   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34,h34);

    % Curve Inside Condensor 
    % 4 to 4L (vapor to sat liquid)
    s4L  = XSteam('sL_p',state{ind4}.p);
    s44L = linspace(state{ind4}.s,s4L,step2);
    p44L = linspace(state{ind4}.p,state{ind4}.p,step2);
    h44L = arrayfun(@(p,s) XSteam('h_ps',p,s), p44L, s44L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L1 = linspace(s4L,state{ind1}.s,step2);
    p4L1 = linspace(state{ind4}.p,state{ind1}.p,step2);
    h4L1 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4L1, s4L1);
    
    s41 = [s44L, s4L1];
    h41 = [h44L, h4L1];
    
Sif = [s23, s34, s41];
Hif = [h23, h34, h41];

%% RH ON, FH = 0
elseif alpha == 8 
    % Main Points
    % Curve 2 to 3HP
    s23 = linspace(state{ind2}.s,state{ind3HP}.s,step2);
    p23 = linspace(state{ind2}.p,state{ind3HP}.p,step2);
    h23 = arrayfun(@(p,s) XSteam('h_ps',p,s), p23, s23);
    
    % Curve 3HP to 4HP
    p34HP   = linspace(state{ind3HP}.p,state{ind4HP}.p,step2);
    s34HP_s = linspace(state{ind3HP}.s,state{ind3HP}.s,step2);
    h34HP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34HP,s34HP_s);
    h34HP   = zeros(1,step2);
    for i = 1:step2
        h34HP(i) = state{ind3HP}.h + eta_SiTHP*(h34HP_s(i) - state{ind3HP}.h);
    end
    s34HP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34HP,h34HP);
    % Curve 4HP to 3LP
    s4HP3LP = linspace(state{ind4HP}.s,state{ind3LP}.s,step2);
    p4HP3LP = linspace(state{ind4HP}.p,state{ind3LP}.p,step2);
    h4HP3LP = arrayfun(@(p,s) XSteam('h_ps',p,s), p4HP3LP, s4HP3LP); 
    % Curve 3BP to 4BP
    p34LP   = linspace(state{ind3LP}.p,state{ind4LP}.p,step2);
    s34LP_s = linspace(state{ind3LP}.s,state{ind3LP}.s,step2);
    h34LP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34LP,s34LP_s);
    h34LP   = zeros(1,step2);
    for i = 1:step2
        h34LP(i) = state{ind3LP}.h + eta_SiTLP*(h34LP_s(i) - state{ind3LP}.h);
    end
    s34LP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34LP,h34LP);
    
    s34 = [s34HP, s4HP3LP, s34LP];
    h34 = [h34HP, h4HP3LP, h34LP];

    % Curve Inside Condensor 
    % 4LP to 4L (vapor to sat liquid)
    s4LPL  = XSteam('sL_p',state{ind4LP}.p);
    s4LP4L = linspace(state{ind4LP}.s,s4LPL,step2);
    p4LP4L = linspace(state{ind4LP}.p,state{ind4LP}.p,step2);
    h4LP4L = arrayfun(@(p,s) XSteam('h_ps',p,s), p4LP4L, s4LP4L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L1 = linspace(s4LPL,state{ind1}.s,step2);
    p4L1 = linspace(state{ind4LP}.p,state{ind1}.p,step2);
    h4L1 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4L1, s4L1);
    
    s4LP1 = [s4LP4L, s4L1];
    h4LP1 = [h4LP4L, h4L1];
    
    
    Sif = [s23, s34, s4LP1];
    Hif = [h23, h34, h4LP1];
    
%% RH OFF, FH # 0
elseif alpha == 10
    % Main Points
    % Curve 2 to 3
    s23 = linspace(state{ind2}.s,state{ind3}.s,step2);
    p23 = linspace(state{ind2}.p,state{ind3}.p,step2);
    h23 = arrayfun(@(p,s) XSteam('h_ps',p,s), p23, s23);
    % Curve 3 to 4
    p34   = linspace(state{ind3}.p,state{ind4}.p,step2);
    s34_s = linspace(state{ind3}.s,state{ind3}.s,step2);
    h34_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34,s34_s);
    h34   = zeros(1,step2);
    for i = 1:step2
        h34(i) = state{ind3}.h + eta_SiTLP*(h34_s(i) - state{ind3}.h);
    end
    s34   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34,h34);
    
    % Curve Inside Condensor 
    % 4 to 4L (vapor to sat liquid)
    s4L = XSteam('sL_p',state{ind4}.p);
    s44L = linspace(state{ind4}.s,s4L,step2);
    p44L = linspace(state{ind4}.p,state{ind4}.p,step2);
    h44L = arrayfun(@(p,s) XSteam('h_ps',p,s), p44L, s44L);
    % 4L to 5 sat liq to undercooled liquid at exit of condensor
    s4L5 = linspace(s4L,state{ind5}.s,step2);
    p4L5 = linspace(state{ind4}.p,state{ind5}.p,step2);
    h4L5 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4L5, s4L5);
    
    s45 = [s44L, s4L5];
    h45 = [h44L, h4L5];
    
    % Condensor Pump 5 to 6
    s56 = linspace(state{ind5}.s,state{ind6}.s,step2);
    p56 = linspace(state{ind5}.p,state{ind6}.p,step2);
    h56 = arrayfun(@(p,s) XSteam('h_ps',p,s), p56, s56);
    % Curves from 6 to 1
    if FH > 4
        index_bleed_bache = alpha + 1;
        position_bache = 1;
        while state{index_bleed_bache}.p < 4.6
            position_bache = position_bache + 1;
            index_bleed_bache = index_bleed_bache + beta;
        end
        % from 6 to bache
        s6bache = linspace(state{ind6}.s,state{index_bleed_bache + 2}.s,step2);
        p6bache = linspace(state{ind6}.p,state{index_bleed_bache + 2}.p,step2);
        h6bache = arrayfun(@(p,s) XSteam('h_ps',p,s), p6bache, s6bache);
        
        % bache pump
        ind6BP = alpha + beta*(position_bache + 1) - 2;
        sBP = linspace(state{index_bleed_bache + 2}.s,state{ind6BP}.s,step2);
        pBP = linspace(state{index_bleed_bache + 2}.p,state{ind6BP}.p,step2);
        hBP = arrayfun(@(p,s) XSteam('h_ps',p,s), pBP, sBP);
        
        % from 6BP to 1
        s6BP1 = linspace(state{ind6BP}.s,state{ind1}.s,step2);
        p6BP1 = linspace(state{ind6BP}.p,state{ind1}.p,step2);
        h6BP1 = arrayfun(@(p,s) XSteam('h_ps',p,s), p6BP1, s6BP1);
        
        s61 = [s6bache, sBP, s6BP1];
        h61 = [h6bache, hBP, h6BP1];
    else
        s61 = linspace(state{ind6}.s,state{ind1}.s,step2);
        p61 = linspace(state{ind6}.p,state{ind1}.p,step2);
        h61 = arrayfun(@(p,s) XSteam('h_ps',p,s), p61, s61);
    end

    Sif = [s23, s34, s45, s56, s61];
    Hif = [h23, h34, h45, h56, h61];
    
    % Additional Points
    % Curves 4i to 7i 
    s47i = zeros(FH,step2);
    p47i = zeros(FH,step2);
    h47i = zeros(FH,step2);
    for i = 1:FH
        ind4i = alpha + beta*(i-1) + 1;
        ind7i = alpha + beta*i - 1;
        if (FH > 4) && (i == position_bache)
            p_bleedbache = state{ind7i}.p; % Pressure behind the isenthalpic valve 
            h_bleedbache = state{ind4i}.h; 
            s_bleedbache = XSteam('s_ph',p_bleedbache,h_bleedbache);
            s47i(i,:) = linspace(s_bleedbache,state{ind7i}.s,step2);
            p47i(i,:) = linspace(p_bleedbache,state{ind7i}.p,step2);
            h47i(i,:) = arrayfun(@(p,s) XSteam('h_ps',p,s), p47i(i,:), s47i(i,:)); 
            
            s4valve = linspace(state{ind4i}.s,s_bleedbache,step2);
            p4valve = linspace(state{ind4i}.p,p_bleedbache,step2);
            h4valve = arrayfun(@(p,s) XSteam('h_ps',p,s), p4valve, s4valve);
        else
        s47i(i,:) = linspace(state{ind4i}.s,state{ind7i}.s,step2);
        p47i(i,:) = linspace(state{ind4i}.p,state{ind7i}.p,step2);
        h47i(i,:) = arrayfun(@(p,s) XSteam('h_ps',p,s), p47i(i,:), s47i(i,:));
        end
    end
    % Additional curves were intentionally left out because they would not
    % bring much interesting information and would render the calculations
    % heavier.
    
    
%% RH ON, FH # 0
elseif alpha == 12
    % Main Points
    % Curve 2 to 3HP
    s23 = linspace(state{ind2}.s,state{ind3HP}.s,step2);
    p23 = linspace(state{ind2}.p,state{ind3HP}.p,step2);
    h23 = arrayfun(@(p,s) XSteam('h_ps',p,s), p23, s23);
    
    % Curve 3HP to 4HP
    p34HP   = linspace(state{ind3HP}.p,state{ind4HP}.p,step2);
    s34HP_s = linspace(state{ind3HP}.s,state{ind3HP}.s,step2);
    h34HP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34HP,s34HP_s);
    h34HP   = zeros(1,step2);
    for i = 1:step2
        h34HP(i) = state{ind3HP}.h + eta_SiTHP*(h34HP_s(i) - state{ind3HP}.h);
    end
    s34HP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34HP,h34HP);
    % Curve 4HP to 3LP
    s4HP3LP = linspace(state{ind4HP}.s,state{ind3LP}.s,step2);
    p4HP3LP = linspace(state{ind4HP}.p,state{ind3LP}.p,step2);
    h4HP3LP = arrayfun(@(p,s) XSteam('h_ps',p,s), p4HP3LP, s4HP3LP);  
    % Curve 3BP to 4BP
    p34LP   = linspace(state{ind3LP}.p,state{ind4LP}.p,step2);
    s34LP_s = linspace(state{ind3LP}.s,state{ind3LP}.s,step2);
    h34LP_s = arrayfun(@(p,s) XSteam('h_ps',p,s),p34LP,s34LP_s);
    h34LP   = zeros(1,step2);
    for i = 1:step2
        h34LP(i) = state{ind3LP}.h + eta_SiTLP*(h34LP_s(i) - state{ind3LP}.h);
    end
    s34LP   = arrayfun(@(p,h) XSteam('s_ph',p,h),p34LP,h34LP);
        
    s34 = [s34HP, s4HP3LP, s34LP];
    h34 = [h34HP, h4HP3LP, h34LP];
    
    % Curve Inside Condensor 
    % 4LP to 4L (vapor to sat liquid)
    s4LPL  = XSteam('sL_p',state{ind4LP}.p);
    s4LP4L = linspace(state{ind4LP}.s,s4LPL,step2);
    p4LP4L = linspace(state{ind4LP}.p,state{ind4LP}.p,step2);
    h4LP4L = arrayfun(@(p,s) XSteam('h_ps',p,s), p4LP4L, s4LP4L);
    % 4L to 1 sat liq to undercooled liquid at exit of condensor
    s4L5 = linspace(s4LPL,state{ind5}.s,step2);
    p4L5 = linspace(state{ind4LP}.p,state{ind5}.p,step2);
    h4L5 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4L5, s4L5);
    
    s45 = [s4LP4L, s4L5];
    h45 = [h4LP4L, h4L5];
    
    % Condensor Pump 1
    s56 = linspace(state{ind5}.s,state{ind6}.s,step2);
    p56 = linspace(state{ind5}.p,state{ind6}.p,step2);
    h56 = arrayfun(@(p,s) XSteam('h_ps',p,s), p56, s56);
    % Curves from 6 to 1
    if FH > 4
        % from 6 to bache
        s6bache = linspace(state{ind6}.s,state{index_bleed_bache + 2}.s,step2);
        p6bache = linspace(state{ind6}.p,state{index_bleed_bache + 2}.p,step2);
        h6bache = arrayfun(@(p,s) XSteam('h_ps',p,s), p6bache, s6bache);
        
        % bache pump
        ind6BP = alpha + beta*(position_bache + 1) - 2;
        sBP = linspace(state{index_bleed_bache + 2}.s,state{ind6BP}.s,step2);
        pBP = linspace(state{index_bleed_bache + 2}.p,state{ind6BP}.p,step2);
        hBP = arrayfun(@(p,s) XSteam('h_ps',p,s), pBP, sBP);
        
        % from 6BP to 1
        s6BP1 = linspace(state{ind6BP}.s,state{ind1}.s,step2);
        p6BP1 = linspace(state{ind6BP}.p,state{ind1}.p,step2);
        h6BP1 = arrayfun(@(p,s) XSteam('h_ps',p,s), p6BP1, s6BP1);
        
        s61 = [s6bache, sBP, s6BP1];
        h61 = [h6bache, hBP, h6BP1];
    else
        s61 = linspace(state{ind6}.s,state{ind1}.s,step2);
        p61 = linspace(state{ind6}.p,state{ind1}.p,step2);
        h61 = arrayfun(@(p,s) XSteam('h_ps',p,s), p61, s61);
    end
    
    Sif = [s23, s34, s45, s56, s61];
    Hif = [h23, h34, h45, h56, h61];
    
    % Additional Points
    % Curve 4HP to 7HP
    s4HP7 = linspace(state{ind4HP}.s,state{n - beta + 3}.s,step2);
    p4HP7 = linspace(state{ind4HP}.p,state{n - beta + 3}.p,step2);
    h4HP7 = arrayfun(@(p,s) XSteam('h_ps',p,s), p4HP7, s4HP7);

    % Curves 4i to 7i 
    s47i = zeros(FH + 1,step2);
    p47i = zeros(FH + 1,step2);
    h47i = zeros(FH + 1,step2);
    for i = 1:FH
        ind4i = alpha + beta*(i-1) + 1;
        ind7i = alpha + beta*i - 1;
        if (FH > 4) && (i == position_bache)
            p_bleedbache = state{ind7i}.p; % Pressure behind the isenthalpic valve 
            h_bleedbache = state{ind4i}.h; 
            s_bleedbache = XSteam('s_ph',p_bleedbache,h_bleedbache);
            s47i(i,:) = linspace(s_bleedbache,state{ind7i}.s,step2);
            p47i(i,:) = linspace(p_bleedbache,state{ind7i}.p,step2);
            h47i(i,:) = arrayfun(@(p,s) XSteam('h_ps',p,s), p47i(i,:), s47i(i,:));
            
            s4valve = linspace(state{ind4i}.s,s_bleedbache,step2);
            p4valve = linspace(state{ind4i}.p,p_bleedbache,step2);
            h4valve = arrayfun(@(p,s) XSteam('h_ps',p,s), p4valve, s4valve);
        else
        s47i(i,:) = linspace(state{ind4i}.s,state{ind7i}.s,step2);
        p47i(i,:) = linspace(state{ind4i}.p,state{ind7i}.p,step2);
        h47i(i,:) = arrayfun(@(p,s) XSteam('h_ps',p,s), p47i(i,:), s47i(i,:));
        end
    end
    s47i(FH+1,:) = s4HP7;
    h47i(FH+1,:) = h4HP7;
    % Additional curves were intentionally left out because they would not
    % bring much interesting information and would render the calculations
    % heavier.
    
end


%% Individual Cycle Points
A = size(state);
if FH > 0
    H_Cycle = ones(1,A(1) - 1);
    S_Cycle = ones(1,A(1) - 1);
    for i=1:A(1) - 1
        if i < alpha + 2
            H_Cycle(i) = state{i}.h;
            S_Cycle(i) = state{i}.s;
        elseif i >= alpha + 2
            H_Cycle(i) = state{i + 1}.h;
            S_Cycle(i) = state{i + 1}.s;
        end
    end
else
    H_Cycle = ones(1,A(1));
    S_Cycle = ones(1,A(1));
    for i = 1:A(1)
        H_Cycle(i) = state{i}.h;
        S_Cycle(i) = state{i}.s;
    end
end



Smain = [s12, Sif];
Hmain = [h12, Hif];


axes(handles.axes4);
plot(S_Cycle,H_Cycle,'*r',S,Hsat,'-b',Smain,Hmain,'-r'); hold on;
if FH > 0
    size_s47i = size(s47i);
    for i = 1:size_s47i(1)
        plot(s47i(i,:),h47i(i,:),'-.g'); hold on;
    end
    if FH > 4
        plot(s4valve,h4valve,'-.r'); hold on;
    end
end
xlabel('s [kJ/kgK]'); hold on;
ylabel('h [kJ/kg]'); hold on;
grid on;

% --- Executes on button press in pushbuttonpieex.
function pushbuttonpieex_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonpieex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pieexplot(handles)

function pieexplot(handles)
ExergyDistribution = handles.exergydistribution;
axes(handles.axes3);
cla reset;

FH   = handles.metricdata.fh;

formatSpec = '%s %.2f';
str1   = ExergyDistribution{1,1};
value1 = ExergyDistribution{2,1};
label1 = sprintf(formatSpec,str1,value1);
str2   = ExergyDistribution{1,2};
value2 = ExergyDistribution{2,2};
label2 = sprintf(formatSpec,str2,value2);
str3   = ExergyDistribution{1,3};
value3 = ExergyDistribution{2,3};
label3 = sprintf(formatSpec,str3,value3);
str4   = ExergyDistribution{1,4};
value4 = ExergyDistribution{2,4};
label4 = sprintf(formatSpec,str4,value4);  
str5   = ExergyDistribution{1,5};
value5 = ExergyDistribution{2,5};
label5 = sprintf(formatSpec,str5,value5);
str6   = ExergyDistribution{1,6};
value6 = ExergyDistribution{2,6};
label6 = sprintf(formatSpec,str6,value6);
str7   = ExergyDistribution{1,7};
value7 = ExergyDistribution{2,7};
label7 = sprintf(formatSpec,str7,value7);

if FH > 0
    str8   = ExergyDistribution{1,8};
    value8 = ExergyDistribution{2,8};
    label8 = sprintf(formatSpec,str8,value8);  
    str9   = ExergyDistribution{1,9};
    value9 = ExergyDistribution{2,9};
    label9 = sprintf(formatSpec,str9,value9);  
     
    Labels = {label1 label2 label3 label4 label5 label6 label7 label8 label9};
else
    Labels = {label1 label2 label3 label4 label5 label6 label7};
end

PrimaryExSum = 0;
for i = 1:length(ExergyDistribution)
    PrimaryExSum = PrimaryExSum + ExergyDistribution{2,i};
end

h = pie( [ExergyDistribution{2,:}],Labels);

hText = findobj(h,'Type','text'); % text object handles

oldExtents_cell = get(hText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell); % numeric array

signValuesx = sign(oldExtents(:,1));
signValuesy = sign(oldExtents(:,2));

textPositions_cell = get(hText,{'Position'}); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array

% X coordinate offset
if FH > 0
    offsetx = 0.64*signValuesx;
    textPositions(1,1) = textPositions(1,1) + offsetx(1);
    textPositions(2,1) = textPositions(2,1) + 1.2*offsetx(2);
    textPositions(3,1) = textPositions(3,1) + offsetx(3);
    textPositions(4,1) = textPositions(4,1) + offsetx(1);
    textPositions(5,1) = textPositions(5,1);
    textPositions(6,1) = textPositions(6,1) + 0.8*offsetx(6);
    textPositions(7,1) = textPositions(7,1) + 0.9*offsetx(7);
    textPositions(8,1) = textPositions(8,1) + offsetx(8);
    textPositions(9,1) = textPositions(9,1) + offsetx(9);
else
    offsetx = 0.64*signValuesx;    
    textPositions(1,1) = textPositions(1,1) + offsetx(1);
    textPositions(2,1) = textPositions(2,1) + offsetx(2);
    textPositions(3,1) = textPositions(3,1) + offsetx(3);
    textPositions(4,1) = textPositions(4,1) + offsetx(1);
    textPositions(5,1) = textPositions(5,1);
    textPositions(6,1) = textPositions(6,1) + 0.8*offsetx(6);
    textPositions(7,1) = textPositions(7,1) + 0.9*offsetx(7);
end

% Y Coordinate offset
if FH > 0
    offsety = 0.1*signValuesy;
    textPositions(1,2) = textPositions(1,2) + offsety(1);
    textPositions(2,2) = textPositions(2,2) - offsety(2);
    textPositions(3,2) = textPositions(3,2);
    textPositions(4,2) = textPositions(4,2) + offsety(4);
    diff = abs(textPositions(4,2) - textPositions(5,2));
    if  diff < 0.5
        textPositions(5,2) = textPositions(4,2) - 0.2;
    else
        textPositions(5,2) = textPositions(5,2) + 2.5*offsety(5);
    end
    textPositions(6,2) = textPositions(6,2) + 3*offsety(6);
    textPositions(7,2) = textPositions(7,2) + offsety(7);
    textPositions(8,2) = textPositions(8,2) - offsety(8);
    textPositions(9,2) = textPositions(9,2) + offsety(9);
else
    offsety = 0.1*signValuesy;
    textPositions(1,2) = textPositions(1,2) + 2*offsety(1);
    textPositions(2,2) = textPositions(2,2);
    textPositions(3,2) = textPositions(3,2) + offsety(3);
    textPositions(4,2) = textPositions(4,2) + offsety(4);
    diff = abs(textPositions(4,2) - textPositions(5,2));
    if  diff < 0.5
        textPositions(5,2) = textPositions(4,2) - 0.2;
    else
        textPositions(5,2) = textPositions(5,2) + 2*offsety(5);
    end
    textPositions(6,2) = textPositions(6,2);
    textPositions(7,2) = textPositions(7,2) + 2*offsety(7);
end




hText(1).Position = textPositions(1,:);
hText(2).Position = textPositions(2,:);
hText(3).Position = textPositions(3,:);
hText(4).Position = textPositions(4,:);
hText(5).Position = textPositions(5,:);
hText(6).Position = textPositions(6,:);
hText(7).Position = textPositions(7,:);
if FH > 0
    hText(8).Position = textPositions(8,:);
    hText(9).Position = textPositions(9,:);
end
PrimSum = strcat({'Primary flux'},{' '},{num2str(round(PrimaryExSum*10)/10)},...
            {' '},{'MW'});
t = title(PrimSum,'FontSize',12,'FontWeight','bold');
        pos = get(t,'position');
        set(t, 'position', pos+[0 0.15 0]);
        
set(handles.axes3,'CameraViewAngle',10.9074493116677,'DataAspectRatio',[1 1 1],...
    'PlotBoxAspectRatio',[1.3 1.3 1]);

% --- Executes on button press in pushbuttonpiieen.
function pushbuttonpiieen_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonpiieen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

EnergyDistribution = handles.energydistribution;
pieenplot(handles);

function pieenplot(handles)
% Function make the Energy Pie plot

EnergyDistribution = handles.energydistribution;
axes(handles.axes3);
cla reset;

formatSpec = '%s %.2f';
str1   = EnergyDistribution{1,1};
value1 = EnergyDistribution{2,1};
label1 = sprintf(formatSpec,str1,value1);
str2   = EnergyDistribution{1,2};
value2 = EnergyDistribution{2,2};
label2 = sprintf(formatSpec,str2,value2);
str3   = EnergyDistribution{1,3};
value3 = EnergyDistribution{2,3};
label3 = sprintf(formatSpec,str3,value3);
str4   = EnergyDistribution{1,4};
value4 = EnergyDistribution{2,4};
label4 = sprintf(formatSpec,str4,value4);  

PrimaryEnSum = 0;
for i = 1:length(EnergyDistribution) - 1
    PrimaryEnSum = PrimaryEnSum + EnergyDistribution{2,i};
end

        
Labels = {label1 label2 label3 label4};
h = pie( [EnergyDistribution{2,1:4}],Labels);

hText = findobj(h,'Type','text'); % text object handles

oldExtents_cell = get(hText,'Extent'); % cell array
oldExtents = cell2mat(oldExtents_cell); % numeric array

signValues = sign(oldExtents(:,1));
offset = 0.6*signValues;

textPositions_cell = get(hText,{'Position'}); % cell array
textPositions = cell2mat(textPositions_cell); % numeric array
textPositions(:,1) = textPositions(:,1) + offset; % add offset

hText(1).Position = textPositions(1,:);
hText(2).Position = textPositions(2,:);
hText(3).Position = textPositions(3,:);
hText(4).Position = textPositions(4,:);

PrimSum = strcat({'Primary flux'},{' '},{num2str(round(PrimaryEnSum*10)/10)},...
            {' '},{'MW'});
t = title(PrimSum,'FontSize',12,'FontWeight','bold');
        pos = get(t,'position');
        set(t, 'position', pos+[0 0.15 0]);
 
   
set(handles.axes3,'CameraViewAngle',10.9074493116677,'DataAspectRatio',[1 1 1],...
    'PlotBoxAspectRatio',[1.2 1.2 1]);


        
        
        
        
