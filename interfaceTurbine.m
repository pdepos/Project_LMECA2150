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

% Last Modified by GUIDE v2.5 04-Dec-2016 15:29:45

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
    initialize_Recup(hObject, handles);
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


function NTUv_Callback(hObject, eventdata, handles)
% hObject    handle to NTUv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NTUv as text
%        str2double(get(hObject,'String')) returns contents of NTUv as a double
NTU = str2double(get(hObject, 'String'));
handles.metricdata.NTU = NTU;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function NTUv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NTUv (see GCBO)
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
    
    fuel = get(handles.fuel,'SelectedObject');
    if fuel == handles.gas
        x = 1;
        y = 4;
        z = 0;
        comb = 'CH4';
    elseif fuel == handles.diesel
        x = 12;
        y = 23;
        z = 0;
        comb = 'C12H23';
    end
    
    %%%% Fixed data %%%%
    etaPiC = 0.9;
    etaPiT = 0.9;
    
    %%%% Simulation %%%%
    [state,Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten,Exergy_losses,labels_Ex,etaRotex,etaCyclex,etaCombex,etaTotex,ma,mc,mg,lambda] = mainTurbineGaz(T1,r,etaPiC,kcc,T3,etaPiT,Pe,x,y,z,comb);
 
    %%%% plot Results %%%%
    
    % mass flow
    set(handles.mav,'String',ma);
    set(handles.mgv,'String',mg);
    set(handles.mcv,'String',mc);
    
    % efficiencies
    set(handles.mecenv,'String',etaMec);
    set(handles.cyclenv,'String',etaCyclen);
    set(handles.totenv,'String',etaToten);
    
    set(handles.mecexv,'String',etaMec);
    set(handles.rotexv,'String',etaRotex);
    set(handles.cyclexv,'String',etaCyclex);
    set(handles.combexv,'String',etaCombex);
    set(handles.totexv,'String',etaTotex);
    
    % states
    stateMat = [state{1}.T-273.15 state{1}.p state{1}.h state{1}.s state{1}.e;...
        state{2}.T-273.15 state{2}.p state{2}.h state{2}.s state{2}.e;...
        state{3}.T-273.15 state{3}.p state{3}.h state{3}.s state{3}.e;...
        state{4}.T-273.15 state{4}.p state{4}.h state{4}.s state{4}.e];
    set(handles.states,'Data',stateMat);
    
    % Diagrams
    graph1 = get(handles.choix1,'SelectedObject');
    if graph1 == handles.Benergy1
        axes(handles.graph1);
        En = pie(Energy_losses);
        hText = findobj(En,'Type','text');
        oldExtents_cell = get(hText,'Extent'); % cell array
        oldExtents = cell2mat(oldExtents_cell); % numeric array
        
        set(hText,{'String'},labels_Energy');
        
        newExtents_cell = get(hText,'Extent'); % cell array
        newExtents = cell2mat(newExtents_cell); % numeric array
        width_change = newExtents(:,3)-oldExtents(:,3);
        signValues = sign(oldExtents(:,1));
        offset = signValues.*(width_change/2);
        textPositions_cell = get(hText,{'Position'}); % cell array
        textPositions = cell2mat(textPositions_cell); % numeric array
        textPositions(:,1) = textPositions(:,1) + offset; % add offset
        set(hText,{'Position'},num2cell(textPositions,[3,2])); % set new position
        
        titre = strcat({'Primary flux'},{' '},{num2str(0.1*round(10*sum(Energy_losses)))},...
            {' '},{'MW'});
        t = title(titre,'FontSize',12,'FontWeight','bold');
        pos = get(t,'position');
        set(t, 'position', pos+[0 0.15 0]);
        
    elseif graph1 == handles.Bexergy1
        axes(handles.graph1);
        Ex = pie(Exergy_losses);
        hText = findobj(Ex,'Type','text');
        oldExtents_cell = get(hText,'Extent'); % cell array
        oldExtents = cell2mat(oldExtents_cell); % numeric array
        
        set(hText,{'String'},labels_Ex');
        
        newExtents_cell = get(hText,'Extent'); % cell array
        newExtents = cell2mat(newExtents_cell); % numeric array
        width_change = newExtents(:,3)-oldExtents(:,3);
        signValues = sign(oldExtents(:,1));
        offset = signValues.*(width_change/2);
        textPositions_cell = get(hText,{'Position'}); % cell array
        textPositions = cell2mat(textPositions_cell); % numeric array
        textPositions(:,1) = textPositions(:,1) + offset; % add offset
        set(hText,{'Position'},num2cell(textPositions,[3,2])) % set new position
        
        titre = strcat({'Primary flux'},{' '},{num2str(0.1*round(10*sum(Exergy_losses)))},...
            {' '},{'MW'});
        t = title(titre,'FontSize',12,'FontWeight','bold');
        pos = get(t,'position');
        set(t, 'position', pos+[0 0.15 0]);
    elseif graph1 == handles.BTS1
        
        
    elseif graph1 == handles.BHS1
        
    end
    
elseif get(handles.recuperator,'Value')==1
    
    %%%% Parameters %%%%
    T1 = handles.metricdata.T1; % [°C]
    r = handles.metricdata.r;
    Pe = handles.metricdata.Pe; % [MW]
    kcc = handles.metricdata.kcc;
    T3 = handles.metricdata.T3; % [°C]
    NTU = handles.metricdata.NTU;
    fuel = get(handles.fuel,'SelectedObject');
    if fuel == handles.gas
        x = 1;
        y = 4;
        z = 0;
        comb = 'CH4';
    elseif fuel == handles.diesel
        x = 12;
        y = 23;
        z = 0;
        comb = 'C12H23';
    end
    
    %%%% Fixed data %%%%
    etaPiC = 0.9;
    etaPiT = 0.9;
    
    %%%% Simulation %%%%
    [state,Energy_losses,labels_Energy,etaMec,etaCyclen,etaToten,Exergy_losses,labels_Ex,etaRotex,etaCyclex,etaCombex,etaTotex,ma,mc,mg,lambda] = mainTurbineGazRecup(T1,r,etaPiC,kcc,T3,etaPiT,Pe,x,y,z,NTU,comb);
 
    %%%% plot Results %%%%
    
    % mass flow
    set(handles.mav,'String',ma);
    set(handles.mgv,'String',mg);
    set(handles.mcv,'String',mc);
    
    % efficiencies
    set(handles.mecenv,'String',etaMec);
    set(handles.cyclenv,'String',etaCyclen);
    set(handles.totenv,'String',etaToten);
    
    set(handles.mecexv,'String',etaMec);
    set(handles.rotexv,'String',etaRotex);
    set(handles.cyclexv,'String',etaCyclex);
    set(handles.combexv,'String',etaCombex);
    set(handles.totexv,'String',etaTotex);
    
    % states
    stateMat = [state{1}.T-273.15 state{1}.p state{1}.h state{1}.s state{1}.e;...
        state{2}.T-273.15 state{2}.p state{2}.h state{2}.s state{2}.e;...
        state{3}.T-273.15 state{3}.p state{3}.h state{3}.s state{3}.e;...
        state{4}.T-273.15 state{4}.p state{4}.h state{4}.s state{4}.e;...
        state{5}.T-273.15 state{5}.p state{5}.h state{5}.s state{5}.e;...
        state{6}.T-273.15 state{6}.p state{6}.h state{6}.s state{6}.e];
    set(handles.states,'Data',stateMat);
    
    % Diagrams
    graph1 = get(handles.choix1,'SelectedObject');
    if graph1 == handles.Benergy1
        axes(handles.graph1);
        En = pie(Energy_losses);
        hText = findobj(En,'Type','text');
        oldExtents_cell = get(hText,'Extent'); % cell array
        oldExtents = cell2mat(oldExtents_cell); % numeric array
        
        set(hText,{'String'},labels_Energy');
        
        newExtents_cell = get(hText,'Extent'); % cell array
        newExtents = cell2mat(newExtents_cell); % numeric array
        width_change = newExtents(:,3)-oldExtents(:,3);
        signValues = sign(oldExtents(:,1));
        offset = signValues.*(width_change/2);
        textPositions_cell = get(hText,{'Position'}); % cell array
        textPositions = cell2mat(textPositions_cell); % numeric array
        textPositions(:,1) = textPositions(:,1) + offset; % add offset
        set(hText,{'Position'},num2cell(textPositions,[3,2])); % set new position
        
        titre = strcat({'Primary flux'},{' '},{num2str(0.1*round(10*sum(Energy_losses)))},...
            {' '},{'MW'});
        t = title(titre,'FontSize',12,'FontWeight','bold');
        pos = get(t,'position');
        set(t, 'position', pos+[0 0.15 0]);
        
    elseif graph1 == handles.Bexergy1
        axes(handles.graph1);
        Ex = pie(Exergy_losses);
        hText = findobj(Ex,'Type','text');
        oldExtents_cell = get(hText,'Extent'); % cell array
        oldExtents = cell2mat(oldExtents_cell); % numeric array
        
        set(hText,{'String'},labels_Ex');
        
        newExtents_cell = get(hText,'Extent'); % cell array
        newExtents = cell2mat(newExtents_cell); % numeric array
        width_change = newExtents(:,3)-oldExtents(:,3);
        signValues = sign(oldExtents(:,1));
        offset = signValues.*(width_change/2);
        textPositions_cell = get(hText,{'Position'}); % cell array
        textPositions = cell2mat(textPositions_cell); % numeric array
        textPositions(:,1) = textPositions(:,1) + offset; % add offset
        set(hText,{'Position'},num2cell(textPositions,[3,2])) % set new position
        
        titre = strcat({'Primary flux'},{' '},{num2str(0.1*round(10*sum(Exergy_losses)))},...
            {' '},{'MW'});
        t = title(titre,'FontSize',12,'FontWeight','bold');
        pos = get(t,'position');
        set(t, 'position', pos+[0 0.15 0]);
    elseif graph1 == handles.BTS1
        
    elseif graph1 == handles.BHS1
        
    end
end

function initialize_gui(fig_handle, handles)

% Cycle figure
axes(handles.cycle);
img = imread('TurbineGas.PNG');
image(img);
set(gca,'Visible','off');

% graph 1
axes(handles.graph1);
img2 = imread('epl.JPG');
image(img2);
set(gca,'Visible','off');

% default parameters
handles.metricdata.T1 = 15;
handles.metricdata.r = 18;
handles.metricdata.Pe = 230;
handles.metricdata.kcc = 0.95;
handles.metricdata.T3 = 1400;

set(handles.T1v,'String',handles.metricdata.T1);
set(handles.rv,'String',18);
set(handles.Pev,'String',230);
set(handles.kccv,'String',0.95);
set(handles.T3v,'String',1400);

% NTU
set(handles.NTU,'String',' ');
set(handles.NTUv,'Visible','off');

% mass flow
set(handles.mav,'String',0);
set(handles.mgv,'String',0);
set(handles.mcv,'String',0);
    
% efficiencies
set(handles.mecenv,'String',0);
set(handles.cyclenv,'String',0);
set(handles.totenv,'String',0);
    
set(handles.mecexv,'String',0);
set(handles.rotexv,'String',0);
set(handles.cyclexv,'String',0);
set(handles.combexv,'String',0);
set(handles.totexv,'String',0);

% state table
set(handles.states,'Data',zeros(4,5));
set(gca,'Visible','off');

% cla(handles.graph1);
guidata(handles.figure1, handles);

function initialize_Recup(fig_handle, handles)

% Cycle figure
axes(handles.cycle);
img = imread('TurbineGasRecup.PNG');
image(img);
set(gca,'Visible','off');

% graph 1
axes(handles.graph1);
img2 = imread('epl.jpg');
image(img2);
set(gca,'Visible','off');

% default parameters
handles.metricdata.T1 = 15;
handles.metricdata.r = 18;
handles.metricdata.Pe = 230;
handles.metricdata.kcc = 0.95;
handles.metricdata.T3 = 1400;

set(handles.T1v,'String',15);
set(handles.rv,'String',18);
set(handles.Pev,'String',230);
set(handles.kccv,'String',0.95);
set(handles.T3v,'String',1400);


% NTU
set(handles.NTU,'String','Number of transfer units (NTU)');
set(handles.NTUv,'Visible','on');

handles.metricdata.NTU = 2.5;
set(handles.NTUv,'String',2.5);
% mass flow
set(handles.mav,'String',0);
set(handles.mgv,'String',0);
set(handles.mcv,'String',0);
    
% efficiencies
set(handles.mecenv,'String',0);
set(handles.cyclenv,'String',0);
set(handles.totenv,'String',0);
    
set(handles.mecexv,'String',0);
set(handles.rotexv,'String',0);
set(handles.cyclexv,'String',0);
set(handles.combexv,'String',0);
set(handles.totexv,'String',0);

% state table
set(handles.states,'Data',zeros(6,5));

guidata(handles.figure1, handles);
