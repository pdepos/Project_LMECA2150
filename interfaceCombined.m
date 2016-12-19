function varargout = interfaceCombined(varargin)
% INTERFACECOMBINED MATLAB code for interfaceCombined.fig
%      INTERFACECOMBINED, by itself, creates a new INTERFACECOMBINED or raises the existing
%      singleton*.
%
%      H = INTERFACECOMBINED returns the handle to a new INTERFACECOMBINED or the handle to
%      the existing singleton*.
%
%      INTERFACECOMBINED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFACECOMBINED.M with the given input arguments.
%
%      INTERFACECOMBINED('Property','Value',...) creates a new INTERFACECOMBINED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before interfaceCombined_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to interfaceCombined_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help interfaceCombined

% Last Modified by GUIDE v2.5 07-Dec-2016 22:33:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @interfaceCombined_OpeningFcn, ...
                   'gui_OutputFcn',  @interfaceCombined_OutputFcn, ...
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


% --- Executes just before interfaceCombined is made visible.
function interfaceCombined_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to interfaceCombined (see VARARGIN)

% Choose default command line output for interfaceCombined
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
initialize_gui(hObject, handles);

% UIWAIT makes interfaceCombined wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = interfaceCombined_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function T1gv_Callback(hObject, eventdata, handles)
% hObject    handle to T1gv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T1gv as text
%        str2double(get(hObject,'String')) returns contents of T1gv as a double
T1g = str2double(get(hObject, 'String'));
handles.metricdata.T1g = T1g;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function T1gv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T1gv (see GCBO)
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



function Pegv_Callback(hObject, eventdata, handles)
% hObject    handle to Pegv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pegv as text
%        str2double(get(hObject,'String')) returns contents of Pegv as a double
Peg = str2double(get(hObject, 'String'));
handles.metricdata.Peg = Peg;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Pegv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pegv (see GCBO)
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



function T3gv_Callback(hObject, eventdata, handles)
% hObject    handle to T3gv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T3gv as text
%        str2double(get(hObject,'String')) returns contents of T3gv as a double
T3g = str2double(get(hObject, 'String'));
handles.metricdata.T3g = T3g;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function T3gv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T3gv (see GCBO)
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

%%%% record data %%%%
T1g = handles.metricdata.T1g;
r = handles.metricdata.r;
Peg = handles.metricdata.Peg;
kcc = handles.metricdata.kcc;
T3g = handles.metricdata.T3g;
DTa = handles.metricdata.DTa;
DTlp = handles.metricdata.DTlp;
DThp = handles.metricdata.DThp;
pLP = handles.metricdata.pLP;
pHP = handles.metricdata.pHP;
Tw = handles.metricdata.Tw;

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
 
%%%% simulation %%%%

[stateV,stateTG,Energy_lossesTGV,labels_EnergyTGV,Exergy_lossesTGV,labels_ExTGV,massFlow,energyEff,exergyEff] = ...
    mainCombined2(Peg,pHP,pLP,Tw,DTa,DTlp,DThp,x,y,z,T1g,r,kcc,T3g,comb);

%%%% plot data %%%%

% Gas states
stateTGmat = zeros(5,5);
for i =1:5
   stateTGmat(i,:) = [stateTG{i}.T stateTG{i}.p stateTG{i}.h stateTG{i}.s stateTG{i}.e];  
end
set(handles.TabTG,'Data',stateTGmat);

% Steam states
stateVmat = cell(10,6);
%stateVmat = zeros(10,6);
for i = 1:10
   %stateVmat{i,:} = [stateV{i}.T stateV{i}.p stateV{i}.h stateV{i}.s stateV{i}.e stateV{i}.x];
   stateVmat{i,1} = stateV{i}.T;
   stateVmat{i,2} = stateV{i}.p;
   stateVmat{i,3} = stateV{i}.h;
   stateVmat{i,4} = stateV{i}.s;
   stateVmat{i,5} = stateV{i}.e;
   if isfinite(stateV{i}.x)
    stateVmat{i,6} = stateV{i}.x;
   else
       stateVmat{i,6} = '                     -';
   end
end
set(handles.TabV,'Data',stateVmat);

% Mass flows
ma = massFlow(1);
mg = massFlow(3);
mc = massFlow(2);
mvLP = massFlow(5);
mvHP = massFlow(4);

set(handles.mav,'String',ma);
set(handles.mcv,'String',mc);
set(handles.mgv,'String',mg);
set(handles.mHPv,'String',mvHP);
set(handles.mLPv,'String',mvLP);
set(handles.mTotv,'String',mvLP+mvHP);

% Energy analysis
etaMecV = energyEff(1);
etaCyclenV = energyEff(2);
etaTotenV = energyEff(3);
etaTotenTGV = energyEff(4);
etaMecTG = energyEff(5);
etaCyclenTG = energyEff(6);
etaTotenTG = energyEff(7);

set(handles.etaMecenVv,'String',etaMecV);
set(handles.etaCyclenVv,'String',etaCyclenV);
set(handles.etaTotenVv,'String',etaTotenV);
set(handles.etaTotenTGVv,'String',etaTotenTGV);
set(handles.etaMecenTGv,'String',etaMecTG);
set(handles.etaCyclenTGv,'String',etaCyclenTG);
set(handles.etaTotenTGv,'String',etaTotenTG);

% Exergy analysis

etaMecTG = exergyEff(1);
etaRotexTG = exergyEff(2);
etaCyclexTG = exergyEff(3);
etaCombexTG = exergyEff(4);
etaTotexTG = exergyEff(5);
etaRotexV = exergyEff(6);
etaCyclexV = exergyEff(7);
etaTransexV = exergyEff(8);
etaTotexV = exergyEff(9);
etaTotexTGV = exergyEff(10);

set(handles.etaMecexTGv,'String',etaMecTG);
set(handles.etaRotexTGv,'String',etaRotexTG);
set(handles.etaCyclexTGv,'String',etaCyclexTG);
set(handles.etaCombexTGv,'String',etaCombexTG);
set(handles.etaTotexTGv,'String',etaTotexTG);
set(handles.etaTotenTGv,'String',etaTotenTG);
set(handles.etaMecexVv,'String',etaMecV);
set(handles.etaRotexVv,'String',etaRotexV);
set(handles.etaCyclexVv,'String',etaCyclexV);
set(handles.etaTransexVv,'String',etaTransexV);
set(handles.etaTotexVv,'String',etaTotexV);
set(handles.etaTotexTGVv,'String',etaTotexTGV);

% Graph 1
graph1 = get(handles.choix1,'SelectedObject');
    if graph1 == handles.energy1
        axes(handles.graph1);
        En = pie(Energy_lossesTGV);
        hText = findobj(En,'Type','text');
        oldExtents_cell = get(hText,'Extent'); % cell array
        oldExtents = cell2mat(oldExtents_cell); % numeric array
        
        set(hText,{'String'},labels_EnergyTGV');
        
        newExtents_cell = get(hText,'Extent'); % cell array
        newExtents = cell2mat(newExtents_cell); % numeric array
        width_change = newExtents(:,3)-oldExtents(:,3);
        signValues = sign(oldExtents(:,1));
        offset = signValues.*(width_change/2);
        textPositions_cell = get(hText,{'Position'}); % cell array
        textPositions = cell2mat(textPositions_cell); % numeric array
        textPositions(:,1) = textPositions(:,1) + offset; % add offset
        set(hText,{'Position'},num2cell(textPositions,[3,2])); % set new position
        
        titre = strcat({'Primary flux'},{' '},{num2str(0.1*round(10*sum(Energy_lossesTGV)))},...
            {' '},{'MW'});
        t = title(titre,'FontSize',12,'FontWeight','bold');
        pos = get(t,'position');
        set(t, 'position', pos+[0 0.01 0]);
        
    elseif graph1 == handles.exergy1
        axes(handles.graph1);
        Ex = pie(Exergy_lossesTGV);
        hText = findobj(Ex,'Type','text');
        oldExtents_cell = get(hText,'Extent'); % cell array
        oldExtents = cell2mat(oldExtents_cell); % numeric array
        
        set(hText,{'String'},labels_ExTGV');
        
        newExtents_cell = get(hText,'Extent'); % cell array
        newExtents = cell2mat(newExtents_cell); % numeric array
        width_change = newExtents(:,3)-oldExtents(:,3);
        signValues = sign(oldExtents(:,1));
        offset = signValues.*(width_change/2);
        textPositions_cell = get(hText,{'Position'}); % cell array
        textPositions = cell2mat(textPositions_cell); % numeric array
        textPositions(:,1) = textPositions(:,1) + offset; % add offset
        set(hText,{'Position'},num2cell(textPositions,[3,2])) % set new position
        
        titre = strcat({'Primary flux'},{' '},{num2str(0.1*round(10*sum(Exergy_lossesTGV)))},...
            {' '},{'MW'});
        t = title(titre,'FontSize',12,'FontWeight','bold');
        pos = get(t,'position');
        set(t, 'position', pos+[0 0.15 0]);
        
    elseif graph1 == handles.TS1
        
        [sL,sV,TL,TV,s110,T110,s39,T39] = combined2TS(stateV);
        
        n = length(stateV);
        stateT = zeros(n,1);
        stateS = zeros(n,1);
        for j = 1:n
            stateT(j) = stateV{j}.T;
            stateS(j) = stateV{j}.s;
        end

        axes(handles.graph1);
        plot(sL,TL,'b',sV,TV,'b',s110,T110,'g',s39,T39,'g',stateS,stateT,'r*');
        xlabel('s [kJ/(kg*K)]');
        ylabel('T [°C]');
        
    elseif graph1 == handles.hs1
        
        [sL,sV,hL,hV,s110,h110,s39,h39] = combined2hs(stateV);
        
        n = length(stateV);
        stateh = zeros(n,1);
        stateS = zeros(n,1);
        for j = 1:n
            stateh(j) = stateV{j}.h;
            stateS(j) = stateV{j}.s;
        end
        
        axes(handles.graph1);
        plot(sL,hL,'b',sV,hV,'b',s110,h110,'g',s39,h39,'g',stateS,stateh,'r*');
        xlabel('s [kJ/(kg*K)]');
        ylabel('h [kJ/kg]');
        
    end

function pHPv_Callback(hObject, eventdata, handles)
% hObject    handle to pHPv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pHPv as text
%        str2double(get(hObject,'String')) returns contents of pHPv as a double
pHP = str2double(get(hObject, 'String'));
handles.metricdata.pHP = pHP;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function pHPv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pHPv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pLPv_Callback(hObject, eventdata, handles)
% hObject    handle to pLPv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pLPv as text
%        str2double(get(hObject,'String')) returns contents of pLPv as a double
pLP = str2double(get(hObject, 'String'));
handles.metricdata.pLP = pLP;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function pLPv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pLPv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Twv_Callback(hObject, eventdata, handles)
% hObject    handle to Twv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Twv as text
%        str2double(get(hObject,'String')) returns contents of Twv as a double
Tw = str2double(get(hObject, 'String'));
handles.metricdata.Tw = Tw;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Twv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Twv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTav_Callback(hObject, eventdata, handles)
% hObject    handle to DTav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTav as text
%        str2double(get(hObject,'String')) returns contents of DTav as a double
DTa = str2double(get(hObject, 'String'));
handles.metricdata.DTa = DTa;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DTav_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTav (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DTlpv_Callback(hObject, eventdata, handles)
% hObject    handle to DTlpv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DTlpv as text
%        str2double(get(hObject,'String')) returns contents of DTlpv as a double
DTlp = str2double(get(hObject, 'String'));
handles.metricdata.DTlp = DTlp;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DTlpv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DTlpv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DThpv_Callback(hObject, eventdata, handles)
% hObject    handle to DThpv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DThpv as text
%        str2double(get(hObject,'String')) returns contents of DThpv as a double
DThp = str2double(get(hObject, 'String'));
handles.metricdata.DThp = DThp;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function DThpv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DThpv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function initialize_gui(fig_handle, handles)
 
% default parameters
handles.metricdata.T1g = 15;
handles.metricdata.r = 18;
handles.metricdata.Peg = 230;
handles.metricdata.kcc = 0.95;
handles.metricdata.T3g = 1400;
handles.metricdata.pLP = 5.8;
handles.metricdata.pHP = 78;
handles.metricdata.Tw = 15;
handles.metricdata.DTa = 50;
handles.metricdata.DTlp = 10;
handles.metricdata.DThp = 10;

set(handles.T1gv,'String',handles.metricdata.T1g);
set(handles.rv,'String',handles.metricdata.r );
set(handles.Pegv,'String',handles.metricdata.Peg);
set(handles.kccv,'String',handles.metricdata.kcc);
set(handles.T3gv,'String',handles.metricdata.T3g);
set(handles.pLPv,'String',handles.metricdata.pLP);
set(handles.pHPv,'String',handles.metricdata.pHP);
set(handles.Twv,'String',handles.metricdata.Tw);
set(handles.DTav,'String',handles.metricdata.DTa);
set(handles.DTlpv,'String',handles.metricdata.DTlp);
set(handles.DThpv,'String',handles.metricdata.DThp);

% mass flow
set(handles.mav,'String',0);
set(handles.mgv,'String',0);
set(handles.mcv,'String',0);
set(handles.mHPv,'String',0);
set(handles.mLPv,'String',0);
set(handles.mTotv,'String',0);

% efficiencies
set(handles.etaMecenTGv,'String',0);
set(handles.etaCyclenTGv,'String',0);
set(handles.etaTotenTGv,'String',0);
    
set(handles.etaMecexTGv,'String',0);
set(handles.etaRotexTGv,'String',0);
set(handles.etaCyclexTGv,'String',0);
set(handles.etaCombexTGv,'String',0);
set(handles.etaTotexTGv,'String',0);

set(handles.etaMecenVv,'String',0);
set(handles.etaCyclenVv,'String',0);
set(handles.etaTotenVv,'String',0);

set(handles.etaTotenTGVv,'String',0);

set(handles.etaMecexVv,'String',0);
set(handles.etaRotexVv,'String',0);
set(handles.etaCyclexVv,'String',0);
set(handles.etaTransexVv,'String',0);
set(handles.etaTotexVv,'String',0);

set(handles.etaTotexTGVv,'String',0);

% state table
set(handles.TabTG,'Data',zeros(5,5));
set(handles.TabV,'Data',zeros(10,6));

set(gca,'Visible','off');
cla(handles.graph1);

guidata(handles.figure1, handles);
