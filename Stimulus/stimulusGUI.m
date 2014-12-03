function varargout = stimulusGUI(varargin)
% STIMULUSGUI M-file for stimulusGUI.fig
%      STIMULUSGUI, by itself, creates a new STIMULUSGUI or raises the existing
%      singleton*.
%
%      H = STIMULUSGUI returns the handle to a new STIMULUSGUI or the handle to
%      the existing singleton*.
%
%      STIMULUSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STIMULUSGUI.M with the given input arguments.
%
%      STIMULUSGUI('Property','Value',...) creates a new STIMULUSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stimulusGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stimulusGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stimulusGUI

% Last Modified by GUIDE v2.5 26-Jan-2014 05:30:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @stimulusGUI_OpeningFcn, ...
    'gui_OutputFcn',  @stimulusGUI_OutputFcn, ...
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


% --- Executes just before stimulusGUI is made visible.
function stimulusGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stimulusGUI (see VARARGIN)

% Choose default command line output for stimulusGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stimulusGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stimulusGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% tspan=edit1_Callback(hObject, eventdata, handles);
% fs=edit2_Callback(hObject, eventdata, handles);
% wave=popupmenu1_Callback(hObject, eventdata, handles);
% freqs=edit3_Callback(hObject, eventdata, handles);
% amps=edit4_Callback(hObject, eventdata, handles);
% Ths=edit5_Callback(hObject, eventdata, handles);
% rTime=edit6_Callback(hObject, eventdata, handles);
% rExp=edit7_Callback(hObject, eventdata, handles);
% s=stimulusMake('fcn',tspan,fs,wave,freqs,amps,Ths,'ramp',rTime,rExp);
s=makeEverything(hObject,eventdata,handles);
range=edit9_Callback(hObject, eventdata, handles);
ind=round(s.fs*range(1))+1:round(s.fs*range(2));
figure(321);plot(s.t(ind),s.x(ind),'color','k');axis tight;zoom xon;drawnow;



% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% tspan=edit1_Callback(hObject, eventdata, handles);
% fs=edit2_Callback(hObject, eventdata, handles);
% wave=popupmenu1_Callback(hObject, eventdata, handles);
% freqs=edit3_Callback(hObject, eventdata, handles);
% amps=edit4_Callback(hObject, eventdata, handles);
% Ths=edit5_Callback(hObject, eventdata, handles);
% rTime=edit6_Callback(hObject, eventdata, handles);
% rExp=edit7_Callback(hObject, eventdata, handles);
perc=edit8_Callback(hObject, eventdata, handles);
% s=stimulusMake('fcn',tspan,fs,wave,freqs,amps,Ths,'ramp',rTime,rExp);
s=makeEverything(hObject,eventdata,handles);
allMyFreqs(s.x,8192,s.fs,perc);



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% tspan=mat2str(edit1_Callback(hObject, eventdata, handles));
% fs=num2str(edit2_Callback(hObject, eventdata, handles));
% wave=popupmenu1_Callback(hObject, eventdata, handles);
% wave=wave{:};
% freqs=mat2str(edit3_Callback(hObject, eventdata, handles));
% amps=mat2str(edit4_Callback(hObject, eventdata, handles));
% Ths=mat2str(edit5_Callback(hObject, eventdata, handles));
% rTime=mat2str(edit6_Callback(hObject, eventdata, handles));
% rExp=mat2str(edit7_Callback(hObject, eventdata, handles));
str=getEverything(hObject, eventdata, handles);
% str=['s=stimulusMake(''fcn''' ',' tspan ',' fs ',' '{''' wave '''}' ','...
%     freqs ',' amps ',' Ths ',' '''ramp'',' rTime ',' rExp ')'];
callStr=edit10_Callback(hObject, eventdata, handles);
disp(str);
assignin('base',callStr,str);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% tspan=edit1_Callback(hObject, eventdata, handles);
% fs=edit2_Callback(hObject, eventdata, handles);
% wave=popupmenu1_Callback(hObject, eventdata, handles);
% freqs=edit3_Callback(hObject, eventdata, handles);
% amps=edit4_Callback(hObject, eventdata, handles);
% Ths=edit5_Callback(hObject, eventdata, handles);
% rTime=edit6_Callback(hObject, eventdata, handles);
% rExp=edit7_Callback(hObject, eventdata, handles);
% s=stimulusMake('fcn',tspan,fs,wave,freqs,amps,Ths,'ramp',rTime,rExp);
s=makeEverything(hObject,eventdata,handles);
structStr=edit11_Callback(hObject, eventdata, handles);
assignin('base',structStr,s);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
global playingCurrently
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% tspan=edit1_Callback(hObject, eventdata, handles);
% fs=edit2_Callback(hObject, eventdata, handles);
% wave=popupmenu1_Callback(hObject, eventdata, handles);
% freqs=edit3_Callback(hObject, eventdata, handles);
% amps=edit4_Callback(hObject, eventdata, handles);
% Ths=edit5_Callback(hObject, eventdata, handles);
% rTime=edit6_Callback(hObject, eventdata, handles);
% rExp=edit7_Callback(hObject, eventdata, handles);
% s=stimulusMake('fcn',tspan,fs,wave,freqs,amps,Ths,'ramp',rTime,rExp);
s=makeEverything(hObject,eventdata,handles);
% soundsc(s.x,s.fs);
playingCurrently=audioplayer(real(s.x*.999/max(abs(s.x))),s.fs);
play(playingCurrently);




% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
global playingCurrently
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% clear playsnd
stop(playingCurrently);





function tspan=edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
temp=sortMeOut(get(handles.edit1,'String'));
if isscalar(temp.val)
    tspan.val=[0 temp.val];
    tspan.str=['[' '0' ' ' temp.str ']'];
else
    tspan.val=temp.val;
    tspan.str=temp.str;
end


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fs=edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
temp=get(handles.edit2,'String');
fs=sortMeOut(temp);
if numel(fs.val)~=1
    error('Sampling frequency must be a scalar');
end


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function wave=popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
val=get(handles.popupmenu1,'Value');
switch val
    case 1
        wave={'sin'};
    case 2
        wave={'cos'};
    case 3
        wave={'exp'};
    case 4
        wave={'saw'};
    case 5
        wave={'squ'};
    case 6
        wave={'bmp'};
    case 7
        wave={'noi'};
end



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freqs=edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
temp=get(handles.edit3,'String');
freqs=sortMeOut(temp);



% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function amps=edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
temp=get(handles.edit4,'String');
amps=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ths=edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
temp=get(handles.edit5,'String');
Ths=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rTime=edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
temp=get(handles.edit6,'String');
rTime=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rExp=edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
temp=get(handles.edit7,'String');
rExp=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function freqPerc=edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
temp=str2num(get(handles.edit8,'String'));
if isempty(temp)
    freqPerc=[0 100];
else
    if numel(temp)~=2
        error('Range of freqs must be two numbers');
    end
    fs=edit2_Callback(hObject, eventdata, handles);
    nyq=fs.val/2;
    if max(temp)>nyq
        error('Specified frequency range to display is beyond the Nyquist frequency');
    end
    freqPerc=[100*temp(1)/nyq 100*temp(2)/nyq];
end


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function timeRange=edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
temp=str2num(get(handles.edit9,'String'));
tspan=edit1_Callback(hObject, eventdata, handles);
if isempty(temp)
    timeRange=[min(min(tspan.val)) max(max(tspan.val))];
else
    if numel(temp)~=2
        error('Range of time must be two numbers');
    end
    if max(temp)>max(max(tspan.val))
        error('Specified time range to display is outside bounds of stimulus');
    end
    timeRange=temp;
end


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function callStr=edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
callStr=get(handles.edit10,'String');


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function structStr=edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
structStr=get(handles.edit11,'String');


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu2.
function wave=popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
val=get(handles.popupmenu2,'Value');
switch val
    case 1
        wave=[];
    case 2
        wave={'sin'};
    case 3
        wave={'cos'};
    case 4
        wave={'saw'};
    case 5
        wave={'squ'};
end


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fAM=edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
temp=get(handles.edit12,'String');
fAM=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function aAM=edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
temp=get(handles.edit13,'String');
aAM=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function wave=popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
val=get(handles.popupmenu3,'Value');
switch val
    case 1
        wave=[];
    case 2
        wave={'cos'};
    case 3
        wave={'saw'};
    case 4
        wave={'squ'};
end


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fFM=edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
temp=get(handles.edit14,'String');
fFM=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function aFM=edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
temp=get(handles.edit15,'String');
aFM=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iterDelay=edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
temp=get(handles.edit16,'String');
iterDelay=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iterNum=edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double
temp=get(handles.edit17,'String');
iterNum=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maskSNR=edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double
temp=get(handles.edit18,'String');
maskSNR=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filtstim=edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double
temp=get(handles.edit19,'String');
filtstim=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filtmask=edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
temp=get(handles.edit20,'String');
filtmask=sortMeOut(temp);


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function s=makeEverything(hObject,eventdata,handles)

temp=edit1_Callback(hObject, eventdata, handles);
tspan=temp.val;

temp=edit2_Callback(hObject, eventdata, handles);
fs=temp.val;

wave=popupmenu1_Callback(hObject, eventdata, handles);
% wave=wave{:};

temp=edit3_Callback(hObject, eventdata, handles);
freqs=temp.val;

temp=edit4_Callback(hObject, eventdata, handles);
amps=temp.val;

temp=edit5_Callback(hObject, eventdata, handles);
Ths=temp.val;

if isempty(Ths)
    Ths=0;
end

temp=edit6_Callback(hObject, eventdata, handles);
rTime=temp.val;

if isempty(rTime)
    rTime=0;
end    

temp=edit7_Callback(hObject, eventdata, handles);
rExp=temp.val;

if isempty(rExp)
    rExp=0;
end    

waveAM=popupmenu2_Callback(hObject, eventdata, handles);
% waveAM=waveAM{:};

temp=edit12_Callback(hObject, eventdata, handles);
fAM=temp.val;

if isempty(fAM)
    fAM=0;
end

temp=edit13_Callback(hObject, eventdata, handles);
aAM=temp.val;

if isempty(aAM)
    aAM=0;
end

if isempty(waveAM)
    waveAM={'sin'};
    fAM=0;
    aAM=0;
end

waveFM=popupmenu3_Callback(hObject, eventdata, handles);
% waveFM=waveFM{:};

temp=edit14_Callback(hObject, eventdata, handles);
fFM=temp.val;

if isempty(fFM)
    fFM=0;
end

temp=edit15_Callback(hObject, eventdata, handles);
aFM=temp.val;

if isempty(aFM)
    aFM=0;
end

if isempty(waveFM)
    waveFM={'cos'};
    fFM=0;
    aFM=0;
end

temp=edit16_Callback(hObject, eventdata, handles);
iterDelay=temp.val;

if isempty(iterDelay)
    iterDelay=0;
end

temp=edit17_Callback(hObject, eventdata, handles);
iterNum=temp.val;

if isempty(iterNum)
    iterNum=0;
end

temp=edit18_Callback(hObject, eventdata, handles);
maskSNR=temp.val;

if isempty(maskSNR)
    maskSNR=inf;
end

temp=edit19_Callback(hObject, eventdata, handles);
filtstim=temp.val;

if isempty(filtstim)
    filtstim={[] []};
end

temp=edit20_Callback(hObject, eventdata, handles);
filtmask=temp.val;

if isempty(filtmask)
    filtmask={[] []};
end

s=stimulusMake('fcn',tspan,fs,wave,freqs,amps,Ths,'ramp',rTime,rExp,...
    'am',waveAM,fAM,aAM,1,'fm',waveFM,fFM,aFM,'iter',iterDelay,iterNum,...
    'mask',maskSNR,'filtstim',filtstim,'filtmask',filtmask);


function str=getEverything(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

temp=edit1_Callback(hObject, eventdata, handles);
tspan=temp.str;

temp=edit2_Callback(hObject, eventdata, handles);
fs=temp.str;

wave=popupmenu1_Callback(hObject, eventdata, handles);
wave=wave{:};

temp=edit3_Callback(hObject, eventdata, handles);
freqs=temp.str;

temp=edit4_Callback(hObject, eventdata, handles);
amps=temp.str;


str=['s=stimulusMake(''fcn''' ',' tspan ',' fs ',' '{''' wave '''}' ','...
    freqs ',' amps];


temp=edit5_Callback(hObject, eventdata, handles);
Ths=temp.str;

if ~isempty(Ths)
    str=[str ',' Ths];
end

temp=edit6_Callback(hObject, eventdata, handles);
rTime=temp.str;

temp=edit7_Callback(hObject, eventdata, handles);
rExp=temp.str;

if ~isempty(rTime) && ~isempty(rExp)
    str=[str ',' '''ramp'',' rTime ',' rExp];
end

waveAM=popupmenu2_Callback(hObject, eventdata, handles);

temp=edit12_Callback(hObject, eventdata, handles);
fAM=temp.str;

temp=edit13_Callback(hObject, eventdata, handles);
aAM=temp.str;


if ~isempty(waveAM) && ~isempty(aAM) && ~isempty(fAM)
    waveAM=waveAM{:};
    str=[str ',' '''am'',' '{''' waveAM '''}' ',' fAM ',' aAM ',' '1'];
end


waveFM=popupmenu3_Callback(hObject, eventdata, handles);

temp=edit14_Callback(hObject, eventdata, handles);
fFM=temp.str;

temp=edit15_Callback(hObject, eventdata, handles);
aFM=temp.str;


if ~isempty(waveFM) && ~isempty(aFM) && ~isempty(fFM)
    waveFM=waveFM{:};
    str=[str ',' '''fm'',' '{''' waveFM '''}' ',' fFM ',' aFM];
end

temp=edit16_Callback(hObject, eventdata, handles);
iterDelay=temp.str;

temp=edit17_Callback(hObject, eventdata, handles);
iterNum=temp.str;


if ~isempty(iterDelay) && ~isempty(iterNum)
    str=[str ',' '''iter'',' iterDelay ',' iterNum];
end

temp=edit18_Callback(hObject, eventdata, handles);
maskSNR=temp.str;

if ~isempty(maskSNR)
    str=[str ',' '''mask'',' maskSNR];
end

temp=edit19_Callback(hObject, eventdata, handles);
filtstim=temp.str;

if ~isempty(filtstim)
    str=[str ',' '''filtstim'',' filtstim];
end

temp=edit20_Callback(hObject, eventdata, handles);
filtmask=temp.str;

if ~isempty(filtmask)
    str=[str ',' '''filtmask'',' filtmask];
end

str=[str ')'];



function struct=sortMeOut(temp)

if ~isempty(temp)
    struct.val=str2num(temp); %#ok<*ST2NM>
    if isempty(struct.val)
        if strcmpi(temp(1),'{')
            struct.val=eval(temp);
        else            
            try
                struct.val=evalin('base',temp);
            catch err
                error('Variable "%s" does not exist in workspace',temp);
            end
        end
        struct.str=temp;
    else
        struct.str=mat2str(struct.val);
    end
else
    struct.str=[];
    struct.val=[];
end
