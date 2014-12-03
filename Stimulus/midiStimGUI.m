function varargout = midiStimGUI(varargin)
% MIDISTIMGUI MATLAB code for midiStimGUI.fig
%      MIDISTIMGUI, by itself, creates a new MIDISTIMGUI or raises the existing
%      singleton*.
%
%      H = MIDISTIMGUI returns the handle to a new MIDISTIMGUI or the handle to
%      the existing singleton*.
%
%      MIDISTIMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MIDISTIMGUI.M with the given input arguments.
%
%      MIDISTIMGUI('Property','Value',...) creates a new MIDISTIMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before midiStimGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to midiStimGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help midiStimGUI

% Last Modified by GUIDE v2.5 07-Sep-2014 20:44:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @midiStimGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @midiStimGUI_OutputFcn, ...
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


% --- Executes just before midiStimGUI is made visible.
function midiStimGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to midiStimGUI (see VARARGIN)

% Choose default command line output for midiStimGUI
handles.output = hObject;

handles.file_name = [];
handles.vararg.time_span = [];
handles.vararg.sample_rate = [];
handles.vararg.synth_type_opt = [];
handles.vararg.num_chan_opt = [];
handles.vararg.note_range = [];
handles.vararg.tempo_mod_opt = [];
handles.vararg.mod_func = [];
handles.vararg.mod_depth = [];
handles.vararg.beat_per = [];
handles.vararg.mod_per  = [];

handles.vararg.synth_type_opt = [];

handles.stim_wrkspc_name = 'stim_str';
handles.stim_samples = [];
%even need the following?
%handles.stim_sample_rate = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes midiStimGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = midiStimGUI_OutputFcn(hObject, eventdata, handles) 
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
callStimMake(hObject, handles);


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
handles.stim_wrkspc_name = get(hObject,'String');
guidata(hObject, handles);



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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
switch get(hObject,'Value')
    case 1
        %default
        if ~isempty(handles.vararg.num_chan_opt)
            handles.vararg.num_chan_opt = [];
        end
        set(handles.edit4, 'String', '')
        set(handles.edit4, 'Enable', 'off')
        handles.vararg.note_range = [];
    case 2
        handles.vararg.num_chan_opt = {'chan_per_note'};
        set(handles.edit4, 'Enable', 'on')
    case 3
        handles.vararg.num_chan_opt = {'chan_per_chan'};
        set(handles.edit4, 'String', '')
        set(handles.edit4, 'Enable', 'off')
        handles.vararg.note_range = [];
end
guidata(hObject, handles);


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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
handles.vararg.note_range = str2num(get(hObject,'String'));
guidata(hObject, handles);


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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if get(hObject, 'Value') == 1 %checked
    if isempty(handles.vararg.beat_per) %it can only be empty if all are empty
        set(handles.edit5, 'String', '.5') %bper
        set(handles.edit6, 'String', '.2') %mdepth
        set(handles.edit7, 'String', '16') %mper
    end
    set(handles.edit5, 'Enable', 'on')
    set(handles.edit6, 'Enable', 'on')
    set(handles.edit7, 'Enable', 'on')
    set(handles.popupmenu2, 'Enable', 'on')
    handles.vararg.tempo_mod_opt = {'tempo_mod'};
else
    set(handles.edit5, 'Enable', 'off')
    set(handles.edit6, 'Enable', 'off')
    set(handles.edit7, 'Enable', 'off')
    set(handles.popupmenu2, 'Enable', 'off')
    if ~isempty(handles.vararg.tempo_mod_opt)
        handles.vararg.tempo_mod_opt = [];
    end
end
guidata(hObject, handles);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.vararg.time_span = str2num(get(hObject,'String'));
guidata(hObject, handles);


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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.vararg.sample_rate = str2num(get(hObject,'String'));
guidata(hObject, handles);


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



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.file_name = get(hObject,'String');
guidata(hObject, handles);


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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
handles.vararg.beat_per = {str2num(get(hObject,'String'))};
guidata(hObject, handles);


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



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
handles.vararg.mod_depth = {str2num(get(hObject,'String'))};

if isempty(handles.vararg.beat_per)
    handles.vararg.beat_per = {str2num(get(handles.edit5, 'String'))};
end
guidata(hObject, handles);


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



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
handles.vararg.mod_per = str2num(get(hObject,'String'));

if isempty(handles.vararg.mod_func)
    handles.vararg.mod_func = {getFuncName(get(handles.popupmenu2,'Value'))};
    handles.vararg.mod_depth = {str2num(get(handles.edit6, 'String'))};
    handles.vararg.beat_per = {str2num(get(handles.edit5, 'String'))};
end

guidata(hObject, handles);


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


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
handles.vararg.mod_func = {getFuncName(get(hObject,'Value'))};

if isempty(handles.vararg.mod_depth)
    handles.vararg.mod_depth = {str2num(get(handles.edit6, 'String'))};
    handles.vararg.beat_per = {str2num(get(handles.edit5, 'String'))};
end

guidata(hObject, handles);

function fname = getFuncName(val)
switch val
    case 1
        fname = 'sin';
    case 2
        fname = 'cos';
    case 3
        fname = 'linear';
    case 4
        fname = 'sawtooth';
    case 5
        fname = 'square';
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


% --------------------------------------------------------------------
% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue, 'Tag')
    case 'radiobutton1'
        if ~isempty(handles.vararg.synth_type_opt)
            handles.vararg.synth_type_opt = [];
            set(handles.edit3, 'String', '160')
            handles.vararg.sample_rate = []; %go back to default
            set(handles.pushbutton2, 'Enable', 'off')
            set(handles.pushbutton3, 'Enable', 'off')
        end
    case 'radiobutton2'
        if isempty(handles.vararg.synth_type_opt)
            handles.vararg.synth_type_opt = {'melody'};
            set(handles.edit3, 'String', '')
            handles.vararg.sample_rate = []; %in case has a value
            set(handles.pushbutton2, 'Enable', 'on')
            set(handles.pushbutton3, 'Enable', 'on')
        end
end
guidata(hObject, handles);


function callStimMake(hObject, handles)
    if isempty(handles.file_name)
        error('Enter a valid MIDI file path')
    end
    
    if ~isempty(handles.vararg.time_span)
       init_vararg = [{handles.vararg.time_span} handles.vararg.sample_rate];
    else
        init_vararg = [];
    end
    
    tempo_args = [];
    if ~isempty(handles.vararg.tempo_mod_opt)
        if ~isempty(handles.vararg.beat_per)
            tempo_args = [handles.vararg.beat_per, handles.vararg.mod_depth, handles.vararg.mod_func, handles.vararg.mod_per]; 
        end
    end
    
    m_vararg = [init_vararg, handles.vararg.num_chan_opt, handles.vararg.note_range, handles.vararg.synth_type_opt, handles.vararg.tempo_mod_opt];
    if ~isempty(tempo_args)
        m_vararg{end+1} = tempo_args;
    end
    
    if ~isempty(m_vararg)
        s = stimulusMake('mid', handles.file_name, m_vararg{:});
    else
        s = stimulusMake('mid', handles.file_name);
    end
    
    handles.stim_samples = s.x;
    handles.stim_sample_rate = s.fs;
    
    assignin('base',handles.stim_wrkspc_name,s);
    guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.stim_samples)
    handles.audio_samples=audioplayer(real(handles.stim_samples*.999/max(abs(handles.stim_samples))),handles.stim_sample_rate);
    guidata(hObject, handles);
    play(handles.audio_samples);
else
    error('Must create a stimulus to play');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
stop(handles.audio_samples);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName, ~] = uigetfile('*.mid');
set(handles.edit1, 'String', horzcat(PathName, FileName))
handles.file_name = horzcat(PathName, FileName);
guidata(hObject, handles);