function varargout = LFM_popup(varargin)
%
% PROJECT2GUI_popup M-file for LFM_popup.fig used in conjunction
% with PROJECT2GUI
%
%Written by Gregory Heim as part of a graduate special problems course for
%Dr. Mark Richards.
% 
% This work is licensed under the Creative Commons Attribution-ShareAlike
% License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-sa/2.5/ 
% or send a letter to:
% Creative Commons,
% 543 Howard Street, 5th Floor, 
% San Francisco, California, 94105, USA.
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LFM_popup_OpeningFcn, ...
                   'gui_OutputFcn',  @LFM_popup_OutputFcn, ...
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


% --- Executes just before LFM_popup is made visible.
function LFM_popup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LFM_popup (see VARARGIN)

% Choose default command line output for LFM_popup
handles.output = hObject;
%Pass the handle of the calling figure to this one and check
%wether it is the figure or just the structure

handles.dist=1000;%default distance
handles.plot_check=0;%sets the output plot to log value
handles.normal_check=0; %normalizes the output plot 
handles.time_scale=100; %adjusts the amount of time shown on the output plot
handles.max_db=40; %maximum db on the plot
handles.min_db=-10; %minimum db on the plot

guidata(hObject, handles);
time_text_CreateFcn(hObject, eventdata, handles)
new_plots(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LFM_popup wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LFM_popup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%plots the matched filter output when there are two scatterers
function new_plots(hObject, eventdata, handles)

%get global variables from Project2GUI
global TT
global XX
global XS
global HH
T=TT;
x = XX;
xs = XS;
h=HH;

%calculate time difference between scatteres and offsetters in samples
 lag=handles.dist/3e8;
 offset = round(lag/T);
 
 %create input if there were two scatters
 rc=zeros(1,length(x)+offset);
 rc(1:length(x)) =rc(1:length(x))+x(1:length(x));
 rc(offset:(offset+length(x)-1))=rc(offset:(offset+length(x)-1))+x(1:length(x));
 rs = [xs,zeros(1,offset)] + [zeros(1,offset),xs];
 
 %filter with matched filter
 rss = xcorr(rs,xs);
 rcc = xcorr(rc,h);

 %set time
 N=floor(length(rcc)/2);
 time = (-N:+N)*(T);

if(length(time)<101)
   handles.time_scale=100;
end
 
 
half_length=round((N*(handles.time_scale/100))-1);

%set y limits
max_db=handles.max_db;
min_db=handles.min_db;
max_lin=10^(max_db/20);
min_lin=10^(min_db/20);

%plot matched filter output
axes(handles.two_pulse);
if(handles.plot_check==0)
    if(handles.normal_check==1)
        plot(time,abs(rss)/max(abs(rss+eps)), 'g', time, abs(rcc)/max(abs(rcc+eps)),'b');
        axis([time(N-half_length),time(N+half_length),min_lin,max_lin]);
    else
        plot(time,abs(rss), 'g', time, abs(rcc),'b');
        axis([time(N-half_length),time(N+half_length),min_lin,max_lin]);
    end
    xlabel('Amplitude')
else
    if(handles.normal_check==1)       
        plot(time,(20*log10(abs(rss)+eps))-max(20*log10(abs(rss)+eps)), 'g', time, (20*log10(abs(rcc)+eps))-max(20*log10(abs(rcc)+eps)),'b');
        axis([time(N-half_length),time(N+half_length),min_db,max_db]);
    else
        plot(time,20*log10(abs(rss)+eps), 'g', time, 20*log10(abs(rcc)+eps),'b');
        axis([time(N-half_length),time(N+half_length),min_db,max_db]);
    end
    xlabel('Amplitude (db)')
end
title('Received Signal with Matched Filter with 2 Scatterers')
ylabel('Delay (sec)')

 


guidata(hObject, handles);

%create contour plot of ambiguity function
function ambig(hObject, eventdata, handles)
global TT
global XX
global XS
global HH
T=TT;
x = XX;
xs = XS;
h=HH;


N=length(x);
M= 2*N-1;
AF = zeros(2*N+1,M);
for n=-N+2:N
  ll = max(1,n);
  ul = min(N,n+N);
  if (n <= 1)
      AF(n+N,1:n+N-1) = x(1:n+N-1).*conj(x(2-n:N));
  else
      AF(n+N,n:N) = x(n:N).*conj(x(1:N-n+1));
  end
end
for n = 1:2*N+1
 AF(n,:) = fftshift(M*ifft(AF(n,:)));
end
AF = abs(AF);
peak=max(max(AF));
f_abs = ((-M/2:M/2-1)/M)*(1/T);
t_abs = (-N:N)*T;
c=[0.994*peak,0.708*peak,0.316*peak,0.1*peak];
axes(handles.contour_plot);
contour(f_abs,t_abs,AF,c);
axis([min(f_abs),max(f_abs),min(t_abs),max(t_abs)])
title('Ambiguity Function')
xlabel('Doppler shift (Hz)')
ylabel('Delay (sec)')

%updates the two-scatterers matched filter output graph
% --- Executes on button press in update_button.
function update_button_Callback(hObject, eventdata, handles)
% hObject    handle to update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_plots(hObject, eventdata, handles)

%plots ambiguity function
%this is not done automatically because of the amount of time it takes
% --- Executes on button press in ambig_button.
function ambig_button_Callback(hObject, eventdata, handles)
% hObject    handle to ambig_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ambig(hObject, eventdata, handles)

%sets the distance between scatterers
function dist_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dist_edit as text
%        str2double(get(hObject,'String')) returns contents of dist_edit as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.dist = NewVal;
guidata(hObject, handles);
new_plots(hObject, eventdata, handles)

%sets output plot to log scale
% --- Executes on button press in plot_check.
function plot_check_Callback(hObject, eventdata, handles)
% hObject    handle to plot_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_check

if(handles.plot_check==0)
    handles.plot_check=1;
else
    handles.plot_check=0;
end
guidata(hObject, handles);
new_plots(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function dist_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%changes time scale when the slider moves
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to time_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
time_scale= get(hObject, 'value');
time_scale=round(time_scale);
if(time_scale==0)
    time_scale=1;
end
handles.time_scale=time_scale;
time_text_CreateFcn(hObject, eventdata, handles);
guidata(hObject, handles);
new_plots(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%changes time_text when slider changes
% --- Executes during object creation, after setting all properties.
function time_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hfs = findobj('tag', 'slider1');
time_scale = get(hfs,'value');
if (iscell(time_scale))
    time_scale=time_scale{1};
end
time_scale=round(time_scale);
if(time_scale==0)
    time_scale=1;
end
hfstext = findobj('tag','time_text');
set(hfstext, 'string', time_scale);
handles.time_scale=time_scale;
guidata(hObject, handles);

%sets maximum y to maxdb
function max_db_Callback(hObject, eventdata, handles)
% hObject    handle to max_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_db as text
%        str2double(get(hObject,'String')) returns contents of max_db as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.max_db = NewVal;
guidata(hObject, handles);
new_plots(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function max_db_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%sets minimium y to mindb
function min_db_Callback(hObject, eventdata, handles)
% hObject    handle to min_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_db as text
%        str2double(get(hObject,'String')) returns contents of min_db as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.min_db = NewVal;
guidata(hObject, handles);
new_plots(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function min_db_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_db (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%nomalizes output graph
% --- Executes on button press in normal_check.
function normal_check_Callback(hObject, eventdata, handles)
% hObject    handle to normal_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normal_check

if(handles.normal_check==0)
    handles.normal_check=1;
else
    handles.normal_check=0;
end
guidata(hObject, handles);
new_plots(hObject, eventdata, handles)
