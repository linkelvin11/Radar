function varargout = LFM(varargin)
% 
% LFM - GUI-based demonstration of linear FM ("chirp") waveforms and
% matched filters, and comparison to simple pulse waveform and matched
% filter.
% 
% Permission to use, copy, modify and distribute, including the right
% to grant others rights to distribute at any tier, this software and
% its documentation for any purpose and without fee or royalty is hereby
% granted, provided that the asscoiated copyright and attribution, and
% disclaimer are retained in ALL copies and derivative works of the
% software and documentation, including modifications that you make
% for internal use or for distribution.
% 
% All software is provided "as is".  It is intended for tutorial
% and demonstration use only and is provided as a convenience and
% courtesy to the user.  No user support is available.
% 
% The developer, the Georgia Institute of Technology, and the
% Distance Learning and Professional Education division of the
% Georgia Institute of Technology make absolutely no warranty of
% merchantability or fitness for any use for any particular purpose,
% or that the use of the software will not infringe on any third
% party patents, copyrights, trademarks, or other rights.
% 
% Except as otherwise noted, this software was developed by
% Dr. Mark A. Richards and/or Gregory Heim of the Georgia Institute
% of Technology, and was provided by the Georgia Institute of Technology
% as part of the professional education course “Fundamentals of Radar
% Signal Processing”, 2006 and subsequent offerings, or by Dr. Mark A.
% Richards to accompany the book "Fundamentals of Radar Signal
% Processing", Mark A. Richards (author), published by McGraw-Hill,
% New York, 2005.
% 
% Copyright 2006-2009, Mark A. Richards and/or Gregory Heim.
% All Rights Reserved.
% 

%
% LFM M-file for LFM.fig
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
                   'gui_OpeningFcn', @LFM_OpeningFcn, ...
                   'gui_OutputFcn',  @LFM_OutputFcn, ...
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


% --- Executes just before LFM is made visible.
function LFM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LFM (see VARARGIN)

% Choose default command line output for LFM
handles.output = hObject;

%default values
handles.BW=10; %initial bandwidth (MHz)
handles.pulse=10; %initial pulse length (microseconds)
handles.k=2; %initial oversampling ratio
handles.BT=(handles.BW*1e6)*(handles.pulse*1e-6); %time-bandwidth product
handles.plot_check=1; %sets the output plot to log value
handles.normal_check=0; %initializes the output plot to unnormalized
handles.phase_check=0; %initializes to wrapped phase plot
handles.time_scale=100; %initializes to full time extent (no zoom)
handles.max_db=46; %intializes maximum db on output plot
handles.min_db=-10; %initializes minimum db on the plot

%set BT text
hfstext = findobj('tag','BTtext');
set(hfstext, 'string', handles.BT);
time_text_CreateFcn(hObject, eventdata, handles);


handles.buttonVal = 1; %upchirp/downchirp radio buttons
set(handles.radiobutton1, 'Value',1);


% Update handles structure
guidata(hObject, handles);
plotchirp(hObject, eventdata, handles)


BTtext_CreateFcn(hObject, eventdata, handles);
% UIWAIT makes LFM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LFM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%plotchirp plots all four graphs, it is called when the window is openened
%and when changes to the signal are made
function plotchirp(hObject, eventdata, handles)

 BW=handles.BW*1e6; %convert to megahertz
 pulse=handles.pulse*1e-6; %convert to microseconds
 k=handles.k;
 %create chirp function
if(BW~=0)
    T=1/(handles.k*BW);
    t=(-pulse/2):T:(pulse/2);
    if(handles.buttonVal==1);
        x=exp((j*pi*BW/pulse)*t.^2);
    else
        x=exp(-(j*pi*BW/pulse)*t.^2);
    end
    change_back=0;
elseif(BW==0)
    handles.BW=1/pulse;
    BW=handles.BW;
    T=1/(handles.k*BW);
    t=(-pulse/2):T:(pulse/2);
    x=ones(1,length(t));
    change_back=1;%flag to change bandwidth back to zero at the end of the function
end

handles.t=t;
handles.x=x;
%plot chirp
axes(handles.chirp_graph);
plot(handles.t,real(handles.x));
title('Output Chirp Signal')
xlabel('Time (sec)')
ylabel('Amplitude')

%convert to frequency domain
Nfft = 2*length(x);
if(Nfft<256)
    Nfft=256;
end
Nf2 = Nfft/2;
X = fftshift(fft(x,Nfft));
freq = ((0:Nfft-1)-Nf2)/Nfft;
nspec = round(Nfft/k);


%find window
hfs = findobj('tag', 'window_menu');
win_val = get(hfs,'value');

if(win_val~=1)
    if(win_val==2)
       window = hamming(nspec)';
    elseif(win_val==3)
        window = hann(nspec)';
    elseif(win_val==4)
       window = blackman(nspec)';
    end
    window=[zeros(1,Nf2-1-floor(nspec/2)),window,zeros(1,Nf2+1-ceil(nspec/2))];
end

%plot frequency domain
axes(handles.FFT_graph);
if(win_val==1)
    plot(linspace(-.5*(handles.BW*1e6)*handles.k,.5*(handles.BW*1e6)*handles.k,Nfft),abs(X));
    axis([-.5*(handles.BW*1e6)*handles.k,.5*(handles.BW*1e6)*handles.k,min(abs(X)),max(abs(X))]);
    title('Magnitude of Chirp Frequency Response')
    xlabel('Frequency (Hz)/(Radians)')
    ylabel('Magnitude')
else
    plot(linspace(-.5*(handles.BW*1e6)*handles.k,.5*(handles.BW*1e6)*handles.k,Nfft),abs(X)/max(abs(X)),'b',linspace(-.5*(handles.BW*1e6)*handles.k,.5*(handles.BW*1e6)*handles.k,Nfft),window,'g');
    axis([-.5*(handles.BW*1e6)*handles.k,.5*(handles.BW*1e6)*handles.k,0,1.2]);
    title('Normalized Magnitude of Chirp Frequency Response and Window')
    xlabel('Frequency (Hz)/(Radians)')
    ylabel('Normalized Magnitude')
end    
            
axes(handles.FFT_graph2);
NT=round(Nfft/(2*handles.k))-1;
start=round(Nfft/2)-NT;
ending=round(Nfft/2)+NT;
freq=linspace(-.5*(handles.BW*1e6)*handles.k,.5*(handles.BW*1e6)*handles.k,Nfft);
if(handles.phase_check==0)
    plot(freq,angle(X));
else
    plot(freq(start:ending),unwrap(angle(X(start:ending))));
end
title('Phase of Chirp Frequency Response')
xlabel('Frequency (Hz)/(Radians)')
ylabel('Phase')


%alter filter by windowing frequency domain
if(win_val~=1)
    H=window.*conj(X);
    h=fftshift(ifft(fftshift(H),Nfft));
    h=conj(h(1:Nf2));

    if(change_back~=1)
        h(1:round(Nf2/2))=h(length(h):-1:round(Nf2/2));
    end
    xcc=xcorr(x,h);
   
% do matched filter output (autocorrelation function) of each
% waveform and compare on same scale
    
else
    h=x;
    xcc=xcorr(x,h);
end

N=length(x);

% create simple pulse of same duration and sampling rate
xs=ones(1,length(t));

%set time
time = (-N+1:+N-1)*(T);

%find matched filter outputs
xcc=xcc((length(xcc)-length(time)+1:length(xcc)));
xss = xcorr(xs);


if(N<101)
    handles.time_scale=100;
    %do not allow small plots to be scaled down
end

%scale time for x limits
half_length=round((length(time)*(handles.time_scale/100))/2)-1;

%set y limits
max_db=handles.max_db;
min_db=handles.min_db;
max_lin=10^(max_db/20);
min_lin=10^(min_db/20);

%plot matched filter output    
axes(handles.matchfilt_output);
if(handles.plot_check==0)
    if(handles.normal_check==1)
        plot(time,abs(xss)/max(abs(xss+eps)),'g',time, abs(xcc)/max(abs(xcc+eps)),'b'); 
        axis([time(N-half_length),time(N+half_length),min_lin,max_lin]);
    else
        plot(time,abs(xss),'g',time, abs(xcc),'b'); 
        axis([time(N-half_length),time(N+half_length),min_lin,max_lin]);
    end
    ylabel('Amplitude')
else
    if(handles.normal_check==1)
        plot(time,20*log10(abs(xss)+eps)-(max(20*log10(abs(xss)+eps))),'g',time, 20*log10(abs(xcc+eps))-(max(20*log10(abs(xcc)+eps))),'b');
        axis([time(N-half_length),time(N+half_length),min_db,max_db]);
    else
        plot(time,20*log10(abs(xss)+eps),'g',time, 20*log10(abs(xcc+eps)),'b');
        axis([time(N-half_length),time(N+half_length),min_db,max_db]);
    end
    ylabel('Amplitude (dB)'); 
end

title('Received Signal with Matched Filter')
xlabel('Time (sec)')

%set handle variables
handles.xcc = xcc;
handles.xss = xss;
handles.N=N;
handles.half_length=half_length;
handles.time=time;
handles.X=X;
handles.Nfft=Nfft;

%set global variables so that they may accessed by the pop-up gui
global TT
TT = T;
global XX
XX = x;
global XS
XS = xs;
global HH
HH=h;
global NN
NN=N;

%change bandwidth back to zero
if(change_back==1)
    handles.BW=0;
end

guidata(hObject, handles);


%opens up pop-up window IF there are enough samples to plot the
%two-scatterer output
% --- Executes on button press in add_plots.
function add_plots_Callback(hObject, eventdata, handles)
% hObject    handle to add_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if((handles.BW~=0)||(handles.k>=15))
    LFM_popup (handles)
else
    errordlg('There are not enough samples to plot these graphs. If you wish to plot the current waveform please increase the oversampling rate.')
end


%edits bandwidth
function edit_BW_Callback(hObject, eventdata, handles)
% hObject    handle to edit_BW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_BW as text
%        str2double(get(hObject,'String')) returns contents of edit_BW as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.BW = NewVal;
handles.BT=(handles.BW*1e6)*(handles.pulse*1e-6);
hfstext = findobj('tag','BTtext');
set(hfstext, 'string', handles.BT);
BTtext_CreateFcn(hObject, eventdata, handles)
guidata(hObject, handles);
plotchirp(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_BW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_BW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%edit pulse length
function pulse_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pulse_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pulse_edit as text
%        str2double(get(hObject,'String')) returns contents of pulse_edit as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.pulse = NewVal;
handles.BT=(handles.BW*1e6)*(handles.pulse*1e-6);
hfstext = findobj('tag','BTtext');
set(hfstext, 'string', handles.BT);
BTtext_CreateFcn(hObject, eventdata, handles)

guidata(hObject, handles);
plotchirp(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pulse_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulse_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%edit oversampling rate
function k_edit_Callback(hObject, eventdata, handles)
% hObject    handle to k_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_edit as text
%        str2double(get(hObject,'String')) returns contents of k_edit as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.k = NewVal;
guidata(hObject, handles);
plotchirp(hObject, eventdata, handles)
BTtext_CreateFcn(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function k_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%edits BT text
% --- Executes during object creation, after setting all properties.
function BTtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BTtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%changes radiobutton that controls upchirp vs downchirp
% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
off = [handles.radiobutton2];
mutual_exclude(off)
handles.buttonVal = 1;
guidata(hObject, handles);
plotchirp(hObject, eventdata, handles)

%changes radiobutton that controls upchirp vs downchirp
% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
off = [handles.radiobutton1];
mutual_exclude(off)
handles.buttonVal = -1;
guidata(hObject, handles);
plotchirp(hObject, eventdata, handles)

%mutual_exclude turns off other radio buttons
function mutual_exclude(off)
set(off,'Value',0)

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
plotchirp(hObject, eventdata, handles)

%replots when window value changes
% --- Executes on selection change in window_menu.
function window_menu_Callback(hObject, eventdata, handles)
% hObject    handle to window_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns window_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from window_menu
plotchirp(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function window_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to window_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%changes time scale when the slider moves
% --- Executes on slider movement.
function time_slider1_Callback(hObject, eventdata, handles)
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
plotchirp(hObject, eventdata, handles)

    
% --- Executes during object creation, after setting all properties.
function time_slider1_CreateFcn(hObject, eventdata, handles)
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
hfs = findobj('tag', 'time_slider1');
time_scale = get(hfs,'value');
time_scale=round(time_scale);
if(time_scale==0)
    time_scale=1;
end
hfstext = findobj('tag','time_text');
set(hfstext, 'string', time_scale);
handles.time_scale=time_scale;
guidata(hObject, handles);


%sets phase to wrapped or unwrapped, and replots phase graph
% --- Executes on button press in phase_check.
function phase_check_Callback(hObject, eventdata, handles)
% hObject    handle to phase_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of phase_check
if(handles.phase_check==0)
    handles.phase_check=1;
else
    handles.phase_check=0;
end
guidata(hObject, handles);
 if (handles.BW~=0)
    X=handles.X;
    Nfft=handles.Nfft;
    NT=round(Nfft/(2*handles.k))-1;
    start=round(Nfft/2)-NT;
    ending=round(Nfft/2)+NT;
    axes(handles.FFT_graph2);
    freq=linspace(-.5*(handles.BW*1e6)*handles.k,.5*(handles.BW*1e6)*handles.k,Nfft);
    if(handles.phase_check==0)
        plot(freq,angle(X));
    else
        plot(freq(start:ending),unwrap(angle(X(start:ending))));
    end
    title('Phase of Chirp Frequency Response')
    xlabel('Frequency (Hz)/(Radians)')
    ylabel('Phase')
 else
     plotchirp(hObject, eventdata, handles)
 end


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
plotchirp(hObject, eventdata, handles)


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
plotchirp(hObject, eventdata, handles)

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
plotchirp(hObject, eventdata, handles)

