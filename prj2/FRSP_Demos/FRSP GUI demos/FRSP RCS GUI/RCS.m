function varargout = RCS(varargin)

% 
% RCS - GUI-based demonstration of the RCS statistics observed with a
% multiple-point scatterer "target".
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

% RCS M-file for RCS.fig
%      RCS, by itself, creates a new RCS or raises the existing
%      singleton*.
%
%      H = RCS returns the handle to a new RCS or the handle to
%      the existing singleton*.
%
%      RCS('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in RCS.M with the given input arguments.
%
%      RCS('Property','Value',...) creates a new RCS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RCS_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RCS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help RCS

% Last Modified by GUIDE v2.5 05-Sep-2005 17:56:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RCS_OpeningFcn, ...
                   'gui_OutputFcn',  @RCS_OutputFcn, ...
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


% --- Executes just before RCS is made visible.
function RCS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RCS (see VARARGIN)

% Choose default command line output for RCS
handles.output = hObject;

%This sets many of the default values that aren't set in the .fig file.

handles.checked=0;
handles.buttonVal = 1; %no dominant scatterers.
set(handles.radiobutton1, 'Value',1);%no dominant scatterers. 
handles.rice_val=1; %value in editable text box
%all checkboxes off
handles.check1 = 0;
handles.check2 = 0;
handles.check3 = 0;
% Update handles structure
guidata(hObject, handles);

%create texts
scat_text_CreateFcn(hObject, eventdata, handles);
distance_text_CreateFcn(hObject, eventdata, handles);
RCS_scatter_text_CreateFcn(hObject, eventdata, handles);
x_text_CreateFcn(hObject, eventdata, handles);
y_text_CreateFcn(hObject, eventdata, handles);

%end opening functions

% UIWAIT makes RCS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RCS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%**
%**All the mathmatical computations take place in the next two functions**
%**

% --- Executes on button press in new_scatter. (the button labeled "New
% Scatterers")

%this function creates and plots the scatterers
function new_scatter_Callback(hObject, eventdata, handles)
% hObject    handle to new_scatter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%find the number of scatters
hfs = findobj('tag', 'scatter_slide');
N = get(hfs,'value');
N=floor(N);
if(N==0)
    N=1;
end

%find the length in the x-direction
hfs2 = findobj('tag', 'x_length');
xl = get(hfs2,'value');
xl=xl-mod(xl,0.1);

%if set to zero, set x1 to machine epsilon
if(xl==0)
    xl=eps;
end

%find the length in the y-direction
hfs3 = findobj('tag', 'y_length');
yl = get(hfs3,'value');
yl=yl-mod(yl,0.1);


%if set to zero, set y1 to machine epsilon
if(yl==0)
    yl=eps;
end

%get random values and center around 0
y=(2*yl)*rand(N,1)-yl; 
x=(2*xl)*rand(N,1)-xl;
guidata(hObject, handles);
axes(handles.scatter_plot);


if(handles.buttonVal==1)%no dominant
    plot(x,y,'*');
else %add dominant scatterer
   y(N+1)=(2*yl)*rand(1,1)-yl; x(N+1)=(2*xl)*rand(1,1)-xl;
   plot(x(1:N),y(1:N),'b*'); hold on;
   plot(x(N+1),y(N+1),'ro'); hold off;
end
    
title('Map of Scatterers');
xlabel('X Direction'); ylabel('Y Direction');
axis('equal'); axis([-xl,xl,-yl,yl]);
handles.x = x;
handles.y = y;
handles.N = N;
guidata(hObject, handles);%store data for other plots

%clear other plots
axes(handles.RCS_plot);
cla
axes(handles.histogram);
cla
axes(handles.autocorr_plot);
cla


% --- Executes on button press in plot_data. (button labeled "Compute RCS
% And Plot Results")

%This creates the RCS plot, histogram, and autocorrelation plots.
function plot_data_Callback(hObject, eventdata, handles)
% hObject    handle to plot_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%load the scatterer data
x=handles.x;
y=handles.y;
N=handles.N;

%check the Rician ratio number
rice_num_Callback(hObject, eventdata, handles);


%check RCS value of individual scatters (z)
hfs2 = findobj('tag', 'RCS_scatter_slider');
z = get(hfs2,'value');
z=z-mod(z,0.1);


%if set to zero, set z to machine epsilon
if(z==0)
    z=eps;
end


%check angle increments (inc)
hfs2 = findobj('tag', 'angle_menu');
val = get(hfs2,'value');

if val==1
    inc = 2;
elseif val==2
    inc = 1;
elseif val==3
    inc = 0.5;
elseif val==4
    inc = 0.2;
elseif val==5
    inc = 0.1;
elseif val==6
    inc = 0.05;
elseif val==7
    inc = 0.02;
elseif val==8
    inc = 0.01;
elseif val==9
    inc = 0.005;
elseif val==10
    inc = 0.002;
elseif val==11
    inc = 0.001;
end

M=360/inc;%M = number of angles

%check distance from center (R)
hfs3 = findobj('tag', 'distance_slide');
R = get(hfs3,'value');
R=round(R)*10^3;

%check frequency (f)
hfs4 = findobj('tag', 'freq_menu');
val = get(hfs4,'Value');
if val == 1 %HF
    f=3e6;
elseif val == 2 %VHF
    f=30e6;
elseif val == 3 %UHF
    f=300e6;
elseif val == 4 %L
    f=1e9;
elseif val == 5 %S
    f=3e9;
elseif val == 6 %C
    f=6e9;
elseif val == 7 %X
    f=10e9;
elseif val == 8 %Ku
    f=16e9;
elseif val == 9 %Ka
    f=35e9;
elseif val == 10 %W
    f=95e9;
end

%**Calculate echo values**

% Loop over aspect angle (k) to do complex voltage estimates
q=zeros(M,1); echo=zeros(M,1);
for k=1:M
   q(k)=(2*pi/M)*(k-1); % current aspect angle
   % Loop over individual point scatterers (p)
  for p=1:N
     phasor=z*exp(j*4*pi*f*norm([x(p)-R*cos(q(k)),y(p)-R*sin(q(k))])/3e8);
     echo(k)=echo(k)+phasor;
  end
end


if(handles.buttonVal~=1)%if dominant scatter exists
    if(handles.buttonVal==2)%chi-squared
        w=sqrt((1+sqrt(2))*N*z);
    elseif (handles.buttonVal==3)%rician
        hfs3 = findobj('tag', 'rice_num');
        NewStrVal = get(hfs3, 'string');
        amp = str2double(NewStrVal);
        w=sqrt(amp*N*z);
    end
    %w = RCS of dominant scatter
    %loop over all angles add the calculation for the dominant scatter in
    for k=1:M
        q(k)=(2*pi/M)*(k-1);
        % add in the dominant scatterer
        phasor=w*exp(j*4*pi*f*norm([x(N+1)-R*cos(q(k)),y(N+1)-R*sin(q(k))])/3e8);
        echo(k)=echo(k)+phasor;
    end
end

%**The magnitude-squared of the total complex echo
%is the RCS to within a constant**
RCS=abs(echo).^2;
RCSdB=10*log10(RCS);
axes(handles.RCS_plot);
plot((180/pi)*q,RCSdB); xlabel('Aspect Angle (degrees)');
ylabel('RCS (dB)');
title('RCS vs. Aspect Angle')';
maxDB=max(RCSdB);
maxDB=maxDB-mod(maxDB,10)+10;
minDB=min(RCSdB);
minDB=minDB-mod(minDB,10)-10;
axis([0,360,minDB,maxDB]);

%**Compute histogram of the data**
[n,sig]=hist(RCS,100);
n=n/(sum(n)*(sig(2)-sig(1)));  % normalizes to unit area
sum(n)*(sig(2)-sig(1));
axes(handles.histogram);
bar(sig,n);
hold on

%add fit line
if(handles.check1==1)%exponential
    msig=mean(RCS);
    % plot(sig,exp(-sig/msig))
    plot(sig,exp(-sig/msig)/msig,'r','LineWidth',3);
end
if (handles.check2==1)%chi-squared
    msig=mean(RCS);
    chisquare=(4*sig/(msig^2)).*exp(-2*sig/msig);
    % plot(sig,chisquare/max(chisquare)); old normalization, not used anymore
    plot(sig,chisquare,'g','LineWidth',3);
end
if (handles.check3==1)%rician
    msig=mean(RCS);
    if(handles.buttonVal==1)
        asq=0;
    elseif(handles.buttonVal==2)
        asq=0;
    elseif(handles.buttonVal==3)
        hfs3 = findobj('tag', 'rice_num');
        NewStrVal = get(hfs3, 'string');
        asq = str2double(NewStrVal);
    end
    riceA=(1/msig)*(1+asq)*exp(-asq-(sig/msig)*(1+asq));
    I0=besseli(1,2*sqrt(asq)*sqrt((1+asq)*(sig/msig)));
    rice=riceA.*I0;
    plot(sig,rice,'y','LineWidth',3);
end

xlabel('Radar Cross Section (m^2)')
ylabel('Relative Probability')
title('PDF of RCS')
hold off

%**Compute autocorrelation**

%find section to plot
hfs3 = findobj('tag', 'thetamax');
NewStrVal = get(hfs3, 'string');
thetamax = str2double(NewStrVal);


thetacorr = 3e8/2/5/f*180/pi;
Dtheta = 0.1;
M2 = round(thetamax/Dtheta);

% Loop over aspect angle to do complex voltage estimates
q2=zeros(M2+1,1); echo2=zeros(M2+1,1);
for k=1:M2+1
   q2(k)=(pi/180)*((k-M2/2-1)*Dtheta); % current aspect angle
   % Loop over individual point scatterers
  for p=1:N
     phasor2=z*exp(j*4*pi*f*norm([x(p)-R*cos(q2(k)),y(p)-R*sin(q2(k))])/3e8);
     echo2(k)=echo2(k)+phasor2;
  end
end
s = xcorr(echo2); s = abs(s)/max(abs(s));
lag = Dtheta*(-M2:+M2);

axes(handles.autocorr_plot);

plot(lag,s); xlabel('Aspect Angle (Degrees)');
ylabel('Normalized Magnitude of Autocorrelation');
title('Angular Autocorrelation of RCS')
grid;

%add measurement lines
for(k=M2+1:2*M2+1)
    if(s(k)<0.25)
        hold on
        plot([-lag(k),-lag(k)],[0,1],'-g','LineWidth',2);
        plot([lag(k),lag(k)],[0,1],'-g','LineWidth',2);
        hold off
        break
    end
end

for(k=M2+1:2*M2)
    if(s(k)<s(k+1))
        hold on
        plot([-lag(k),-lag(k)],[0,1],'-k','LineWidth',2);
        plot([lag(k),lag(k)],[0,1],'-k','LineWidth',2);
        hold off
        break
    end
end

if(thetacorr<thetamax)
    hold on
    plot([-thetacorr,-thetacorr],[0,1],'-r','LineWidth',2);
    plot([thetacorr,thetacorr],[0,1],'-r','LineWidth',2);
    hold off
end

axis([-thetamax,thetamax,0,1]);
%**END PLOTTING FUNCTION**






%**The following functions are for visual display rather than calculation**

%**Number of scatterers**
% --- Executes on slider movement.
function scatter_slide_Callback(hObject, eventdata, handles)
% hObject    handle to scatter_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.textaboutnumscatters = get(hObject, 'value');
guidata(hObject, handles);
scat_text_CreateFcn(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function scatter_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scatter_slide (see GCBO)
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

% --- Executes during object creation, after setting all properties.
function scat_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scat_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hfs = findobj('tag', 'scatter_slide');
num_scatter = get(hfs,'value');
num_scatter=floor(num_scatter);
%do not allow to go to zero
if(num_scatter==0)
    num_scatter=1;
end
hfstext = findobj('tag','scat_text');
set(hfstext, 'string', num_scatter);



%**RCS of individual scatterers**

% --- Executes on slider movement.
function RCS_scatter_slider_Callback(hObject, eventdata, handles)
% hObject    handle to RCS_scatter_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.textaboutscatterrcs = get(hObject, 'value');
guidata(hObject, handles);
RCS_scatter_text_CreateFcn(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function RCS_scatter_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RCS_scatter_slider (see GCBO)
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

% --- Executes during object creation, after setting all properties.
function RCS_scatter_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RCS_scatter_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hfs = findobj('tag', 'RCS_scatter_slider');
scatterRCS = get(hfs,'value');
scatterRCS=scatterRCS-mod(scatterRCS,0.1);
hfstext = findobj('tag','RCS_scatter_text');
set(hfstext, 'string', scatterRCS);


%**x and y lengths**

% --- Executes during object creation, after setting all properties.
function x_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hfs = findobj('tag', 'x_length');
xl = get(hfs,'value');
xl=xl-mod(xl,0.1);
hfstext = findobj('tag','x_text');
set(hfstext, 'string', xl*2);

% --- Executes during object creation, after setting all properties.
function y_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hfs = findobj('tag', 'y_length');
yl = get(hfs,'value');
yl=yl-mod(yl,0.1);
hfstext = findobj('tag','y_text');
set(hfstext, 'string', yl*2);

% --- Executes on slider movement.
function x_length_Callback(hObject, eventdata, handles)
% hObject    handle to x_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.x_length = get(hObject, 'value');
guidata(hObject, handles);
x_text_CreateFcn(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function x_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_length (see GCBO)
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


% --- Executes on slider movement.
function y_length_Callback(hObject, eventdata, handles)
% hObject    handle to y_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.y_length = get(hObject, 'value');
guidata(hObject, handles);
y_text_CreateFcn(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function y_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_length (see GCBO)
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



%**Distance from the center**

% --- Executes on slider movement.
function distance_slide_Callback(hObject, eventdata, handles)
% hObject    handle to distance_slide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.distance = get(hObject, 'value');
guidata(hObject, handles);
distance_text_CreateFcn(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function distance_slide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance_slide (see GCBO)
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

% --- Executes during object creation, after setting all properties.
function distance_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hfs = findobj('tag', 'distance_slide');
distance_R = get(hfs,'value');
distance_R=round(distance_R);
hfstext = findobj('tag','distance_text');
set(hfstext, 'string', distance_R);



%**Frequency menu**

% --- Executes on selection change in freq_menu.
function freq_menu_Callback(hObject, eventdata, handles)
% hObject    handle to freq_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns freq_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from freq_menu
val = get(hObject,'Value');
str = get(hObject, 'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function freq_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




%**Radio buttons for dominant scatterer info**

%this function turns off the other two radio buttons when one is selected
function mutual_exclude(off)
set(off,'Value',0)

% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
off = [handles.radiobutton2,handles.radiobutton3];
mutual_exclude(off)
handles.buttonVal = 1;
guidata(hObject, handles);

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
off = [handles.radiobutton1,handles.radiobutton3];
mutual_exclude(off)
handles.buttonVal = 2;
guidata(hObject, handles);

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
off = [handles.radiobutton1,handles.radiobutton2];
mutual_exclude(off)
handles.buttonVal = 3;
guidata(hObject, handles);


%**Rician number editable text box info**

function rice_num_Callback(hObject, eventdata, handles)
% hObject    handle to rice_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rice_num as text
%        str2double(get(hObject,'String')) returns contents of rice_num as a double

NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.rice_val = NewVal;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function rice_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rice_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



%**checkboxes for fit line (tagged as radiobutton4-radiobutton6 for
%parallelism)**

% --- Executes on button press in radiobutton1.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1

if (get(hObject,'Value') == get(hObject,'Max'))
    handles.check1 = 1;
else
    handles.check1 = 0;
end
guidata(hObject, handles);

% --- Executes on button press in radiobutton2.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
if (get(hObject,'Value') == get(hObject,'Max'))
    handles.check2 = 1;
else
    handles.check2 = 0;
end
guidata(hObject, handles);

% --- Executes on button press in radiobutton3.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (get(hObject,'Value') == get(hObject,'Max'))
    handles.check3 = 1;
else
    handles.check3 = 0;
end
guidata(hObject, handles);



%**angle increment menu (here reffered to as popupmenu2)**

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




%**set the maximum angle for plotting the autocorrelation**

function thetamax_Callback(hObject, eventdata, handles)
% hObject    handle to thetamax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thetamax as text
%        str2double(get(hObject,'String')) returns contents of thetamax as a double


% --- Executes during object creation, after setting all properties.
function thetamax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thetamax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


