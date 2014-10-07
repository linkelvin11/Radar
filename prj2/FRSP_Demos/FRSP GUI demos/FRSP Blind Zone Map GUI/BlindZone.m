function varargout = BlindZone(varargin)
% 
% BlindZone - GUI-based demonstration of blind zone maps for an M-of-N
% detection process
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

% BLINDZONE M-file for BlindZone.fig
%      BLINDZONE, by itself, creates a new BLINDZONE or raises the existing
%      singleton*.
%
%      H = BLINDZONE returns the handle to a new BLINDZONE or the handle to
%      the existing singleton*.
%
%      BLINDZONE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLINDZONE.M with the given input arguments.
%
%      BLINDZONE('Property','Value',...) creates a new BLINDZONE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BlindZone_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BlindZone_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help BlindZone

% Last Modified by GUIDE v2.5 21-Nov-2005 00:40:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BlindZone_OpeningFcn, ...
                   'gui_OutputFcn',  @BlindZone_OutputFcn, ...
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


% --- Executes just before BlindZone is made visible.
function BlindZone_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BlindZone (see VARARGIN)

% Choose default command line output for BlindZone
handles.output = hObject;

% Update handles structure
handles.CPI=50e-3;
handles.PRI=100*1e-6;
handles.tau=1e-6;
handles.CV=20;
handles.NearEnd=2;
handles.thresh=1;
handles.buttonVal = 1;

set(handles.radiobutton1, 'Value',1);

tau=handles.tau;
CPI=handles.CPI;
f=10e9;%x-band
c=3e8;

lambda=c/f;
range_res=c*tau/2;
vel_res=(1/CPI)*(lambda/2);

hfstext = findobj('tag','range_text');
set(hfstext, 'string', range_res/1e3);
hfstext = findobj('tag','vel_text');
set(hfstext, 'string', vel_res);

guidata(hObject, handles);


% UIWAIT makes BlindZone wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BlindZone_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function PRI_edit_Callback(hObject, eventdata, handles)
% hObject    handle to PRI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PRI_edit as text
%        str2double(get(hObject,'String')) returns contents of PRI_edit as a double
hfs = findobj('tag', 'PRI_edit');
NewStrVal = get(hfs, 'string');
%if isCell(NewStrVal)
%    temp=NewStrVal{1};
%else
%    temp = NewStrVal;
%end
NewVal = str2num(NewStrVal);
handles.PRI=NewVal*1e-6;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function PRI_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PRI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function tau_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tau_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau_edit as text
%        str2double(get(hObject,'String')) returns contents of tau_edit as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.tau = NewVal*1e-6;

tau=NewVal*1e-6;
CPI=handles.CPI;
f=10e9;%x-band
c=3e8;

lambda=c/f;
range_res=c*tau/2;
vel_res=(1/CPI)*(lambda/2);

hfstext = findobj('tag','range_text');
set(hfstext, 'string', range_res/1e3);
hfstext = findobj('tag','vel_text');
set(hfstext, 'string', vel_res);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tau_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function CPI_edit_Callback(hObject, eventdata, handles)
% hObject    handle to CPI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CPI_edit as text
%        str2double(get(hObject,'String')) returns contents of CPI_edit as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.CPI = NewVal*1e-3;
CPI=NewVal*1e-3;
tau= handles.tau;
f=10e9;%x-band
c=3e8;

lambda=c/f;
range_res=c*tau/2;
vel_res=(1/CPI)*(lambda/2);

hfstext = findobj('tag','range_text');
set(hfstext, 'string', range_res/1e3);
hfstext = findobj('tag','vel_text');
set(hfstext, 'string', vel_res);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function CPI_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CPI_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot_spots(hObject, eventdata, handles)

function plot_spots(hObject, eventdata, handles)
%initialize user-specified values

CPI_edit_Callback(hObject, eventdata, handles)
tau_edit_Callback(hObject, eventdata, handles)
PRI_edit_Callback(hObject, eventdata, handles)

edit_thresh_Callback(hObject, eventdata, handles)
thresh = handles.thresh;

CPI=handles.CPI;
tau=handles.tau;
PRI=handles.PRI;
CV=handles.CV;

if ((thresh>length(PRI))&&(handles.buttonVal==-1))
    errordlg('Your threshold is greater than the number of PRIs.')
    error('Your threshold is greater than the number of PRIs.')
end

%initialize hardcoded values
f=10e9;%x-band
c=3e8;

lambda=c/f;
range_res=c*tau/2;
vel_res=(1/CPI)*(lambda/2);

hfstext = findobj('tag','range_text');
set(hfstext, 'string', range_res/1e3);
hfstext = findobj('tag','vel_text');
set(hfstext, 'string', vel_res);


%POSSIBLE CHANGING VALUES BASED ON CLUTTER
vel_width = ceil(CV/vel_res);
%clutter blinds things within +/- 20 m/s 
%range_width = round(((tau*c)+20)/range_res);
range_width = handles.NearEnd;
%clutter interferes with thing within 20 m

PRF=1./PRI;

y_max=5*c*max(PRI)/2;
x_max=4*max(PRF)*lambda/2;

y=0:range_res:y_max;
x=0:vel_res:x_max;


Y=length(y);
X=length(x);



total_area=zeros(Y,X);
v1_vect=1;
start_vect=1;
for k = 1:length(PRI)

area=ones(Y,X);


area(:,1:vel_width)=0;

v1=PRF(k)*lambda/2;
v=v1;
%v1=abs((f/(f+freq1)-1)*c)
v_index=round(v/vel_res);
v_vect=[];
v_vect=[v_vect,v_index];
n=2;
while (v_index-vel_width<X)
    if((v_index+vel_width)<X)
        area(:,(v_index-vel_width):(v_index+vel_width))=0;
    else
        area(:,(v_index-vel_width):X)=0;
    end
    v=v1*n;
    v_index=round(v/vel_res);
    v_vect=[v_vect,v_index];
    n=n+1;
end

start = 1;
while (start < Y)
     if((start+range_width)<Y)
            area(start:(start+range_width),:)=0;
     else
         area(start:Y,:)=0;
     end
     start = start+round((PRI(k)*c/2)/range_res);
     start_vect=[start_vect,start];
end



total_area=total_area+area;
end

if (handles.buttonVal==1)
    total_area=total_area./length(PRI);
else
    total_area=floor(total_area./handles.thresh);
    total_area=min(total_area,1);
end


%y=y(length(y):-1:1);
y=y/(1e3);

velocity_blind_zones = v_vect*vel_res;
range_blind_zones=start_vect*range_res-range_res;




axes(handles.plot);
axis xy
imagesc(x,y,total_area);
axis xy
%rotate(Z,[0,0,1],180);
xlabel('Velocity (m/s)');
ylabel('Range (km)');
colormap(gray)

%figure
%imagesc(x,y,total_area);



function edit_CV_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CV as text
%        str2double(get(hObject,'String')) returns contents of edit_CV as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.CV = NewVal;



guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function edit_CV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function edit_NearEnd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NearEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NearEnd as text
%        str2double(get(hObject,'String')) returns contents of edit_NearEnd as a double
NewStrVal = get(hObject, 'string');
NewVal = str2double(NewStrVal);
handles.NearEnd = NewVal;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_NearEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NearEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




function edit_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresh as text
%        str2double(get(hObject,'String')) returns contents of edit_thresh as a double
hfs = findobj('tag', 'edit_thresh');
NewStrVal = get(hfs, 'string');
NewVal = str2double(NewStrVal);
handles.thresh = NewVal;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

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

%mutual_exclude turns off other radio buttons
function mutual_exclude(off)
set(off,'Value',0)



