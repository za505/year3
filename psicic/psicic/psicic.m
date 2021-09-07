function varargout = psicic(varargin)
% psicic M-file for psicic.fig
%      psicic, by itself, creates a new psicic or raises the existing
%      singleton*.
%
%      H = psicic returns the handle to a new psicic or the handle to
%      the existing singleton*.
%
%      psicic('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in psicic.M with the given input arguments.
%
%      psicic('Property','Value',...) creates a new psicic or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before psicic_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to psicic_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help psicic

% Last Modified by GUIDE v2.5 10-Nov-2008 15:28:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @psicic_OpeningFcn, ...
    'gui_OutputFcn',  @psicic_OutputFcn, ...
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


% --- Executes just before psicic is made visible.
function psicic_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to psicic (see VARARGIN)


[handles.filename, handles.imdir] = uigetfile('*.tif','Select an image file');
cd(handles.imdir);

    handles.data.filename = fullfile(handles.imdir,handles.filename);
    filename = handles.data.filename;
    handles.img = imread(handles.data.filename);
    curimgadj = handles.img;
    handles.img(:,:,1) = imadjust(handles.img(:,:,1)); %EXPERIMENTAL
    fileinfo = imfinfo(filename);
    % Find fluorescent bit depth
    if size(handles.img,3) > 1
        handles.normfactor = (2^fileinfo.BitsPerSample(2)-1);
        if size(handles.img, 3) == 2
           handles.img(:,:,3) = 0; 
        end
    else
        set(handles.slider1,'Enable','Off');
    end

if size(handles.img,3)>1
    % Set current levels and display adjusted image
    handles.curlevel = (double(.5*max(max(handles.img(:,:,2))))/handles.normfactor);
    handles.autolevel = handles.curlevel;
    imshow(imadjust(handles.img, [0 0 0; 1 max(handles.curlevel,0.01) 1]))

    % Reset slider
    set(handles.slider1,'Value',1-handles.curlevel)
else
    imshow(handles.img)
    
end
set(handles.curfilename,'String',handles.filename);
set(handles.imhi, 'String', ['High: ' num2str(max(curimgadj(:)))]);
set(handles.imlo, 'String', ['Low: ' num2str(min(curimgadj(:)))]);
set(handles.immid, 'String', ['Mid: ' num2str( (min(curimgadj(:)) + max(curimgadj(:)))/2 )]);
%Uncomment next line to have automatic contourlevel picking
%set(handles.contourlevel, 'String', num2str( uint16(.9*((min(curimgadj(:)) + max(curimgadj(:)))/2)) ));

% Choose default command line output for psicic
handles.output = hObject;

handles.button_state = get(handles.togglebutton1,'Value');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes psicic wait for user response (see UIRESUME)
% uiwait(handles.psicic);


% --- Outputs from this function are returned to the command line.
function varargout = psicic_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function psicic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psicic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function mainfig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mainfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate mainfig


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.curlevel = 1-get(hObject,'Value');
redisplay(handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1

handles.button_state = get(hObject,'Value');
redisplay(handles);

guidata(hObject, handles);


% --- Executes on button press in autoadjust.
function autoadjust_Callback(hObject, eventdata, handles)
% hObject    handle to autoadjust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curlevel = handles.autolevel;
set(handles.slider1,'Value',1-handles.curlevel)
redisplay(handles);

guidata(hObject, handles);


% Mine: redisplay image according to togglebutton1 state
function redisplay(handles)
if size(handles.img,3)>1
    if handles.button_state == get(handles.togglebutton1,'Max')
        imshow(imadjust(handles.img(:,:,2), [0; handles.curlevel]))
    else
        imshow(imadjust(handles.img, [0 0 0; 1 max(0.01,handles.curlevel) 1]))
    end
else
    imshow((handles.img))
end


function contourlevel_Callback(hObject, eventdata, handles)
% hObject    handle to contourlevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contourlevel as text
%        str2double(get(hObject,'String')) returns contents of contourlevel as a double


% --- Executes during object creation, after setting all properties.
function contourlevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contourlevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mincellperim_Callback(hObject, eventdata, handles)
% hObject    handle to mincellperim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mincellperim as text
%        str2double(get(hObject,'String')) returns contents of mincellperim as a double


% --- Executes during object creation, after setting all properties.
function mincellperim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mincellperim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxcellperim_Callback(hObject, eventdata, handles)
% hObject    handle to maxcellperim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxcellperim as text
%        str2double(get(hObject,'String')) returns contents of maxcellperim as a double


% --- Executes during object creation, after setting all properties.
function maxcellperim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxcellperim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minarea_Callback(hObject, eventdata, handles)
% hObject    handle to minarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minarea as text
%        str2double(get(hObject,'String')) returns contents of minarea as a double


% --- Executes during object creation, after setting all properties.
function minarea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxarea_Callback(hObject, eventdata, handles)
% hObject    handle to maxarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxarea as text
%        str2double(get(hObject,'String')) returns contents of maxarea as a double


% --- Executes during object creation, after setting all properties.
function maxarea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxarea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in runcload.
function runcload_Callback(hObject, eventdata, handles)
% hObject    handle to runcload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.runcload,'Enable','Off');
pause(0.01)
data = cload(handles.filename,'contourlevel',paramval('contourlevel',handles),'mincellperim',paramval('mincellperim',handles),'maxcellperim',paramval('maxcellperim',handles),'minarea',paramval('minarea',handles),'maxarea',paramval('maxarea',handles),'optmid',logical(get(handles.optmid,'Value')),'filterwidth',logical(get(handles.filterwidth,'Value')),'wdfilter',logical(get(handles.wdfilter,'Value')),'npoints',paramval('npoints',handles),'maxthick',paramval('maxthick',handles),'minthick',paramval('minthick', handles));
if get(handles.filtint,'Value')
    data = filtint(data);
end
set(handles.runcload,'Enable','On');
set(handles.celldelete,'Enable','On');
assignin('base','data',data);
redisplay(handles);
cplot(data);
if get(handles.showmid,'Value')
    cplot(data,'r','mid');
end
set(handles.numcells, 'String', ['Number of cells: ' num2str(data.CellCount)]);
set(handles.runcload,'Enable','On');
handles.celldata = data;
guidata(hObject, handles);


% --- Return user-entered parameter value as number
function x = paramval(param,handles)
x = str2num(get(handles.(param),'String'));


% --- Executes on button press in optmid.
function optmid_Callback(hObject, eventdata, handles)
% hObject    handle to optmid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optmid


% --- Executes on button press in filterwidth.
function filterwidth_Callback(hObject, eventdata, handles)
% hObject    handle to filterwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filterwidth


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function filequit_Callback(hObject, eventdata, handles)
% hObject    handle to filequit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close


% --- Executes during object creation, after setting all properties.
function runcload_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runcload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in wdfilter.
function wdfilter_Callback(hObject, eventdata, handles)
% hObject    handle to wdfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wdfilter


% --- Executes on mouse press over axes background.
function mainfig_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to mainfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in celldelete.
function celldelete_Callback(hObject, eventdata, handles)
% hObject    handle to celldelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of celldelete

button = get(hObject,'Value');
while button == 1 && get(hObject,'Value') == 1
    [xi, yi ,button] = ginput(1);
    if button == 1
        for cellcount = 1:handles.celldata.CellCount
            if inpolygon(xi, yi,handles.celldata.cells(cellcount).border(1,:),handles.celldata.cells(cellcount).border(2,:));
                handles.celldata.cells(cellcount) = [];
                handles.celldata.CellCount = handles.celldata.CellCount - 1;
                redisplay(handles);
                cplot(handles.celldata);
                if get(handles.showmid,'Value')
                    cplot(handles.celldata,'r','mid');
                end
                assignin('base','data',handles.celldata);
                guidata(hObject, handles);
                set(handles.numcells, 'String', ['Number of cells: ' num2str(handles.celldata.CellCount)]);
                break
            end
        end
    end
end
set(hObject, 'Value',0);



function npoints_Callback(hObject, eventdata, handles)
% hObject    handle to npoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of npoints as text
%        str2double(get(hObject,'String')) returns contents of npoints as a double


% --- Executes during object creation, after setting all properties.
function npoints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to npoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in filtint.
function filtint_Callback(hObject, eventdata, handles)
% hObject    handle to filtint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filtint



function minthick_Callback(hObject, eventdata, handles)
% hObject    handle to minthick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minthick as text
%        str2double(get(hObject,'String')) returns contents of minthick as a double


% --- Executes during object creation, after setting all properties.
function minthick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minthick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxthick_Callback(hObject, eventdata, handles)
% hObject    handle to maxthick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxthick as text
%        str2double(get(hObject,'String')) returns contents of maxthick as a double


% --- Executes during object creation, after setting all properties.
function maxthick_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxthick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function curfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in showmid.
function showmid_Callback(hObject, eventdata, handles)
% hObject    handle to showmid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showmid

if isfield(handles,'celldata') == 1
    redisplay(handles);
    cplot(handles.celldata);
    if get(handles.showmid,'Value')
        cplot(handles.celldata,'r','mid');
    end
end
    


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tempimg = handles.img(:,:,1);
figure, hist(double(tempimg(:)),128)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tempimg = handles.img(:,:,1);
set(handles.pushbutton5,'Enable','Off');
pause(0.01)
figure, contour(tempimg(:,:,1)); axis ij; colorbar
set(handles.pushbutton5,'Enable','On');


% --- Executes on button press in zoombutton.
function zoombutton_Callback(hObject, eventdata, handles)
% hObject    handle to zoombutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoombutton
if get(hObject,'Value')
    zoom on
else
    zoom off
end