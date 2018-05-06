function varargout = QSM(varargin)
% QSM MATLAB code for QSM.fig
%      QSM, by itself, creates a new QSM or raises the existing
%      singleton*.
%
%      H = QSM returns the handle to a new QSM or the handle to
%      the existing singleton*.
%
%      QSM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QSM.M with the given input arguments.
%
%      QSM('Property','Value',...) creates a new QSM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before QSM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to QSM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help QSM

% Last Modified by GUIDE v2.5 06-May-2018 14:55:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QSM_OpeningFcn, ...
                   'gui_OutputFcn',  @QSM_OutputFcn, ...
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


% --- Executes just before QSM is made visible.
function QSM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QSM (see VARARGIN)

% Choose default command line output for QSM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes QSM wait for user response (see UIRESUME)
% uiwait(handles.figure1);
addpath('utils')
addpath('MEDI_toolbox');
addpath('view3dgui');



% --- Outputs from this function are returned to the command line.
function varargout = QSM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes on button press in convertDCM2NII.
function convertDCM2NII_Callback(hObject, eventdata, handles)
% hObject    handle to convertDCM2NII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
te1DCMdir = uigetdir(matlabroot,'Where are the first echo dicoms?');
set(handles.edit1, 'String', te1DCMdir);
allNIIdir = uigetdir(matlabroot,'Where are the NII files?');
cmd = ['dcm2nii -n y -o ' allNIIdir ' ' te1DCMdir];
[status,~] = system(cmd);
if(status==0)
    magNIIfile = [allNIIdir '/qsmMag.nii.gz'];
    result = dir([allNIIdir '/o*.nii.gz']);
    movefile([allNIIdir '/' result.name], magNIIfile);
    set(handles.edit2, 'String', magNIIfile);
else
    set(handles.edit2, 'String', 'Fail!');
    magNIIfile = [];
end
% Update handles.structure
handles.te1DCMdir = te1DCMdir;
handles.allNIIdir = allNIIdir;
handles.magNIIfile = magNIIfile;
guidata(hObject, handles); 

% --- Executes on button press in betBrainExtraction.
function betBrainExtraction_Callback(hObject, eventdata, handles)
% hObject    handle to betBrainExtraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
maskNIIfile = [handles.allNIIdir '/qsm_mask.nii.gz'];
cmd = ['bet ' handles.magNIIfile ' ' handles.allNIIdir '/qsm.nii.gz -m'];
[status,~] = system(cmd);
if(status==0)
    temp = load_untouch_nii_gz(maskNIIfile);
    temp= logical(temp.img);
    handles.bigMask = flip(temp,2);
    set(handles.edit5, 'String', maskNIIfile);
else
    set(handles.edit5, 'String', 'Fail!');
    maskNIIfile = [];
end
% Update handles structure
handles.maskNIIfile = maskNIIfile;
guidata(hObject, handles);


% --- Executes on button press in convertDCM2MAT.
function convertDCM2MAT_Callback(hObject, eventdata, handles)
% hObject    handle to convertDCM2MAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
allDCMdir = uigetdir(matlabroot,'Where are the dicoms?');
set(handles.edit3, 'String', allDCMdir);
allMATdir = uigetdir(matlabroot,'Where are the MAT files?');
rawMATfile = [allMATdir '/raw.mat'];
% [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_Philips_DICOM_xm(allDCMdir);
load('/Users/xinmiao/Documents/research/data/QSM_raw_mat/SCD053/raw.mat');
save(rawMATfile,'iField','voxel_size','matrix_size','CF','delta_TE','TE','B0_dir');
set(handles.edit4, 'String', rawMATfile);

Mask = genMask_xm(abs(iField(:,:,:,1)), voxel_size,10,8); 
handles.bigMask = handles.bigMask .* Mask;

% Update handles structure
handles.iField = iField; clear iField;
handles.voxel_size = voxel_size; clear voxel_size;
handles.matrix_size = matrix_size; clear matrix_size;
handles.CF = CF; clear CF;
handles.delta_TE = delta_TE; clear delta_TE;
handles.TE = TE; clear TE;
handles.B0_dir = B0_dir; clear B0_dir;
handles.allDCMdir = allDCMdir;
handles.allMATdir = allMATdir;
guidata(hObject,handles);

view3dgui(abs(handles.iField(:,:,:,1)).*handles.bigMask, handles.voxel_size);
keyboard;
Mask = handles.bigMask;
bigMaskFile = [handles.allMATdir '/bigMask.mat'];
save(bigMaskFile, 'Mask');clear Mask;


% --- Executes on button press in eddyCurrentCorrection.
function eddyCurrentCorrection_Callback(hObject, eventdata, handles)
% hObject    handle to eddyCurrentCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load([handles.allMATdir '/raw.mat']);
iField1 = eddy_current_cor_hyper(handles.iField);
keyboard;
iField = iField1; clear iField1;
save([handles.allMATdir '/raw.mat'],'iField','voxel_size','matrix_size','CF','delta_TE','TE','B0_dir');

% Update handles structure
handles.iField = iField; clear iField;
guidata(hObject,handles);


% --- Executes on button press in fitB0Field.
function fitB0Field_Callback(hObject, eventdata, handles)
% hObject    handle to fitB0Field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[iFreq_raw, N_std] = Fit_ppm_complex(handles.iField.*repmat(handles.bigMask,[1 1 1 4]));
view3dgui(iFreq_raw, handles.voxel_size);
view3dgui(N_std, handles.voxel_size);

% Update handles structure
handles.iFreq_raw = iFreq_raw; clear iFreq_raw;
handles.N_std = N_std; clear N_std;
guidata(hObject,handles);


% --- Executes on button press in unwrapPhase.
function unwrapPhase_Callback(hObject, eventdata, handles)
% hObject    handle to unwrapPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[iFreq] = unwrapLaplacian(handles.iFreq_raw.*handles.bigMask, handles.matrix_size, handles.voxel_size);
view3dgui(iFreq, handles.voxel_size);

% Update handles structure
handles.iFreq = iFreq; clear iFreq;
guidata(hObject,handles);


% --- Executes on button press in extractNoiseMask.
function extractNoiseMask_Callback(hObject, eventdata, handles)
% hObject    handle to extractNoiseMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
threshold = 0.02;
noiseMask = handles.bigMask;
noiseMask(handles.N_std > threshold) = 0;
noiseMask = logical(noiseMask);
Mask = noiseMask;

noiseMaskFile = [handles.allMATdir '/noiseMask.mat'];
save(noiseMaskFile, 'Mask');clear Mask;

% Update handles structure
handles.noiseMask = noiseMask; clear noiseMask;
guidata(hObject,handles);


% --- Executes on button press in removeBGField.
function removeBGField_Callback(hObject, eventdata, handles)
% hObject    handle to removeBGField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[RDF,~] = PDF(handles.iFreq, handles.N_std, handles.bigMask,...
    handles.matrix_size, handles.voxel_size, handles.B0_dir);
view3dgui(RDF, handles.voxel_size);
keyboard;

% Update handles structure
handles.RDF = RDF; clear RDF;
guidata(hObject,handles);

% --- Executes on button press in inversion.
function inversion_Callback(hObject, eventdata, handles)
% hObject    handle to inversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nfm = handles.RDF / (2*pi * handles.delta_TE* handles.CF) * 1e6;
lambda = 4e-4;
QSML1= admm_qsm(nfm, handles.matrix_size, handles.noiseMask, handles.voxel_size, lambda);

view3dgui(QSML1, handles.voxel_size);
keyboard;

QSML1File = [handles.allMATdir '/QSML1.mat'];
save(QSML1File, 'QSML1');

% Update handles structure
handles.QSML1 = QSML1; clear QSML1;
handles.QSML1File = QSML1File;
guidata(hObject,handles);

% --- Executes on button press in registerT12QSM.
function registerT12QSM_Callback(hObject, eventdata, handles)
% hObject    handle to registerT12QSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bse = '/Applications/BrainSuite16a1/bin/bse16a1_x86_64-apple-darwin14';
addpath('/Users/xinmiao/Google Drive/matlab/projects/QSM_tissue/brainsuite_utilities/register_files_affine_v12');

dir = handles.allNIIdir;
foQSM = [dir '/qsmMag.nii.gz'];
fbfc = [dir '/SCD053.bfc.nii.gz'];
fmask = [dir '/SCD053.mask.nii.gz'];
flabels = [dir '/SCD053.svreg.label.nii.gz'];
patientType = 'SCD';

fqsm_mask = skullStripMask_MagnQSM(foQSM, bse, patientType);
[~, fregt1_mat] = regT1_to_MagnQSM(fbfc, fmask, foQSM, fqsm_mask, patientType);
flabel_trans = [dir '/QSM_coord.sverg.label.nii.gz'];
transSegmentLabels_to_QSMcoord(flabels, fbfc, foQSM, fregt1_mat, flabel_trans);


% --- Executes on button press in analyzeROI.
function analyzeROI_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.mat');
roiMaskFile = fullfile(path,file);
load(roiMaskFile);
roi_mean = mean(handles.QSML1(mask(:)));
roi_std = std(handles.QSML1(mask(:)));
set(handles.edit9, 'String', num2str(roi_mean));
set(handles.edit10, 'String', num2str(roi_std));


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



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


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



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


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



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


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



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


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
