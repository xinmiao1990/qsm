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

% Last Modified by GUIDE v2.5 18-May-2018 10:07:29

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
path1 = getenv('PATH')
path1 = [path1 ':/Users/xinmiao/Documents/softwares/mricron']
setenv('PATH', path1)
path1 = getenv('PATH')
path1 = [path1 ':/usr/local/fsl/bin']
setenv('PATH', path1)
setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
handles.iField = [];
handles.voxel_size = [];
handles.matrix_size = [];
handles.CF = [];
handles.delta_TE = [];
handles.B0_dir = [];
handles.TE = [];
handles.bigMask = [];
handles.iFreq_raw = [];
handles.N_std = [];
handles.iFreq = [];
handles.RDF = [];
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = QSM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in choseBaseDir.
function choseBaseDir_Callback(hObject, eventdata, handles)
% hObject    handle to choseBaseDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
baseDCMdir = uigetdir(matlabroot,'Where are the dicoms?');
set(handles.edit_baseDCMdir, 'String', baseDCMdir);
baseNIIdir = uigetdir(matlabroot,'Where are the nifties?');
set(handles.edit_baseNIIdir, 'String', baseNIIdir);
baseMATdir = uigetdir(matlabroot,'Where are the mat files?');
set(handles.edit_baseMATdir, 'String', baseMATdir);

% Update handles.structure
handles.baseDCMdir = baseDCMdir;
handles.baseNIIdir = baseNIIdir;
handles.baseMATdir = baseMATdir;
guidata(hObject, handles); 


% --- Executes on button press in convertDCM2NII.
function convertDCM2NII_Callback(hObject, eventdata, handles)
% hObject    handle to convertDCM2NII (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
te1DCMdir = [handles.baseDCMdir '/echo1'];
cmd = ['dcm2nii -n y -o ' handles.baseNIIdir ' ' te1DCMdir];
[status,~] = system(cmd);
if(status==0)
    magNIIfile = [handles.baseNIIdir '/qsmMag.nii.gz'];
    result = dir([handles.baseNIIdir '/o*.nii.gz']);
    movefile([handles.baseNIIdir '/' result.name], magNIIfile);
    set(handles.edit_magNIIfile, 'String', magNIIfile);
else
    set(handles.edit_magNIIfile, 'String', 'Fail!');
end

% Update handles.structure
guidata(hObject, handles); 

% --- Executes on button press in betBrainExtraction.
function betBrainExtraction_Callback(hObject, eventdata, handles)
% hObject    handle to betBrainExtraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
magNIIfile = [handles.baseNIIdir '/qsmMag.nii.gz'];
maskNIIfile = [handles.baseNIIdir '/qsm_mask.nii.gz'];
cmd = ['bet ' magNIIfile ' ' handles.baseNIIdir '/qsm.nii.gz -m'];
[status,~] = system(cmd);
if(status==0)
    temp = load_untouch_nii_gz(maskNIIfile);
    temp= logical(temp.img);
    Mask = flip(flip(temp,2),1);
    set(handles.edit_maskNIIfile, 'String', maskNIIfile);
    
    maskMATfile = [handles.baseMATdir '/bigMask.mat'];
    save(maskMATfile,'Mask');
    set(handles.edit_maskMATfile, 'String', maskMATfile);
else
    set(handles.edit_maskNIIfile, 'String', 'Fail!');
end

% Update handles structure
handles.bigMask = Mask;
guidata(hObject, handles);


% --- Executes on button press in convertDCM2MAT.
function convertDCM2MAT_Callback(hObject, eventdata, handles)
% hObject    handle to convertDCM2MAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
allDCMdir = [handles.baseDCMdir '/all'];
rawMATfile = [handles.baseMATdir '/raw.mat'];
[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_Philips_DICOM_xm(allDCMdir);
% load('/Users/xinmiao/Documents/research/data/QSM_raw_mat/SCD053/raw.mat');
iField = flip(iField,3);
B0_dir(3) = -B0_dir(3);
view3dgui(abs(iField(:,:,:,1)), voxel_size);
keyboard;

save(rawMATfile,'iField','voxel_size','matrix_size','CF','delta_TE','TE','B0_dir');
set(handles.edit_rawMATfile, 'String', rawMATfile);

% Update handles structure
handles.iField = iField;
handles.voxel_size = voxel_size;
handles.matrix_size = matrix_size;
handles.CF = CF;
handles.delta_TE = delta_TE;
handles.B0_dir = B0_dir;
handles.TE = TE;
guidata(hObject,handles);




% --- Executes on button press in eddyCurrentCorrection.
function eddyCurrentCorrection_Callback(hObject, eventdata, handles)
% hObject    handle to eddyCurrentCorrection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rawMATfile = [handles.baseMATdir '/raw.mat'];
load(rawMATfile);
iField1 = eddy_current_cor_hyper(iField);
keyboard;
iField = iField1; clear iField1;
save(rawMATfile,'iField','voxel_size','matrix_size','CF','delta_TE','TE','B0_dir');
set(handles.edit_rawMATfile3, 'String', 'Saved Corrected version!');
% Update handles structure
guidata(hObject,handles);


% --- Executes on button press in fitB0Field.
function fitB0Field_Callback(hObject, eventdata, handles)
% hObject    handle to fitB0Field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.iField)
    rawMATfile = [handles.baseMATdir '/raw.mat'];
    load(rawMATfile);
    handles.iField = iField;
    handles.voxel_size = voxel_size;
    handles.matrix_size = matrix_size;
    handles.CF = CF;
    handles.delta_TE = delta_TE;
    handles.B0_dir = B0_dir;
    handles.TE = TE;
    guidata(hObject,handles);
end

if isempty(handles.bigMask)
    bigMaskFile = [handles.baseMATdir '/bigMask.mat'];
    load(bigMaskFile);
    handles.bigMask = Mask;
    guidata(hObject,handles);
end

Mask2 = genMask_xm(abs(handles.iField(:,:,:,1)), handles.voxel_size,10,8); 
Mask = logical(handles.bigMask.*Mask2);
view3dgui(abs(handles.iField(:,:,:,1)).*Mask, handles.voxel_size);
keyboard;
handles.bigMask = Mask;
bigMaskFile = [handles.baseMATdir '/bigMask.mat'];
save(bigMaskFile,'Mask');

[iFreq_raw, N_std] = Fit_ppm_complex(handles.iField.*repmat(Mask,[1 1 1 4]));
view3dgui(iFreq_raw, handles.voxel_size);
view3dgui(N_std, handles.voxel_size);
keyboard;
% Update handles structure
ifreqrawMATfile = [handles.baseMATdir '/iFreq_raw.mat'];
save(ifreqrawMATfile,'iFreq_raw', 'N_std');
set(handles.edit_ifreqrawMATfile, 'String', ifreqrawMATfile);

handles.iFreq_raw = iFreq_raw;
handles.N_std = N_std;
guidata(hObject,handles);


% --- Executes on button press in unwrapPhase.
function unwrapPhase_Callback(hObject, eventdata, handles)
% hObject    handle to unwrapPhase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.iFreq_raw)
    ifreqrawMATfile = [handles.baseMATdir '/iFreq_raw.mat'];
    load(ifreqrawMATfile);
    handles.iFreq_raw = iFreq_raw;
    guidata(hObject,handles);
end

if isempty(handles.bigMask)
    bigMaskFile = [handles.baseMATdir '/bigMask.mat'];
    load(bigMaskFile);
    handles.bigMask = Mask;
    guidata(hObject,handles);
end

[iFreq] = unwrapLaplacian(handles.iFreq_raw.*handles.bigMask, handles.matrix_size, handles.voxel_size);
view3dgui(iFreq, handles.voxel_size);

ifreqMATfile = [handles.baseMATdir '/iFreq.mat'];
save(ifreqMATfile,'iFreq');
set(handles.edit_ifreqMATfile,'String', ifreqMATfile);
% Update handles structure
handles.iFreq = iFreq; clear iFreq;
guidata(hObject,handles);


% --- Executes on button press in extractNoiseMask.
function extractNoiseMask_Callback(hObject, eventdata, handles)
% hObject    handle to extractNoiseMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.iFreq_raw)
    ifreqrawMATfile = [handles.baseMATdir '/iFreq_raw.mat'];
    load(ifreqrawMATfile);
    handles.iFreq_raw = iFreq_raw;
    handles.N_std = N_std;
    guidata(hObject,handles);
end

if isempty(handles.bigMask)
    bigMaskFile = [handles.baseMATdir '/bigMask.mat'];
    load(bigMaskFile);
    handles.bigMask = Mask;
    guidata(hObject,handles);
end

threshold = 0.002;
noiseMask = handles.bigMask;
noiseMask(handles.N_std > threshold) = 0;
noiseMask = logical(noiseMask);
Mask = noiseMask;
view3dgui(~Mask,[0.45 0.45 1.3]);

noiseMaskFile = [handles.baseMATdir '/noiseMask.mat'];
save(noiseMaskFile, 'Mask');clear Mask;
set(handles.edit_edit_noisemaskMATfile,'String', noiseMaskFile);

% Update handles structure
handles.noiseMask = noiseMask; clear noiseMask;
guidata(hObject,handles);


% --- Executes on button press in removeBGField.
function removeBGField_Callback(hObject, eventdata, handles)
% hObject    handle to removeBGField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.iField)
    rawMATfile = [handles.baseMATdir '/raw.mat'];
    load(rawMATfile);
    handles.iField = iField;
    handles.voxel_size = voxel_size;
    handles.matrix_size = matrix_size;
    handles.CF = CF;
    handles.delta_TE = delta_TE;
    handles.B0_dir = B0_dir;
    handles.TE = TE;
    guidata(hObject,handles);
end

if isempty(handles.bigMask)
    bigMaskFile = [handles.baseMATdir '/bigMask.mat'];
    load(bigMaskFile);
    handles.bigMask = Mask;
    guidata(hObject,handles);
end

if isempty(handles.iFreq_raw) || isempty(handles.N_std)
    ifreqrawMATfile = [handles.baseMATdir '/iFreq_raw.mat'];
    load(ifreqrawMATfile);
    handles.iFreq_raw = iFreq_raw;
    handles.N_std = N_std;
    guidata(hObject,handles);
end

if isempty(handles.iFreq)
    ifreqMATfile = [handles.baseMATdir '/iFreq.mat'];
    load(ifreqMATfile);
    handles.iFreq = iFreq;
    guidata(hObject,handles);
end

[RDF,~] = PDF(handles.iFreq, handles.N_std, handles.bigMask,...
    handles.matrix_size, handles.voxel_size, handles.B0_dir);
view3dgui(RDF, handles.voxel_size);
keyboard;

RDFFile = [handles.baseMATdir '/RDF.mat'];
save(RDFFile, 'RDF');
set(handles.edit_rdfMATfile,'String',RDFFile);

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

QSML1File = [handles.baseMATdir '/QSML1.mat'];
save(QSML1File, 'QSML1');
set(handles.edit_qsmMATfile,'String',QSML1File);

% Update handles structure
handles.QSML1 = QSML1; clear QSML1;
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
function edit_baseDCMdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_baseDCMdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_te1DCMdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_te1DCMdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_allDCMdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_allDCMdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_rawMATfile2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rawMATfile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes during object creation, after setting all properties.
function edit_magNIIfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_magNIIfile (see GCBO)
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



function edit_baseNIIdir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_baseNIIdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_baseNIIdir as text
%        str2double(get(hObject,'String')) returns contents of edit_baseNIIdir as a double


% --- Executes during object creation, after setting all properties.
function edit_baseNIIdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_baseNIIdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_baseMATdir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_baseMATdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_baseMATdir as text
%        str2double(get(hObject,'String')) returns contents of edit_baseMATdir as a double


% --- Executes during object creation, after setting all properties.
function edit_baseMATdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_baseMATdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_baseDCMdir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_baseDCMdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_baseDCMdir as text
%        str2double(get(hObject,'String')) returns contents of edit_baseDCMdir as a double



function edit_te1DCMdir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_te1DCMdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_te1DCMdir as text
%        str2double(get(hObject,'String')) returns contents of edit_te1DCMdir as a double



function edit_allDCMdir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_allDCMdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_allDCMdir as text
%        str2double(get(hObject,'String')) returns contents of edit_allDCMdir as a double



function edit_magNIIfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_magNIIfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_magNIIfile as text
%        str2double(get(hObject,'String')) returns contents of edit_magNIIfile as a double



function edit_magNIIfile2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_magNIIfile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_magNIIfile2 as text
%        str2double(get(hObject,'String')) returns contents of edit_magNIIfile2 as a double


% --- Executes during object creation, after setting all properties.
function edit_magNIIfile2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_magNIIfile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maskNIIfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maskNIIfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maskNIIfile as text
%        str2double(get(hObject,'String')) returns contents of edit_maskNIIfile as a double


% --- Executes during object creation, after setting all properties.
function edit_maskNIIfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maskNIIfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rawMATfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rawMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rawMATfile as text
%        str2double(get(hObject,'String')) returns contents of edit_rawMATfile as a double


% --- Executes during object creation, after setting all properties.
function edit_rawMATfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rawMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rawMATfile2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rawMATfile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rawMATfile2 as text
%        str2double(get(hObject,'String')) returns contents of edit_rawMATfile2 as a double



function edit_rawMATfile3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rawMATfile3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rawMATfile3 as text
%        str2double(get(hObject,'String')) returns contents of edit_rawMATfile3 as a double


% --- Executes during object creation, after setting all properties.
function edit_rawMATfile3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rawMATfile3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_ifreqrawMATfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ifreqrawMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ifreqrawMATfile as text
%        str2double(get(hObject,'String')) returns contents of edit_ifreqrawMATfile as a double


% --- Executes during object creation, after setting all properties.
function edit_ifreqrawMATfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ifreqrawMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ifreqMATfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ifreqMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ifreqMATfile as text
%        str2double(get(hObject,'String')) returns contents of edit_ifreqMATfile as a double


% --- Executes during object creation, after setting all properties.
function edit_ifreqMATfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ifreqMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maskMATfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maskMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maskMATfile as text
%        str2double(get(hObject,'String')) returns contents of edit_maskMATfile as a double


% --- Executes during object creation, after setting all properties.
function edit_maskMATfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maskMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_noisemaskMATfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_noisemaskMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_noisemaskMATfile as text
%        str2double(get(hObject,'String')) returns contents of edit_noisemaskMATfile as a double


% --- Executes during object creation, after setting all properties.
function edit_noisemaskMATfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_noisemaskMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_RDFMATfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RDFMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_RDFMATfile as text
%        str2double(get(hObject,'String')) returns contents of edit_RDFMATfile as a double


% --- Executes during object creation, after setting all properties.
function edit_RDFMATfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RDFMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_qsmMATfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_qsmMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_qsmMATfile as text
%        str2double(get(hObject,'String')) returns contents of edit_qsmMATfile as a double


% --- Executes during object creation, after setting all properties.
function edit_qsmMATfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_qsmMATfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
