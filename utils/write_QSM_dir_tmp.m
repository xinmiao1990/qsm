%   Last modified by Alexey Dimov on 2016.05.13

function y = write_QSM_dir_tmp(QSM,dicomDir,saveDir, SeriesDescription, SeriesNumber)
warning( 'off', 'all' );

% if ~exist('qsm_DICOM','dir')
%     mkdir('qsm_DICOM')
% end
QSM = permute(QSM,[2,1,3]); %change row/column order due to differences in representations in DICOM and Matlab
qsm1 = int16(QSM*1000);

filelist = dir(dicomDir);
i=1;
while i<=length(filelist)
    if filelist(i).isdir==1
        filelist = filelist([1:i-1 i+1:end]);   % skip folders
    elseif isempty(strfind(filelist(i).name,'.dcm'))
        filelist = filelist([1:i-1 i+1:end]);   % skip other files
    else
        i=i+1;
    end
end
    
sliceNum = size(QSM,3);
sliceIdx = 0;
imIdx = 1;
UID = dicomuid;

while sliceIdx < sliceNum
    fid = fopen([dicomDir '/' filelist(imIdx).name]);
    if fid>0
        sliceIdx=sliceIdx + 1;
        fclose(fid);
        info = dicominfo([dicomDir '/' filelist(imIdx).name]);
        %--
        info.SeriesDescription = SeriesDescription;
        info.SeriesInstanceUID = UID;
        info.SeriesNumber = SeriesNumber;
        info.EchoNumber = 1;
        disp(['TE=' num2str(info.EchoTime)]);
        info.EchoTime = 0.0;
        info.InstanceNumber = sliceIdx;
        %--
        dicomwrite(qsm1(:,:,sliceIdx),[saveDir '/' num2str(sliceIdx) '.dcm'], info);
    end
    imIdx = imIdx + 1;
end
addpath(saveDir);
warning( 'on', 'all' )