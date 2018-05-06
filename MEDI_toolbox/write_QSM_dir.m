function y = write_QSM_dir(QSM,dicomDir,saveDir)
% if ~exist('qsm_DICOM','dir')
%     mkdir('qsm_DICOM')
% end
qsm1 = int16(QSM*1000);
fileList = dir(dicomDir);
sliceNum = size(QSM,3);
sliceIdx = 0;
imIdx = 1;
UID = dicomuid;
saveDir = [saveDir '/' 'DICOM'];
mkdir(saveDir);
while sliceIdx < sliceNum
    fid = fopen([dicomDir '/' fileList(imIdx).name]);
    if fid>0
        sliceIdx=sliceIdx + 1;
        fclose(fid);
        info = dicominfo([dicomDir '/' fileList(imIdx).name]);
        %--
        info.SeriesDescription = 'QSM recon';
        info.SeriesInstanceUID = UID;
        info.SeriesNumber = 99;
        info.EchoNumber = 1;
        info.EchoTime = 0.0;
        %--
        dicomwrite(qsm1(:,:,sliceIdx),[saveDir '/' num2str(sliceIdx) '.dcm'], info);
    end
    imIdx = imIdx + 1;
end
addpath(saveDir);