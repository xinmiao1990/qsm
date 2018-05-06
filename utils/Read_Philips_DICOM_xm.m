% for PHILIPS DICOM images consisting of magnitude and phase parts
%   [iField,voxel_size,matrix_size,CF,delta_TE,B0_dir] =Read_Philips_DICOM(DicomFolder);
%
%   output
%   iField - the multi-echo complex MR image
%   voxel_size - size of the voxel in mm
%   matrix_size - dimension of the field of view
%   CF - central frequency in Hz
%   delta_TE - TE spacing in s
%   TE - echo times in s
%   B0_dir - direction of the B0 field
%
%   input
%   DicomFloder: This folder should NOT contain anything else 
%                except the DICOM images belonging to a particular
%                series that you want to process.
%
%   Apdated from Tian Liu
%   Created by Shuai Wang on 2013.07.19 
%   Last modified by Tian Liu on 2013.08.06


function [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_Philips_DICOM_xm(DicomFolder)
    
    filelist = dir(DicomFolder);
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
 
    fnTemp=[DicomFolder '/' filelist(10).name];
    info = dicominfo(fnTemp);
    matrix_size(1) = single(info.Width);
    matrix_size(2) = single(info.Height);
    voxel_size(1,1) = single(info.PixelSpacing(1));
    voxel_size(2,1) = single(info.PixelSpacing(2));
    voxel_size(3,1) = single(info.SliceThickness);
    CF = info.ImagingFrequency *1e6;
    
    minSlice = 1e10;
    maxSlice = -1e10;
    NumEcho = 0;
    
    for i = 1:length(filelist)
        info = dicominfo([DicomFolder '/' filelist(i).name]);
        
        if info.SliceLocation<minSlice
            minSlice = info.SliceLocation;
            minLoc = info.ImagePositionPatient;
        end
        if info.SliceLocation>maxSlice
            maxSlice = info.SliceLocation;
            maxLoc = info.ImagePositionPatient;
        end
        if info.EchoNumber>NumEcho
            NumEcho = info.EchoNumber;
        end
    end

    matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3)) + 1;  
    Affine2D = reshape(info.ImageOrientationPatient,[3 2]);
    Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
    B0_dir = Affine3D\[0 0 1]';

    iMag = single(zeros([matrix_size NumEcho]));
    iPhase = single(zeros([matrix_size NumEcho]));
    TE = single(zeros([NumEcho 1]));


    for i = 1:length(filelist)

        info = dicominfo([DicomFolder '/' filelist(i).name]);
        disp([DicomFolder '/' filelist(i).name]);
        

        slice = int32(round(norm(info.ImagePositionPatient-minLoc)/voxel_size(3)) +1);
        
        
        if TE(info.EchoNumber)==0
            TE(info.EchoNumber)=info.EchoTime*1e-3;
        end
        
        if (info.ImageType(18)=='P')|(info.ImageType(18)=='p') % the last half of dicom files contain phase information
            
            temp  = single(dicomread([DicomFolder '/' filelist(i).name])');%phase
            
            if (size(temp,1)==0)
                    iPhase(:,:,slice,info.EchoNumber) =0;
                    s = sprintf('i=%d, slice=%d, echonumber=%d, TE=%f, sizedata=%d,%d, info.ImageType(18)=%s',i, slice, info.EchoNumber, info.EchoTime*1e-3,size(temp,1),size(temp,2),info.ImageType(18));
                    disp(s);
            else
                    %iPhase(:,:,slice,info.EchoNumber)  = 1e-3*(temp*info.RealWorldValueMappingSequence.Item_1.RealWorldValueSlope+info.RealWorldValueMappingSequence.Item_1.RealWorldValueIntercept);%phase
                    iPhase(:,:,slice,info.EchoNumber)  = 1e-3*(temp*1.535-3142);
            end


        elseif (info.ImageType(18)=='M')|(info.ImageType(18)=='m') % the first half of dicom files contain magnitude information
            
            temp  = single(dicomread([DicomFolder '/' filelist(i).name])');%magnitude
                        
            if (size(temp,1)==0)
                    iMag(:,:,slice,info.EchoNumber) =1;
                    s = sprintf('i=%d, slice=%d, echonumber=%d, TE=%f, sizedata=%d,%d, info.ImageType(18)=%s',i, slice, info.EchoNumber, info.EchoTime*1e-3,size(temp,1),size(temp,2),info.ImageType(18));
                    disp(s);
            else
                    iMag(:,:,slice,info.EchoNumber) = temp;              
            end

            
        end
    end
    
    iField = iMag.*exp(-1i*iPhase);
    clear iMag iPhase;
    if length(TE)==1
        delta_TE = TE;
    else
        delta_TE = TE(2) - TE(1);
    end
    
end