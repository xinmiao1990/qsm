function [iField,TE, delta_TE, matrix_size, voxel_size, B0_dir] = Read_Philips_enhanced_DICOM_xm(path,varargin)
% Written by:  Eamon Doyle, eamon@cornercase.net
%
% A function to read Dicom series for T2 mapping.  This is specifically
% meant to correct the pixel value scaling that is common in Philips
% dicoms.  Other manufacturers have not yet been considered.
%
p = inputParser;
addRequired(p,'path',@ischar);
addParameter(p,'Nx',480);
addParameter(p,'Ny',480);
addParameter(p,'Nz',92);
addParameter(p,'Nt',4);
addParameter(p,'adjustScale',1);
addParameter(p,'rescaleFieldName','realWorld');
addParameter(p,'dispImageAndPause',0);
parse(p,path,varargin{:});

Nx = p.Results.Nx;
Ny = p.Results.Ny;
Nz = p.Results.Nz;
Nt = p.Results.Nt;


if isstruct(path)
    dcmInfo = path;
else
    dcmInfo = dicominfo(path);
end

if isfield( dcmInfo,'Manufacturer')
    if strfind(upper(dcmInfo.Manufacturer),'PHILIPS')
        man = 'PHILIPS';
    end
end

%pick extended or not
switch dcmInfo.MediaStorageSOPClassUID
    case '1.2.840.10008.5.1.4.1.1.4'
        %MR Image Storage	 
        itype = 'standard';
    case '1.2.840.10008.5.1.4.1.1.4.1'
        %Enhanced MR Image Storage	 
        itype = 'enhanced';
    case '1.2.840.10008.5.1.4.1.1.4.2'
        %MR Spectroscopy Storage
        itype = 'spectro';
        error('spectroscopy images not yet supported');
    otherwise
        error('unknown image type');
end

%pick field type
switch p.Results.rescaleFieldName
    case 'private'
        fieldType = 'private';
    case 'realWorld'
        fieldType = 'realWorld';
    otherwise
        error('unknown field name');
end

scaleSlope = zeros(1,Nz*Nt*2);
scaleIntercept = zeros(1,Nz*Nt*2);
TEs = zeros(1,Nz*Nt*2);
% all of this is to set the scale factor correctly
if p.Results.adjustScale
    if strcmp(man,'PHILIPS')
        if  strcmp(itype,'enhanced')               
            % TE = dcmInfo.PerFrameFunctionalGroupsSequence.Item_1.Private_2005_140f.Item_1.EchoTime/1000;
            % XM 11/30/2016: can't find TE information.
            minSlice = 1e10;
            maxSlice = -1e10;
            for i = 1:(Nz*Nt*2)
                scaleInfo = ... 
                    eval(sprintf('dcmInfo.PerFrameFunctionalGroupsSequence.Item_%i.RealWorldValueMappingSequence.Item_1',i));
                scaleSlope(i) = scaleInfo.RealWorldValueSlope;
                scaleIntercept(i) = scaleInfo.RealWorldValueIntercept;
                
                sliceInfo = ... 
                    eval(sprintf('dcmInfo.PerFrameFunctionalGroupsSequence.Item_%i.PlanePositionSequence.Item_1.ImagePositionPatient',i));
                slicePos = sliceInfo(3);
                TEs(i) = ...
                    eval(sprintf('dcmInfo.PerFrameFunctionalGroupsSequence.Item_%i.MREchoSequence.Item_1.EffectiveEchoTime',i)); 
                
                if slicePos<minSlice
                    minSlice = slicePos;
                    minLoc = sliceInfo;
                end
                if slicePos>maxSlice
                    maxSlice = slicePos;
                    maxLoc = sliceInfo;
                end
                
            end
            
            pixelInfo = ... 
                    eval(sprintf('dcmInfo.PerFrameFunctionalGroupsSequence.Item_%i.PixelMeasuresSequence.Item_1',i));
            voxel_size(1)=pixelInfo.PixelSpacing(1);
            voxel_size(2)=pixelInfo.PixelSpacing(2);
            voxel_size(3)=pixelInfo.SliceThickness;  
            
            matrix_size(1) = Nx;
            matrix_size(2) = Ny;
            matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3)) + 1;
            
            Affine2D = dcmInfo.PerFrameFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient;              
            Affine2D = reshape(Affine2D,[3 2]);
            Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
            B0_dir = Affine3D\[0 0 1]';
        end
        
    end
end

TE = unique(TEs);
delta_TE = TE(2)-TE(1);

images = dicomread(dcmInfo);
images = squeeze(images);
im = zeros(size(images));

if p.Results.adjustScale
    if strcmp(fieldType, 'private')
        im=((double(images))-scaleIntercept)/scaleSlope;
    elseif strcmp(fieldType, 'realWorld')
        for i = 1:(Nz*Nt*2)
            im(:,:,i)=1e-3*(double(images(:,:,i))*scaleSlope(i)+scaleIntercept(i)) ;
        end
    else
        error('something very bad happened');
    end
end

im = reshape(im,[Nx,Ny,Nt,Nz,2]);
    
if p.Results.dispImageAndPause == 1
    imagesc(im);
    title(dcmInfo.Filename);
    pause;
    disp('press key to continue');
    
im = permute(im,[1 2 4 3 5]);
iField = im(:,:,:,:,1).*exp(-1i*im(:,:,:,:,2));
end
