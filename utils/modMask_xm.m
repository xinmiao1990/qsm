% Mask generation
%   Mask = genMask(iField, voxel_size)
% 
%   output
%   Mask - the biggest contiguous region that has decent SNR
%
%   input
%   iField - the complex MR image
%   voxel_size - the size of a voxel
%
%   Created by Tian Liu in 20013.07.24
%   Last modified by Tian Liu on 2013.07.24

function Mask = modMask_xm(mask, voxel_size,eroIn,eroOut)
    matrix_size = size(mask);
    m1 = SMV(mask,matrix_size, voxel_size, eroIn)>0.999;   % erode the boundary by 10mm
    l = bwlabeln(m1,6);                       % find the biggest contiguous region
    Mask = (l==mode(l(l~=0)));
    Mask = SMV(Mask, matrix_size, voxel_size, eroOut)>0.001; % restore the enrosion by 4m

end
