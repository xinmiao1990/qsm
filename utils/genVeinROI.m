function mask = genVeinROI(roimaskDir)

[mask,~]=ReadData3D([roimaskDir 'mask_1.mhd']);
mask = (mask~=0);

end