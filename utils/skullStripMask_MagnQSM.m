function fqsm_mask = skullStripMask_MagnQSM(foQSM, bse, patientType)

    [pathstr,name,~] = fileparts(foQSM);
    
    switch patientType
        case 'ACTL'
            patientID = name(1:7);
        case 'QCTL'
            patientID = name(1:7);
        otherwise
            patientID = name(1:6);
    end
    
    fqsm_mask = [pathstr '/' patientID '.qsm.mask.nii.gz'];
    if ~exist(fqsm_mask,'file')
        system([bse, ' -i ', foQSM, ' --mask ',fqsm_mask])
    else
        disp(['Found' fqsm_mask]);
    end

end