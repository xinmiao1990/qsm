function [ft1reg, fregt1_mat] = regT1_to_MagnQSM(fbfc, fmask, foQSM, fqsm_mask, patientType)

    [pathstr,name,~] = fileparts(foQSM); 
    
    switch patientType
        case 'ACTL'
            patientID = name(1:7);
        case 'QCTL'
            patientID = name(1:7);
        otherwise
            patientID = name(1:6);
    end
         
    ft1reg = [pathstr '/' patientID '.t1_reg.nii.gz'];
    fregt1_mat = [pathstr '/' patientID '.t1_reg.rigid_registration_result.mat'];

    opts.similarity = 'sd';
    opts.dof = 6;
    opts.intensity_norm = 'histeq';
    opts.moving_mask = fmask;
    opts.static_mask = fqsm_mask;
    
    %  register_files_affine(moving_filename, static_filename, output_filename, opts);
    if ~exist(ft1reg,'file')
        register_files_affine(fbfc,foQSM,ft1reg,opts);
    else
        disp(['Found' ft1reg]);
    end

end