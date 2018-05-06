function transSegmentLabels_to_QSMcoord(flabels, fbfc, foQSM, fregt1_mat, flabel_trans)
    
%     [pathstr,name,~] = fileparts(foQSM);  
%     patientID = name(1:6);
%     flabel_trans = [pathstr '/' patientID '.qsm_coord.sverg.label.nii.gz'];
    %%
    data_file = flabels;
    data_coord = 'm';
    output_file = flabel_trans;
    reg_moving_file = fbfc;
    reg_static_file = foQSM;
    reg_mat_file = fregt1_mat;
    method = 'nearest';
    
    %% transform_data_affine(data_file, data_coord, output_file,
    %reg_moving_file, reg_static_file, reg_mat_file, method);
    if ~exist(flabel_trans,'file')
        transform_data_affine(data_file, data_coord, output_file, reg_moving_file, reg_static_file, reg_mat_file, method);
    else
        disp(['Found:' flabel_trans]);
    end

end