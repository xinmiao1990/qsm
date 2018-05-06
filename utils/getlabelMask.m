
function [maskMat] = getlabelMask(flabel_trans, structName, patientType)

    %%  label ID
    %%612-617, 640-641
    %right = even, left = odd
    
    if strcmp(structName, 'CN')
        labelID = [612,613];
    elseif strcmp(structName, 'PT')
        labelID = [614,615];
    elseif strcmp(structName, 'GP')  
        labelID = [616,617];
    elseif strcmp(structName, 'CSF')
        labelID = [720,721];
    end
    
    [pathstr,name,~] = fileparts(flabel_trans);
    
    switch patientType
        case 'ACTL'
            patientID = name(1:7);
        case 'QCTL'
            patientID = name(1:7);
        otherwise
            patientID = name(1:6);
    end
    
    flabelmask = [pathstr '/' patientID '.' structName '.qsmMask.nii'];
    

    lbl = load_untouch_nii_gz(flabel_trans);

    mask = lbl;
    mask.img=zeros(size(lbl.img));
    
    for n=1:length(labelID)
        mask.img(lbl.img==labelID(n))=1;
    end

%     displayNII(mask);
    
    save_untouch_nii(mask, flabelmask);
    
    maskMat = logical(mask.img);

end