function D6c_normalization(Images, struct_dir, vox_size, prefix)

spm('defaults','fmri');
spm_jobman('initcfg');

warning off

f1 = spm_select('List', struct_dir, '^y.*\.nii$');
numVols = size(f1,1);
parameter_file = cellstr([repmat([struct_dir filesep], numVols, 1) f1 ]);

%--------------------------------------------------------------------------
% --------------------------- Normalize  ----------------------------
%--------------------------------------------------------------------------

for i = 1:size(Images, 2)


    matlabbatch{1}.spm.spatial.normalise.write.subj.def = parameter_file;
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = Images{i};
    
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70;78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox_size;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = prefix;
    
    spm_jobman('run', matlabbatch)
    clear jobs
end