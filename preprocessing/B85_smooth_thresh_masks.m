function B85_smooth_thresh_masks(mask1, mask2, kernel, thresh_csf, thresh_wm)

% A CompCor method as described in Behzadi et al., 2007

spm('defaults', 'FMRI');

% Mask1=CSF; Mask2=WM

fileset = {[mask1 ',1']; [mask2 ',1']}

matlabbatch{1}.spm.spatial.smooth.data      = fileset;
matlabbatch{1}.spm.spatial.smooth.fwhm      = [kernel kernel kernel];
matlabbatch{1}.spm.spatial.smooth.dtype     = 0;
matlabbatch{1}.spm.spatial.smooth.im        = 0;
matlabbatch{1}.spm.spatial.smooth.prefix    = sprintf('s%d', kernel);

spm_jobman('run', matlabbatch);
clear matlabbatch

cc_prefix = [sprintf('s%d', kernel)];

% CSF MASK
[pathstr, fname, ext] = fileparts(mask1);
matlabbatch{1}.spm.util.imcalc.input = {[pathstr '/' cc_prefix fname '.nii,1']};
matlabbatch{1}.spm.util.imcalc.expression       = ['(i1>' sprintf('%.2f)', thresh_csf)];
matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
matlabbatch{1}.spm.util.imcalc.options.mask     = 1;
matlabbatch{1}.spm.util.imcalc.options.interp   = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype    = 4;
matlabbatch{1}.spm.util.imcalc.output = [pathstr '/' sprintf('tCSF%dtWM%d%s', thresh_csf*100, thresh_wm*100, cc_prefix) fname '.nii'];

spm_jobman('run', matlabbatch);
clear matlabbatch

% WM MASK
[pathstr, fname, ext] = fileparts(mask2);
matlabbatch{1}.spm.util.imcalc.input = {[pathstr '/' cc_prefix fname '.nii,1']};
matlabbatch{1}.spm.util.imcalc.expression       = ['(i1>' sprintf('%.2f)', thresh_wm)];
matlabbatch{1}.spm.util.imcalc.options.dmtx     = 0;
matlabbatch{1}.spm.util.imcalc.options.mask     = 1;
matlabbatch{1}.spm.util.imcalc.options.interp   = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype    = 4;
matlabbatch{1}.spm.util.imcalc.output = [pathstr '/' sprintf('tCSF%dtWM%d%s', thresh_csf*100, thresh_wm*100, cc_prefix) fname '.nii'];

spm_jobman('run', matlabbatch);
clear matlabbatch


end