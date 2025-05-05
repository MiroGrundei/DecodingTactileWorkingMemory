function D2_FIR(data_dir, onsets, TR, WMdelay, currPrefix, runs, fir_out, condnames, hm, cc, rpPrefix, ccPrefix)

% This script implements a FIR model as commonly used in Time-resolved
% Decoding.
% Each time-bin (bin width = TR) of a delay phase is modelled individually

% SPM defaults
spm('defaults','fmri');
global defaults;
global UFp;
spm_jobman('initcfg');
% OUTPUT Directory (as subdirectory of the SJ directory)
tgt_dir = fullfile(data_dir, ['FIR_' fir_out]);
if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
end

%******************************************************************
%               Setting FIR parameter
% *****************************************************************
% Output directory
jobs{1}.stats{1}.fmri_spec.dir = cellstr(tgt_dir);
% timing parameters
jobs{1}.stats{1}.fmri_spec.timing.units     = 'secs';
jobs{1}.stats{1}.fmri_spec.timing.RT        = TR;
jobs{1}.stats{1}.fmri_spec.timing.fmri_t    = 16;
jobs{1}.stats{1}.fmri_spec.timing.fmri_t0   = 1;
% FIR-SPECIFICATION
jobs{1}.stats{1}.fmri_spec.bases.fir.length = WMdelay;      % seconds
jobs{1}.stats{1}.fmri_spec.bases.fir.order  = WMdelay/TR;   % time-bins
% Other Specifications
jobs{1}.stats{1}.fmri_spec.fact              = struct('name', {}, 'levels', {});
jobs{1}.stats{1}.fmri_spec.volt              = 1;
jobs{1}.stats{1}.fmri_spec.global            = 'None';
jobs{1}.stats{1}.fmri_spec.mask              = {''};
jobs{1}.stats{1}.fmri_spec.cvi               = 'None'; 
jobs{1}.stats{1}.fmri_spec.mthresh 	         = 0.1;


% *****************************************************************
%            Specify the Desing/Conditions/Onsets
% *****************************************************************
% get all runIDs from file names in 'func' folder
nifti_dir = fullfile(data_dir, 'func');
runIDs = 1:4;

% Loop over runs, as the logfiles (with onset data) are saved
% separately per session
for s = 1:size(runs,2)
    runID = runIDs(s);
    % high-pass cut-off 
    jobs{1}.stats{1}.fmri_spec.sess(s).hpf     = 192; 
    % Allocation of Data (EPIs/Images) for the current Session
    % TODO get rid of hardcoded file names here
    filename = fullfile(nifti_dir, [currPrefix runs{runID}]);

    fprintf('FIR: %s \n', filename)

    % use 'expand' to read 4d nifti file
    f = spm_select('expand', filename); 
    jobs{1}.stats{1}.fmri_spec.sess(s).scans   = cellstr(f);

    % Motion
    cov=1;
    multi_reg={''};
    if hm                
        k = strfind(runs{runID}, '.nii')
        f1=spm_select('List', nifti_dir, ['^rp_' rpPrefix runs{runID}(1:k-1) '.txt']);            
        if ~isempty(f1)
            multi_reg(cov,1)={[nifti_dir filesep f1]};
            cov=cov+1;
        else
	    error('No Head Motion params found: %s',nifti_dir )
        end
    end
    % CompCorr
    if cc                
        k = strfind(runs{runID}, '.nii')
        f1=spm_select('List', nifti_dir, ['^' ccPrefix runs{runID}(1:k-1) '_CompCorPCs.txt']);            
        if ~isempty(f1)
            multi_reg(cov,1)={[nifti_dir filesep f1]};
            cov=cov+1;
        else
	    error('No Comp Corr params found: %s',nifti_dir )
        end
    end

    fprintf('FIR - Addtional regressor files: \n')

    for m = 1:numel(multi_reg)
        fprintf('%s \n', multi_reg{m})
    end
    jobs{1}.stats{1}.fmri_spec.sess(s).multi_reg = multi_reg; % multiple regressors

    for c = 1:numel(condnames)
        onsets_temp = onsets{s,c};
        jobs{1}.stats{1}.fmri_spec.sess(s).cond(c).name     = condnames{c};
        jobs{1}.stats{1}.fmri_spec.sess(s).cond(c).onset    = onsets_temp;
        jobs{1}.stats{1}.fmri_spec.sess(s).cond(c).duration = 0;
        jobs{1}.stats{1}.fmri_spec.sess(s).cond(c).tmod     = 0;
        jobs{1}.stats{1}.fmri_spec.sess(s).cond(c).pmod     = struct('name', {}, 'param', {}, 'poly', {});
    end

end

% Create model
fprintf(['Creating GLM\n'])
spm_jobman('run', jobs);
clear jobs

%  Model Estimation
load(fullfile(tgt_dir, 'SPM.mat'));
fprintf(['Estimating GLM \n']);
cd(tgt_dir);
SPM = spm_spm(SPM);
clear SPM;

end