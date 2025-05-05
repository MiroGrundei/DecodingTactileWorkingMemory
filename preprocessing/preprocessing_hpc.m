function preprocessing_hpc(subject_id, session_id, analysis_steps, prefix)
% Preprocessing HPC
% -------------------------------------------------------------------------
% subject_id: double

% Change details here:
ntask = session_id;  
analysis_switch = analysis_steps; % [4, 5, 1, 6, 9]; 
start_prefix=prefix; 
%''; % if totally raw data, then keep empty, otherwise add prefix, e.g. 's8wra'
corrPrefix = ''; 
%prefix; %''; % if you perform differnt kinds of preprocessing, there will be 
% multiple 'mean...nii' files -> check which one you want and if theres anything 
% in between the 'mean' and 'sub-00...' then put that there (probably nothing or an a?)

compcorr_sj_space = 1;
do_cc_reslice = 1;
do_cc_smooth_thresh = 0;
cc_kernel = 4;
cc_thresh_wm = 0;
cc_thresh_csf = 0;

% anatomical masks (for comp corr without individual masks)
% mni_wm_mask=['C:\Users\...\wm_mask_eroded.nii']; %white matter mask file
% mni_csf_mask=['C:\Users\...\csf_mask_eroded.nii']; %csf mask file
% full_brain_mask=['C:\Users\...\full_brain_mask.nii']; %full brain mask file

% Set paths
% -------------------------------------------------------------------------
src_dir      = '/scratch/grum90/memoslap/analysis/';
addpath(genpath('/home/grum90/memoslap/preproc/')); 
addpath(genpath('/scratch/grum90/toolboxes/hMRI-toolbox-0.6.0'));
addpath(genpath('/scratch/grum90/toolboxes/DPABI-master'));
SPM_path = '/scratch/grum90/toolboxes/spm12/';          
addpath(SPM_path);

% Get data
% -------------------------------------------------------------------------
cd(src_dir) 
SJ = sprintf('sub-2%03d', subject_id);

% session & run identifiers
session = sprintf('%s_ses-task%d', SJ, ntask);
ses_dir = fullfile(src_dir, SJ, session);

% unzip
zip_files = dir(fullfile(ses_dir, '**', ['sub-', '*.gz']));
if ~isempty(zip_files) 
    fprintf('Unzipping:\n')
	for z = 1:size(zip_files, 1)
        fprintf('%s\n', [zip_files(z).folder filesep zip_files(z).name])
		gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
	end
end

% get runs
run_dir = fullfile(ses_dir, 'func');
rd = dir(fullfile(run_dir, 'sub*.nii'));
for r = 1:length(rd)
    runs{r} = rd(r).name;
end

fprintf('Analysing Data\n')
fprintf('Subject: %s \n', SJ)
fprintf('Sesssion: %s \n', session)
fprintf('Runs: \n')
for r = 1:size(runs,2)
    fprintf('%s \n', runs{r})
end

% anatomy identifier
ana = ['anat'];
struct_dir = fullfile(ses_dir, ana);

nifti_files = dir(fullfile(run_dir, ['sub-', '*bold*.nii'])); %look for all functional nifti files
anat_files = dir(fullfile(struct_dir, ['sub-', '*T1w.nii'])); %look for all anat nifti files

fprintf('Anatomy files: \n')
for i = 1:size(anat_files,1)
    fprintf('%s \n', anat_files(i).name)
end

% get the data from the json file 
json_files = (dir(fullfile(src_dir, '**', ['task', '*json']))); %extract all json files, althoguh they should have the same info 
%to check if the first command returns an empty structure. If yes, it means the json files
%have a different naming, starting with subject 
if isequal(size(json_files), [0, 1]) 
    json_files = (dir(fullfile(src_dir, '**', ['sub-', '*bold*.json'])));
    % json_files = (dir(fullfile(src_dir, '**', ['sub-', '*.json'])));
end

json_file = [json_files(1).folder, filesep, json_files(1).name]; %we select the first json file to extract metadata from 
TR_json = get_metadata_val(json_file,'RepetitionTime') / 1000; % repetition time in sec
slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
n_slices_json = height(slice_timing); %compute number of slices from slice timing
[~,y]= sort(slice_timing); %compute slice order 
slice_order = y';

% get the same info from nifti header  
nifti_file_metadata = [nifti_files(1).folder, filesep, nifti_files(1).name] ;
info = niftiinfo(nifti_file_metadata);
TR_nifti = info.PixelDimensions(4); 
n_slices_nifti = info.ImageSize(3);

% compare json and nifti header 
if round(TR_nifti, 4) ~= round(TR_json, 4)
    warning ("TR does not match between json file and nifti") 
end 
if n_slices_json ~= n_slices_nifti
    warning ("Number of slices does not match between json file and nifti") 
end 

TR = TR_json;
n_slices = n_slices_json; 

% Additional input
% -------------------------------------------------------------------------
% ------ NOTE: spike removal (e.g. "artrepair") should be performed as first step
% 1)  Segmentation
% ------ Create nuisance masks on your own or take the provided ones
% 2) --> remove first x scans                       --> prefix: x(number of cut volumes)
x=0;
% 3) --> slice time correction                      --> prefix: a
%  for interleaved slice order: do slice time correction, then realignment
%  otherwise do first realignment, then slice time correction (in analysis_switch 4 before 3)
% n_slices = 37; % number of slices
% slice_order=[1:n_slices];
refslice = slice_order(round(length(slice_order)/2)); % reference slice
% TR=2; % repetition time in sec.
%# step 4  Realignment                                --> prefix: r
% realign over all runs
%# step 5  Coregister (estimate) mean-epi 2 anatomy (DEFAULT)
%# step 5b  Coregister (estimate & resclice) mean-epi 2 anatomy --> prefix c
Co_er = 0; %default: 0, if 1, then estimate & reslice
%# step 6  Normalization                              --> prefix: w
vox_size=[2 2 2]; % preferred voxel size after Normailzation (in mm)
%vox_size=repmat(info.PixelDimensions(1),1,3); % Voxel size from JSON-file (not changed)
%# step 7  Scrubbing: calculate, interpolate outliers --> prefix: m(scrub_thresh)
scrub_thresh=0.4; % threshhold FD for scrubbing
%# step 8 Calculate WM and CSF Nuisance Signal
numComp = 5; % number of principle components
%# step 9 Smoothing                                   --> prefix: s
kernel_size=[8 8 8]; %FWHM kernel size
%# step 10 Detrending                                 --> prefix: d

%% Analyses
% -------------------------------------------------------------------------
currPrefix=start_prefix;

for n = analysis_switch
    
    switch n
        
        case 1 % Segmentation
        % -----------------------------------------------------------------
        warning off
        fprintf('\n\nStep 1, segmentation: %s', session)
        B1_segmentation(struct_dir, SJ, SPM_path, '^s.*\.nii');

        case 2 % Delete first X scans
        % -----------------------------------------------------------------
        fprintf('\n\n')
        if x>0
            for r = 1:size(runs, 2)
                fprintf('Step 2, delete first %s volumes %s %s', num2str(x), SJ, runs{r})
                B2_delete_scans(run_dir, ['^' currPrefix runs{r}], x);
            end
            currPrefix = ['x' num2str(x) currPrefix];
        end   

        case 3 % Slice time correction
        % -----------------------------------------------------------------
        fprintf('\n\n')
        for r = 1:size(runs, 2)
            fprintf('Step 3, slice time correction: %s, %s', session, runs{r})
            B3_slice_time_correction(SJ, runs{r}, run_dir, ['^' currPrefix runs{r}], n_slices, slice_order, refslice, TR);
        end
        currPrefix = ['a' currPrefix];   

        case 4 % Realignment
        % -----------------------------------------------------------------
        fprintf('\n\n')
        fprintf('Step 4, realignment: %s\n', session)
        for r = 1:size(runs, 2)
            run_files{r} = spm_select('List',run_dir,['^' currPrefix runs{r}]);
            fprintf('%s\n', run_files{r})
        end
        B4_Realignment_all_runs(run_dir, run_files);
        currPrefix = ['r' currPrefix]; %%%%%%%%%%%%%% fix
            
        case 5 % Coregister (estimate)
        % -----------------------------------------------------------------
        fprintf('\n\n')
		if Co_er ~= 1
            fprintf('Step 5, coregistration (estimate): %s\n', session)
        	B5_coregister_est(run_dir, struct_dir, '^s.*\.nii', runs, corrPrefix);
        else
            fprintf('Step 5, coregistration (estimate & reslice): %s\n', session)
			B5b_coregister_est_re(currPrefix, run_dir, struct_dir, '^s.*\.nii', runs);
		end

        case 6 % Normalization
        % -----------------------------------------------------------------
        fprintf('\n\n')
        fprintf('Step 6, normalization: %s\n', session)
        B6_normalization_run(run_dir, struct_dir, runs, vox_size, currPrefix);

        currPrefix=['w' currPrefix];
            
        case 7 % Scrubbing: calculate outliers
        % -----------------------------------------------------------------
        fprintf('\n\n')
        scrub_prefix=['m' num2str(scrub_thresh)];
        for r = 1:size(runs, 2)
            fprintf('Step 7, scrubbing: %s, %s', session, runs{r})
            %estimate and save motion statistics
            n=1;
            k = strfind(runs{r}, 'd.nii'); %% use prefix without file extension
            f=spm_select('List', run_dir, ['^rp_' currPrefix(n:end) runs{r}(1:k) '.txt']);
            while isempty(f)
                n=n+1;
                f=spm_select('List', run_dir, ['^rp_' currPrefix(n:end) runs{r}(1:k) '.txt']);
            end
            cfg.motionparam=[run_dir filesep f];
            cfg.prepro_suite = 'spm';
            [fwd,rms]=bramila_framewiseDisplacement(cfg);
            outliers=fwd>scrub_thresh;
            percent_out=(sum(outliers)/length(outliers))*100;
            disp(['outliers for ' SJ ', ' runs{r} ': ' num2str(percent_out) '%']);
            save([run_dir filesep scrub_prefix currPrefix runs{r}(1:k) '_FWDstat.mat'],'fwd','rms','outliers','percent_out','scrub_thresh','cfg')
            %srub outliers by replacing them with average of nearest neighbors
            B7_scrub_data(run_dir, ['^' currPrefix runs{r}], outliers,  scrub_prefix);
            all_percent_out(r)=percent_out;
            all_rp{r}=load(cfg.motionparam);            
        end

        currPrefix=[scrub_prefix currPrefix];
        save([src_dir filesep 'all_MOTIONstat_' currPrefix '.mat'],'SJ','runs','scrub_thresh','all_percent_out','all_rp')
            
        case 8 % CompCorr
        % -----------------------------------------------------------------
	    fprintf('\nWM / CSF Component Correction\n')
        struct_dir = fullfile(ses_dir, ana);
	    if compcorr_sj_space
            wm_mask = [struct_dir filesep 'c2' session '_T1w.nii'];
            csf_mask = [struct_dir filesep 'c3' session '_T1w.nii'];

            f2 = spm_select('List', run_dir, ['^mean' corrPrefix runs{1}]);
            numVols = size(f2,1);
            mean_img   = cellstr([repmat([run_dir filesep], numVols, 1) f2 repmat(',1', numVols, 1)]);

            B85_reslice_masks_to_functional(mean_img, wm_mask, csf_mask);

		    cc_prefix = '';

	        if do_cc_smooth_thresh
                wm_mask = [struct_dir filesep 'rc2' session '_T1w.nii'];
                csf_mask = [struct_dir filesep 'rc3' session '_T1w.nii'];
    
		        B85_smooth_thresh_masks(csf_mask, wm_mask, cc_kernel, cc_thresh_csf, cc_thresh_wm);
		        cc_prefix = sprintf('tCSF%dtWM%ds%d', cc_thresh_csf*100, cc_thresh_wm*100, cc_kernel);
	        end 
    
	        wm_mask = [struct_dir filesep cc_prefix 'rc2' session '_T1w.nii'];
            csf_mask = [struct_dir filesep cc_prefix 'rc3' session '_T1w.nii'];

        else
            wm_mask = mni_wm_mask;
            csf_mask = mni_csf_mask;
        end

        fprintf('\n\n')
        for r = 1:size(runs, 2)
            fprintf('Step 8, CompCorr: %s, %s', session, runs{r})
            B8_compcorr_run(run_dir, SJ, ['^' currPrefix runs{r}], numComp, wm_mask, csf_mask, TR, cc_prefix);
        end
   
        case 9 % Smoothing
        % -----------------------------------------------------------------  
        fprintf('\n\n')
        for r = 1:size(runs, 2)
            fprintf('Step 9, smoothing: %s, %s', session, runs{r})
            B9_smoothing_run(run_dir, SJ, ['^' currPrefix runs{r}],kernel_size);
        end
        currPrefix=['s' num2str(unique(kernel_size)) currPrefix];
            
        case 10 % Detrending
        % -----------------------------------------------------------------   
        fprintf('\n\n')
        for r = 1:size(runs, 2)
            fprintf('Step 10, Detrending: %s, %s', session, runs{r})
            f = spm_select('List',run_dir, ['^' currPrefix runs{r}]);
            V=spm_vol([run_dir filesep f(1,:)]);
            files={};
            for i=1:(size(V,1))
                files{i} = [run_dir filesep strtrim(f(1,:)) ',' int2str(i)];
            end
            fileset{r}=char(files);
        end

        if exist('fileset')
            B10_detrending_lmgs(fileset);
            clear fileset files
        end

        otherwise
        % -----------------------------------------------------------------   

        fprintf('######################################################################################')
        fprintf('############################## Case %d does not exsist ##############################', n)
        fprintf('######################################################################################')

    end
end
