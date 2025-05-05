function decoding_hpc(subject_id, session_id, analysis_steps, prefix)
% Decoding HPC
% -------------------------------------------------------------------------
% subject_id: double

analysis_switch = analysis_steps;
currPrefix = prefix; 
rpPrefix = '';
ccPrefix = 'r';

cat_norm = 1;
do_coreg = 0;

% Change details here:
ntask = session_id;  
control_analysis = 0;
hm = 1;                                                                     % head motion parameters from realignment (step 4 in B0_preprocessing)
cc = 1;                                                                     % CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
% here we model the whole trial with 13 Bins
WMdelay = 13;                                
rad = 6;  
s_kernel = [5 5 5];
condnamesFIR = {'WM1', 'WM2', 'WM3', 'WM4'};

dec_type = "searchlight";                                                   % searchlight or roi
dec_task = "classification";                                                % regression (SVR) and classification (SVM) are implemented so far

% selection of analysis steps (1-7) to be performed
% -------------------------------------------------------------------------

% 1: Read log-files, extract onsets
% 2: Create FIR model
% 3: Create GLM 
% 4: Run Decoding
% 5: Average Decoding results
% 6: Normalize accuracy maps
% 7: Smooth accuracy maps

% Set paths
% -------------------------------------------------------------------------
src_dir      = '/scratch/grum90/memoslap/analysis/';
logDir 	     = '/scratch/grum90/memoslap/data/logs_task/';

addpath(genpath('/home/grum90/memoslap/decoding/')); 
addpath(genpath('/scratch/grum90/toolboxes/hMRI-toolbox-0.6.0'));
addpath('/scratch/grum90/toolboxes/tdt_3.999F/decoding_toolbox');
SPM_path = '/scratch/grum90/toolboxes/spm12/';            
addpath(SPM_path);

% Get data
% -------------------------------------------------------------------------
cd(src_dir) 
SJ = sprintf('sub-2%03d', subject_id);

% session & run identifiers
session = sprintf('%s_ses-task%d', SJ, ntask);
ses_dir = fullfile(src_dir, SJ, session);

% get runs
run_dir = fullfile(ses_dir, 'func');
rd = dir(fullfile(run_dir, 'sub*.nii'));
for r = 1:length(rd)
    runs{r} = rd(r).name;
end

% anatomy identifier
struct_dir = fullfile(ses_dir, 'anat');

fprintf('Decoding - Analysing Data (prefix=%s): \n', prefix)
fprintf('========================\n\n')
fprintf('\nSubject: %s \n', SJ)
fprintf('\nSesssion: %s \n', session)
fprintf('\nRuns: \n')
for r = 1:size(runs,2)
    fprintf('%s \n', runs{r})
end

% get the data from the json file 
json_file = fullfile(run_dir, [runs{1}(1:end-3) 'json']); %we select the first json file to extract metadata from 

fprintf('Metadata: %s\n', json_file)

TR_json = get_metadata_val(json_file,'RepetitionTime') / 1000; % repetition time in sec
slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
n_slices_json = height(slice_timing); %compute number of slices from slice timing
[~,y]= sort(slice_timing); %compute slice order 
slice_order = y';

% get the same info from nifti header  
nifti_file_metadata = fullfile(run_dir, runs{1}); 

fprintf('Nifti header: %s\n', nifti_file_metadata)

info = niftiinfo(nifti_file_metadata);
TR_nifti = info.PixelDimensions(4); 
n_slices_nifti = info.ImageSize(3);

vox_size=repmat(info.PixelDimensions(1),1,3);

% compare json and nifti header 
if round(TR_nifti, 4) ~= round(TR_json, 4)
    warning ("TR does not match between json file and nifti") 
end 
if n_slices_json ~= n_slices_nifti
    warning ("Number of slices does not match between json file and nifti") 
end 

n_slices = n_slices_json; % number of slices
TR=TR_json; % repetition time in sec.


%%
isglm = 0;

if control_analysis 
    if hm
        namepart = 'WMcontrol_hm';
    else
        namepart = 'WMcontrol_nohm';
    end
else
    if hm
        namepart = 'WMid_hm';
    else
        namepart = 'WMid_nohm';
    end   
end

if cc
    namepart = [namepart '_cc' ccPrefix ];
end  


% 2) FIR
% -------------------------------------------------------------------------
refslice = slice_order(round(length(slice_order)/2));
fir_out = sprintf('%dBins_%s', WMdelay/TR, namepart);

% 3) GLM
% -------------------------------------------------------------------------
% 1st level glm for NOT-normalized data
% must include variable 'onsets' (cell with sj x runs x conditions) including onset-times in sec
condnames = {''};  
duration = 0;                                                               % epoch duration; for single events set 0
tr = TR;
fmri_t = n_slices;                                                          % Microtime resolution; If you have performed slice-timing correction, change this parameter to match the number of slices specified there; otherwise, set default 16
fmri_t0	= refslice;
% do you want to normalize the realigned data too, to use masks in MNI-space for ROI-based decoding?
% if yes (1) nomalization will be initiated before the construction of the glm

% 4) Decoding
% -------------------------------------------------------------------------
Nbins = WMdelay/TR; 
con_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
for i = 1:Nbins
    for j = 1:max(max(con_pairs))
        labelnames{j, i} = sprintf('WM%d bin %d', j, i);
    end
end
outputDir4 = 'D4';

sbj_level_folder = sprintf('DEC_%dBins_evenOnset_%s_%dvoxRad', Nbins, namepart, rad);

if dec_type == "roi" 
	cd('');                                                                 % dir that includes all roi-masks (binary!!)
	masks = dir(['*.nii']);                                                 % or whatever you use as an identifier
else
    cd(src_dir);
    masks = dir('mask.nii');
end

% 5) Average over conditions
% -------------------------------------------------------------------------
outputDir5 = 'mean';
% (9s at TR = 1,5 equals 6 Intervals; 12s at TR = 2 equals 6 Intervals)

% 6) Normalize
% -------------------------------------------------------------------------
% T1 MUST HAVE BEEN SEGMENTED AND FUNCTIONAL DATA MUST HAVE
% BEEN REALIGNED (see preprocessing Batch B0 for these steps) 

% 7) Smooth
% -------------------------------------------------------------------------

%% Run
% -------------------------------------------------------------------------
for n = analysis_switch
    switch n	

		case 1 % Logfiles / Onsets
        % -----------------------------------------------------------------

        % get logfiles
        oruns=dir(fullfile(logDir, sprintf('%s_ses-task', SJ), sprintf('%s_tDCS_TWMD_ses-%02d_run*.tsv', SJ, ntask)));
        rmlogs =[];
        for itemp = 1:numel(oruns)
            try
                temp = tdfread([oruns(itemp).folder filesep oruns(itemp).name]);
            catch err
                if err.identifier == 'MATLAB:repmat:invalidReplications';
                    temp.Subject_Number = [];
                else
                    warning('Some Logfiles are not recognized. Skipping..')
                end
            end
            if size(temp.Subject_Number,1) < 48
                warning('Some Logfiles have short runs. Skipping..')
                rmlogs = [rmlogs, itemp];
            end
        end
        oruns(rmlogs) = [];
	
        % sanity checks
        % 4 logfiles expected
        if numel(oruns) ~= size(runs,2)
           error('Found %d full logfiles but expected 4.', numel(oruns))
        end
        % see if remaining logs are numbered 1-X
        temp = [];
        for itemp = 1:numel(oruns)
            temp = [temp, str2num(oruns(itemp).name(31:32))];
        end 

		use_n_runs = [1, 2, 3, 4]; 

	    if ~isempty(runs)
		    for r=1:numel(use_n_runs)
			    rr = use_n_runs(r);
        	    fprintf('Logfile - taking file: %s \n', [oruns(rr).folder filesep oruns(rr).name])
               
                % read logfile
		        this_log = tdfread([oruns(rr).folder filesep oruns(rr).name]);
    
                % Extract WM conditions
                % condition
                conditions = nan(1,48);
                for i = 1:48
                    target_label = sprintf('Stimulus_%d', this_log.Target_Stimulus(i));
                    nonTarget_label = sprintf('Stimulus_%d', 3-this_log.Target_Stimulus(i));
                    conditions(i) = this_log.(target_label)(i);
                    control_conditions(i) = this_log.(nonTarget_label)(i);
                end
                if sum(conditions==1)~=12 | sum(conditions==2)~=12 | sum(conditions==3)~=12 | sum(conditions==4)~=12
                    error('Wrong number of conditions!')
                end
                if sum(control_conditions==1)~=12 | sum(control_conditions==2)~=12 | sum(control_conditions==3)~=12 | sum(control_conditions==4)~=12
                    error('Wrong number of non-conditions!')
                end
    
	            use_onset = this_log.Trial_Onset./1000;
                uneven_index = find(mod(use_onset,1)==.5);
                use_onset(uneven_index) = use_onset(uneven_index)+0.5;
    
                % WM Stimuli:
                % Target Conditions
                % -------------------------------------------------------------
		        if control_analysis
    
    	            % cond 1: WM stim 1
	                onsets{r,1} = use_onset(control_conditions==1); 
	                % cond 2: WM stim 2
	                onsets{r,2} = use_onset(control_conditions==2);
	                % cond 3: WM stim 3
	                onsets{r,3} = use_onset(control_conditions==3);
	                % cond 4: WM stim 4
	                onsets{r,4} = use_onset(control_conditions==4);   
    
                else
		            % cond 1: WM stim 1
		            onsets{r,1} = use_onset(conditions==1); 
		            % cond 2: WM stim 2
		            onsets{r,2} = use_onset(conditions==2);
		            % cond 3: WM stim 3
		            onsets{r,3} = use_onset(conditions==3);
		            % cond 4: WM stim 4
		            onsets{r,4} = use_onset(conditions==4);   
                end
	        end 
	    end 

		case 2 % FIR
        % -----------------------------------------------------------------
	    fprintf('FIR dir: %s\n', fir_out)
        D2_FIR(ses_dir, onsets, TR, WMdelay, currPrefix, runs, fir_out, condnamesFIR, hm, cc, rpPrefix, ccPrefix)
            
        case 3 % GLM
        % -----------------------------------------------------------------
		
		beta_dir = ['1st_level_D0_' currPrefix];                            % folder that will contain the created job.mat file and SPM file         
        display(['Step 3, 1st level glm: ' SJ ])
	    
        D3_1st_level_glm(ses_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, runs, condnames, onsets, duration, hm, cc);

        case 4 % Decoding
		% -----------------------------------------------------------------
	    fprintf('DEC dir: %s\n', sbj_level_folder)
        if isglm
            beta_dir = ['1st_level_D0_' currPrefix];
        else
            beta_dir = [ses_dir filesep ['FIR_' fir_out]];
        end

        betas = beta_dir;
        dec_folder = [ses_dir filesep sbj_level_folder];
        if ~exist(dec_folder, 'dir')
            mkdir(dec_folder);
        end
           
        for b = 1:Nbins
            if dec_task == "classification"
                for cp = 1:size(con_pairs,1)                                % maybe adjust to allow restriction of used betas beforehand
                    label1 = labelnames(con_pairs(cp,1),b);
                    label2 = labelnames(con_pairs(cp,2),b);
                    these_labelnames = [label1, label2];
                    display(['Step 4, Decoding: ' SJ ', Conditions: ' cell2mat(label1) ', ' cell2mat(label2)])
                    beta_path    = fullfile(betas);
                    output_path  = fullfile(dec_folder, [outputDir4 '_' cell2mat(label1) '_' cell2mat(label2)]);
                    D4class_Decoding_batch(beta_path, output_path, these_labelnames, dec_type, rad, masks); 
                    clear label1 label2 these_labelnames
                end
            elseif dec_task == "regression"
                beta_path    = fullfile(betas);
                output_path  = fullfile(dec_folder, [outputDir4 '_bin_' num2str(b)]);
                these_labelnames = labelnames(:,b)';
                D4regress_Decoding_batch(beta_path, output_path, these_labelnames, dec_type, rad, masks); 
                clear these_labelnames
            end
        end
        
        case 5 % Average
        % -----------------------------------------------------------------
        dec_folder = [ses_dir filesep sbj_level_folder];

        for b = 1:Nbins 

            bin_folders = spm_select('List', dec_folder, 'dir', [' bin ' num2str(b) '_']);
            
            matlabbatch{1}.spm.util.imcalc.input = {strcat(dec_folder, filesep, bin_folders(1,:), filesep, 'res_accuracy_minus_chance.nii')
                                                    strcat(dec_folder, filesep, bin_folders(2,:), filesep, 'res_accuracy_minus_chance.nii')
                                                    strcat(dec_folder, filesep, bin_folders(3,:), filesep, 'res_accuracy_minus_chance.nii')
                                                    strcat(dec_folder, filesep, bin_folders(4,:), filesep, 'res_accuracy_minus_chance.nii')
                                                    strcat(dec_folder, filesep, bin_folders(5,:), filesep, 'res_accuracy_minus_chance.nii')
                                                    strcat(dec_folder, filesep, bin_folders(6,:), filesep, 'res_accuracy_minus_chance.nii')};

            matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2 + i3 + i4 + i5 + i6) / 6';
            matlabbatch{1}.spm.util.imcalc.output = 'res_accuracy_minus_chance.nii';

   	        if ~exist(strcat(dec_folder, filesep, outputDir5, ['_bin_' num2str(b)]), 'dir')
                mkdir(strcat(dec_folder, filesep, outputDir5, ['_bin_' num2str(b)]));
            end
                    
            matlabbatch{1}.spm.util.imcalc.outdir = {strcat(dec_folder, filesep, outputDir5, ['_bin_' num2str(b)])};

            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;        

            fprintf('Creating average image from \n')
            matlabbatch{1}.spm.util.imcalc.input

            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);

        end

        case 6 % Normalize
        % -----------------------------------------------------------------
        if isglm
            data_dir = [ses_dir filesep '1st_level_D0_' currPrefix];
            f3 = spm_select('List', data_dir, '.*con.*.nii');
            numVols = size(f3,1);
            Images{1} = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);
        else
            for b = 1:Nbins
                data_dir = [ses_dir filesep sbj_level_folder filesep 'mean_bin_' num2str(b)];
                f3 = spm_select('List', data_dir, '^res_accuracy_minus_chance.nii');
                if isempty(f3)
                    data_dir = [ses_dir filesep sbj_level_folder filesep [outputDir4 '_bin_' num2str(b)]];
                    f3 = spm_select('List', data_dir, '^res_zcorr.nii');
                end
                numVols = size(f3,1);
                Images{b} = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);
            end
        end
        
        if do_coreg
            fprintf('\nCoregistration & segmentation \n')
            D6a_coregister_est(run_dir, struct_dir, Images);
            D6b_segmentation(struct_dir, SPM_path, '^s.*\.nii');
        end

		norm_prefix = 'w';

        fprintf('\nNormalizing images \n')
        for i = 1:numel(Images)
            fprintf('%s\n', Images{i}{:})
        end

        D6c_normalization(Images, struct_dir, vox_size, norm_prefix);
        currPrefix = norm_prefix;

        case 7 % Smooth
        % -----------------------------------------------------------------

        if isglm
            data_dir = [ses_dir filesep '1st_level_D0_' currPrefix];
            D7_smoothing_run(data_dir, ['^' currPrefix 'con_0001.nii'], s_kernel);
        else
            for b = 1:Nbins
                if dec_task == "classification"
                    data_dir = [ses_dir filesep sbj_level_folder filesep 'mean_bin_' num2str(b)];
                    D7_smoothing_run(data_dir, ['^' currPrefix 'res_accuracy_minus_chance.nii'], s_kernel);
                elseif dec_task == "regression"
                    data_dir = [ses_dir filesep sbj_level_folder filesep [outputDir4 '_bin_' num2str(b)]];
                    D7_smoothing_run(data_dir, ['^' currPrefix 'res_zscore.nii'], s_kernel);
                end
            end
        end

    end
end
