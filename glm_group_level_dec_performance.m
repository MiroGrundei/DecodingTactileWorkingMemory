clc;clear
addpath('C:\Users\nnu04\code\MATLAB\toolboxes\spm12')
ana_dir = 'F:\decoding';

res_dir = fullfile(ana_dir, 'group_level', 'dec_performance_groups');

if ~exist(res_dir,'dir')
    mkdir(res_dir) 
end 

% note: running subject number include four study drop-outs and two exclusions based on
% non-performance (chance level)
% low / high performer based on behav. performance
subnums_low = [2, 5, 6, 11, 13, 17, 18, 19, 21, 24];
subnums_high = [0, 3, 7, 9, 10, 12, 15, 20, 23];

sjs = [subnums_low, subnums_high];

sesnums = [1, 2];

sjs_g1 = cell(numel(subnums_low), numel(sesnums));
for i = 1:numel(sesnums)
    for j = 1:numel(subnums_low)
        sjs_g1{j, i} = sprintf('sub-2%03d_ses-task%d',subnums_low(j),sesnums(i));
    end
end

sjs_g2 = cell(numel(subnums_high), numel(sesnums));
for i = 1:numel(sesnums)
    for j = 1:numel(subnums_high)
        sjs_g2{j, i} = sprintf('sub-2%03d_ses-task%d',subnums_high(j),sesnums(i));
    end
end

dec_dir = 'DEC_13Bins_WMid_hm_cc';
dec_bins = [5:10];
bin_dir_prefix = 'mean_bin_';
map_name = 's5cres_accuracy_minus_chance.nii,1';

ng1 = size(subnums_low,2);
ng2 = size(subnums_high,2);

design_group1 = [ones(1,6)];
design_group2 = [2*ones(1,6)];

connames = {'Low < High', 'Low > High'};
convecs = [-1 1, ones(1,ng1)*-1/ng1, ones(1,ng2)*1/ng2;
           1 -1, ones(1,ng1)*1/ng1, ones(1,ng2)*-1/ng2];

%%

matlabbatch{1}.spm.stats.factorial_design.dir = {res_dir};

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;      
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;  
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Performance';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;

scount=0;
for s = 1:size(sjs_g1,1)
    scount = scount+1;
    imgs = cell(size(sjs_g1,2)*length(dec_bins),1);
    bcount = 0;
    for ses = 1:size(sjs_g1,2)
        sub_dir = fullfile(ana_dir, sjs_g1{s}(1:8), sjs_g1{s, ses}, dec_dir);
        for b = dec_bins
            bcount = bcount+1;
            map = fullfile(sub_dir, sprintf('%s%d',bin_dir_prefix,b), map_name);
            imgs{bcount} = map;
        end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(scount).scans = imgs;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(scount).conds = design_group1; 
end

for s = 1:size(sjs_g2,1)
    scount = scount+1;
    imgs = cell(size(sjs_g2,2)*length(dec_bins),1);
    bcount = 0;
    for ses = 1:size(sjs_g2,2)
        sub_dir = fullfile(ana_dir, sjs_g2{s}(1:8), sjs_g2{s, ses}, dec_dir);
        for b = dec_bins
            bcount = bcount+1;
            map = fullfile(sub_dir, sprintf('%s%d',bin_dir_prefix,b), map_name);
            imgs{bcount} = map;
        end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(scount).scans = imgs;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(scount).conds = design_group2; 
end

matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 2;

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'E:\MeMoSLAP\masks\mask.nii'};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE GLM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Creating GLM\n')

spm_jobman('run', matlabbatch);
clear matlabbatch

%%%%%%%%%%%%%%%%%%%%%%%%%% ESTIMATE GLM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(fullfile(res_dir, 'SPM.mat'));

fprintf('Estimating GLM \n');
cd(res_dir);
SPM = spm_spm(SPM);
clear SPM;

%%%%%%%%%%%%%%%%%%%%%%%%%% Contrasts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Contrasts \n');

matlabbatch{1}.spm.stats.con.spmmat = {fullfile(res_dir, 'SPM.mat')};
% Cycle over contrast specifications
for c = 1:numel(connames)
    % Allocate t-contrast structure
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name    = connames{c};       
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = convecs(c,:);         
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';             
end 
% Delete existing contrasts (1=yes)
matlabbatch{1}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch);
clear matlabbatch
