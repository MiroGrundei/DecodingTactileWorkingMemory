clc;clear
addpath('C:\Users\nnu04\code\MATLAB\toolboxes\spm12')
ana_dir = 'F:\decoding';

analysis = 'WMid'; 
% analysis = 'WMcontrol';

res_dir = fullfile(ana_dir, 'group_level', sprintf('decoding_%s', analysis);
if ~exist(res_dir,'dir')
    mkdir(res_dir) 
end 

% note: running subject number include four study drop-outs and two exclusions based on
% non-performance (chance level)
subnums = [0, 2, 3, 5, 6, 7, 9, 10, 11, 12, 13, 15, 17, 18, 19, 20, 21, 23, 24];

sesnums = [1, 2];

sjs = cell(numel(subnums), numel(sesnums));
for i = 1:numel(sesnums)
    for j = 1:numel(subnums)
        sjs{j, i} = sprintf('sub-2%03d_ses-task%d',subnums(j),sesnums(i));
    end
end

dec_dir = sprintf('DEC_13Bins_%s_hm_cc', analysis);
dec_bins = [5:10];
bin_dir_prefix = 'mean_bin_';
map_name = 's5cres_accuracy_minus_chance.nii,1';

design = [ones(1,6), 2*ones(1,6); 1:6, 1:6];

connames = {'Early WM', 'Late WM'};
convecs = [1 1 1 0 0 0, 1 1 1 0 0 0, ones(1,size(sjs,1))*6/size(sjs,1);
           0 0 0 1 1 1, 0 0 0 1 1 1, ones(1,size(sjs,1))*6/size(sjs,1)];

%%

matlabbatch{1}.spm.stats.factorial_design.dir = {res_dir};

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;       
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;  
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Session';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Bin';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;

scount=0;
for s = 1:size(sjs,1)
    scount = scount+1;

    imgs = cell(size(sjs,2)*length(dec_bins),1);
    bcount = 0;
    for ses = 1:size(sjs,2)
        sub_dir = fullfile(ana_dir, sjs{s}(1:8), sjs{s, ses}, dec_dir);
        for b = dec_bins
            bcount = bcount+1;
            map = fullfile(sub_dir, sprintf('%s%d',bin_dir_prefix,b), map_name);
            imgs{bcount} = map;
        end
    end

    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(scount).scans = imgs;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(scount).conds = design; 

end

matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = [2; 3];

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''}; 
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