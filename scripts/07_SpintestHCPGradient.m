% require software download from: https://github.com/spin-test/spin-test 
% (by Aaron Alexander-Bloch & Siyuan Liu)
% & installation of FreeSurfer Matlab toolboxes (part of FreeSurfer installation) 

clear variables
close all

%% set up dirs
base_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram';
gradient_data_dir = [base_dir, '/data/Gradients_Margulies2016/fsaverage/nifti'];
cell_type_dir = [base_dir, '/para/data/ahba'];
output_data_dir = [base_dir,'/data/Gradients_Margulies2016'];
schaeffer_dir = [base_dir,'/data/Schaefer'];

fshome = '/Applications/freesurfer/7.2.0';
fsmatlab = '/Applications/freesurfer/7.2.0/matlab';
spintest_dir = [base_dir,'/software/spin-test-master'];
spintest_scripts_dir = [spintest_dir,'/scripts'];

path(path,fsmatlab);
path(path,spintest_dir);
path(path,spintest_scripts_dir);
path(path,[base_dir,'/software/ENIGMA-1.1.3/matlab/shared/parcellations']);
path(path,[base_dir,'/software/ENIGMA-1.1.3/matlab/scripts/useful']);
addpath(genpath('/Users/zhangxihan/Downloads/ENIGMA-1.1.3/matlab/'));


%% prepare the data (in fsaverage5, Kevin's work is in fsaverage6)
% 0. unzip the .nii.gz then convert to csv
%gradient_1_l = [gradient_data_dir,'/hcp.embed.grad_1.L.fsa6.func.nii.gz'];
%gunzip(gradient_1_l);
%gradient_1_l = [gradient_data_dir,'/hcp.embed.grad_1.L.fsa6.func.nii'];
% did this in python

% 1. path to csv files
gradient_1_l = [gradient_data_dir,'/hcp.embed.grad_1.L.fsa5.func.csv']
gradient_1_r = [gradient_data_dir,'/hcp.embed.grad_1.R.fsa5.func.csv']

gradient_2_l = [gradient_data_dir,'/hcp.embed.grad_2.L.fsa5.func.csv']
gradient_2_r = [gradient_data_dir,'/hcp.embed.grad_2.R.fsa5.func.csv']

% 2. load schaefer 400 labels for fsaverage5 vertex:
schaefer_400_fsa5 = csvread([base_dir,'/software/ENIGMA-1.1.3/matlab/shared/parcellations/schaefer_400_fsa5.csv']);
% make mask for 0
schaefer_400_fsa5_l = schaefer_400_fsa5(1:10242);
schaefer_400_fsa5_r = schaefer_400_fsa5(10243:end);
schaefer_400_fsa5_l(schaefer_400_fsa5_l==0) = 500;
schaefer_400_fsa5_r(schaefer_400_fsa5_r==0) = 500;
schaefer_400_fsa5_l(schaefer_400_fsa5_l~=500) = 0;
schaefer_400_fsa5_r(schaefer_400_fsa5_r~=500) = 0;
schaefer_400_fsa5_l(schaefer_400_fsa5_l==500) = 1;
schaefer_400_fsa5_r(schaefer_400_fsa5_r==500) = 1;
writematrix(schaefer_400_fsa5_l,[output_data_dir,'/schaefer_400_fsa5_l_mask.csv']);
writematrix(schaefer_400_fsa5_r,[output_data_dir,'/schaefer_400_fsa5_r_mask.csv']);
% path to mask of vertex labeled as 0 in schaefer 400
schaefer_400_fsa5_l = [output_data_dir,'/schaefer_400_fsa5_l_mask.csv'];
schaefer_400_fsa5_r = [output_data_dir,'/schaefer_400_fsa5_r_mask.csv'];

%% perform spintest

SpinPermuFS(fshome,fsmatlab,gradient_1_l,gradient_1_r,schaefer_400_fsa5_l,schaefer_400_fsa5_r,1000,[output_data_dir,'/rotation.hcp.embed.grad_1.fsa5.mat'])

SpinPermuFS(fshome,fsmatlab,gradient_2_l,gradient_2_r,schaefer_400_fsa5_l,schaefer_400_fsa5_r,1000,[output_data_dir,'/rotation.hcp.embed.grad_2.fsa5.mat'])

%% parcelize the vertex

leftmask=importdata(schaefer_400_fsa5_l);
rightmask=importdata(schaefer_400_fsa5_r);
mask = [leftmask; rightmask];

load([output_data_dir,'/rotation.hcp.embed.grad_1.fsa5.mat']);
CT = [bigrotl bigrotr];

gradient_1_schaefer_400 = zeros(1000,401);
for ii=1:1000
    this_CT = CT(ii,:)';
    CT_new = zeros(size(this_CT));

    % find the index of shuffled vertex that was labeled as 0 in schaefer 400
    index = find(this_CT~=100);
    if size(index,1)<sum(mask==0)
        index_pad = randsample(index,sum(mask==0)-size(index,1));
        index = [index;index_pad];
    elseif size(index,1)>sum(mask==0)
        index = index(1:sum(mask==0));
    end
    
    % fill the values to the none 0-label-vertex in original atalas
    CT_new(mask==0) = this_CT(index);
    
    % parcelize
    CT_new = CT_new';
    CT_schaefer_400 = surface_to_parcel(CT_new, 'schaefer_400_fsa5');
    gradient_1_schaefer_400(ii,:) = CT_schaefer_400;
end

writematrix(gradient_1_schaefer_400,[output_data_dir,'/gradient_1_schaefer_400_null.csv']);
gradient_1_schaefer_400 = csvread([output_data_dir,'/gradient_1_schaefer_400_null.csv']);
gradient_1_schaefer_400(:,1)=[];
bigrotl=gradient_1_schaefer_400(:,1:200);
bigrotr=gradient_1_schaefer_400(:,201:400);
save([output_data_dir,'/gradient_1_schaefer_400_null.mat'],'bigrotl','bigrotr');

load([output_data_dir,'/rotation.hcp.embed.grad_2.fsa5.mat']);
CT = [bigrotl bigrotr];

gradient_2_schaefer_400 = zeros(1000,401);
for ii=1:1000
    this_CT = CT(ii,:)';
    CT_new = zeros(size(this_CT));

    % find the index of shuffled vertex that was labeled as 0 in schaefer 400
    index = find(this_CT~=100);
    if size(index,1)<sum(mask==0)
        index_pad = randsample(index,sum(mask==0)-size(index,1));
        index = [index;index_pad];
    elseif size(index,1)>sum(mask==0)
        index = index(1:sum(mask==0));
    end
    
    % fill the values to the none 0-label-vertex in original atalas
    CT_new(mask==0) = this_CT(index);
    
    % parcelize
    CT_new = CT_new';
    CT_schaefer_400 = surface_to_parcel(CT_new, 'schaefer_400_fsa5');
    gradient_2_schaefer_400(ii,:) = CT_schaefer_400;
end

writematrix(gradient_2_schaefer_400,[output_data_dir,'/gradient_2_schaefer_400_null.csv']);
gradient_2_schaefer_400 = csvread([output_data_dir,'/gradient_2_schaefer_400_null.csv']);
gradient_2_schaefer_400(:,1)=[];
bigrotl=gradient_2_schaefer_400(:,1:200);
bigrotr=gradient_2_schaefer_400(:,201:400);
save([output_data_dir,'/gradient_2_schaefer_400_null.mat'],'bigrotl','bigrotr');
%% calculate p-value
% 1. prepare the cell-type distribution matrix
schaefer_400_LakeDFC = readtable([output_data_dir,'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3.csv'],'ReadRowNames',true);
schaefer_400_LakeVIS = readtable([output_data_dir,'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3.csv'],'ReadRowNames',true);

% list of cell types
cell_type_list_LakeDFC = schaefer_400_LakeDFC.Properties.RowNames;
cell_type_list_LakeVIS = schaefer_400_LakeVIS.Properties.RowNames;

% split the table into left and right halves:
schaefer_400_LakeDFC_l = table2array(schaefer_400_LakeDFC(:,1:200));
schaefer_400_LakeDFC_l_dir = [output_data_dir,'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3_l.csv'];
writematrix(schaefer_400_LakeDFC_l,schaefer_400_LakeDFC_l_dir);

schaefer_400_LakeDFC_r = table2array(schaefer_400_LakeDFC(:,201:400));
schaefer_400_LakeDFC_r_dir = [output_data_dir,'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3_r.csv'];
writematrix(schaefer_400_LakeDFC_r,schaefer_400_LakeDFC_r_dir);

schaefer_400_LakeVIS_l = table2array(schaefer_400_LakeVIS(:,1:200));
schaefer_400_LakeVIS_l_dir = [output_data_dir,'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3_l.csv'];
writematrix(schaefer_400_LakeVIS_l,schaefer_400_LakeVIS_l_dir);

schaefer_400_LakeVIS_r = table2array(schaefer_400_LakeVIS(:,201:400));
schaefer_400_LakeVIS_r_dir = [output_data_dir,'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3_r.csv'];
writematrix(schaefer_400_LakeVIS_r,schaefer_400_LakeVIS_r_dir);
% mask for the NaN
schaefer_400_DFC_mask_l = zeros(1,size(schaefer_400_LakeDFC_l,2));
schaefer_400_DFC_mask_l(isnan(schaefer_400_LakeDFC_l(1,:))) = 1;
schaefer_400_DFC_mask_l = logical(schaefer_400_DFC_mask_l);
save([gradient_data_dir,'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3_l_mask.mat'],'schaefer_400_DFC_mask_l');
load([gradient_data_dir,'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3_l_mask.mat']);

schaefer_400_DFC_mask_r = zeros(1,size(schaefer_400_LakeDFC_r,2));
schaefer_400_DFC_mask_r(isnan(schaefer_400_LakeDFC_r(1,:))) = 1;
schaefer_400_DFC_mask_r = logical(schaefer_400_DFC_mask_r);
save([gradient_data_dir,'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3_r_mask.mat'],'schaefer_400_DFC_mask_r');
load([gradient_data_dir,'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3_r_mask.mat']);

schaefer_400_VIS_mask_l = zeros(1,size(schaefer_400_LakeVIS_l,2));
schaefer_400_VIS_mask_l(isnan(schaefer_400_LakeVIS_l(1,:))) = 1;
schaefer_400_VIS_mask_l = logical(schaefer_400_VIS_mask_l);
save([gradient_data_dir,'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3_l_mask.mat'],'schaefer_400_VIS_mask_l');
load([gradient_data_dir,'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3_l_mask.mat']);

schaefer_400_VIS_mask_r = zeros(1,size(schaefer_400_LakeVIS_r,2));
schaefer_400_VIS_mask_r(isnan(schaefer_400_LakeVIS_r(1,:))) = 1;
schaefer_400_VIS_mask_r = logical(schaefer_400_VIS_mask_r);
save([gradient_data_dir,'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3_r_mask.mat'],'schaefer_400_VIS_mask_r');
load([gradient_data_dir,'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3_r_mask.mat']);

% 2. prepare the real gradient value
gradient_1_l = importdata([gradient_data_dir,'/hcp.embed.grad_1.L.fsa5.func.csv']);
gradient_1_r = importdata([gradient_data_dir,'/hcp.embed.grad_1.R.fsa5.func.csv']);
gradient_1 = [gradient_1_l;gradient_1_r];
gradient_1 = gradient_1';
gradient_1_schaefer_400_fsa5 = surface_to_parcel(gradient_1, 'schaefer_400_fsa5');
gradient_1_schaefer_400_fsa5(1)=[];
writematrix(gradient_1_schaefer_400_fsa5,[output_data_dir,'/gradient_1_schaefer_400_fsa5.csv']);

gradient_1_schaefer_400_fsa5_l = gradient_1_schaefer_400_fsa5(1:200);
gradient_1_schaefer_400_fsa5_l_dir = [output_data_dir,'/gradient_1_schaefer_400_fsa5_l.csv'];
writematrix(gradient_1_schaefer_400_fsa5_l,gradient_1_schaefer_400_fsa5_l_dir);

gradient_1_schaefer_400_fsa5_r = gradient_1_schaefer_400_fsa5(201:400);
gradient_1_schaefer_400_fsa5_r_dir = [output_data_dir,'/gradient_1_schaefer_400_fsa5_r.csv'];
writematrix(gradient_1_schaefer_400_fsa5_r,gradient_1_schaefer_400_fsa5_r_dir);

gradient_2_l = importdata([gradient_data_dir,'/hcp.embed.grad_2.L.fsa5.func.csv']);
gradient_2_r = importdata([gradient_data_dir,'/hcp.embed.grad_2.R.fsa5.func.csv']);
gradient_2 = [gradient_2_l;gradient_2_r];
gradient_2 = gradient_2';
gradient_2_schaefer_400_fsa5 = surface_to_parcel(gradient_2, 'schaefer_400_fsa5');
gradient_2_schaefer_400_fsa5(1)=[];
writematrix(gradient_2_schaefer_400_fsa5,[output_data_dir,'/gradient_2_schaefer_400_fsa5.csv']);

gradient_2_schaefer_400_fsa5_l = gradient_2_schaefer_400_fsa5(1:200);
gradient_2_schaefer_400_fsa5_l_dir = [output_data_dir,'/gradient_2_schaefer_400_fsa5_l.csv'];
writematrix(gradient_2_schaefer_400_fsa5_l,gradient_2_schaefer_400_fsa5_l_dir);

gradient_2_schaefer_400_fsa5_r = gradient_2_schaefer_400_fsa5(201:400);
gradient_2_schaefer_400_fsa5_r_dir = [output_data_dir,'/gradient_2_schaefer_400_fsa5_r.csv'];
writematrix(gradient_2_schaefer_400_fsa5_r,gradient_2_schaefer_400_fsa5_r_dir);

% 3. the 1000 permutated gradient in schaefer 400
gradient_1_schaefer_400_null_dir = [output_data_dir,'/gradient_1_schaefer_400_null.mat'];
gradient_2_schaefer_400_null_dir = [output_data_dir,'/gradient_2_schaefer_400_null.mat'];

% 4. calculate pvalue
% gradient 1
for ii=1:18
    [rho_gradient_1_schaefer_400{ii}, pval_gradient_1_schaefer_400{ii}]=pvalvsNull(gradient_1_schaefer_400_fsa5_l_dir,...
        gradient_1_schaefer_400_fsa5_r_dir,...
        schaefer_400_LakeDFC_l_dir,...
        schaefer_400_LakeDFC_r_dir,...
        ii,1000,...
        gradient_1_schaefer_400_null_dir, ...
        schaefer_400_DFC_mask_l, schaefer_400_DFC_mask_r);
end

% gradient 2
for ii=1:18
    [rho_gradient_2_schaefer_400{ii}, pval_gradient_2_schaefer_400{ii}]=pvalvsNull(gradient_2_schaefer_400_fsa5_l_dir,...
        gradient_2_schaefer_400_fsa5_r_dir,...
        schaefer_400_LakeDFC_l_dir,...
        schaefer_400_LakeDFC_r_dir,...
        ii,1000,...
        gradient_2_schaefer_400_null_dir, ...
        schaefer_400_DFC_mask_l, schaefer_400_DFC_mask_r);
end

pvalues_G1 = cell2table(pval_gradient_1_schaefer_400, "VariableNames", cell_type_list_LakeDFC,'RowNames',{'Gradeint1_pvalue'});
pvalues_G2 = cell2table(pval_gradient_2_schaefer_400, "VariableNames", cell_type_list_LakeDFC,'RowNames',{'Gradeint2_pvalue'});
rhos_G1 = cell2table(rho_gradient_1_schaefer_400, "VariableNames", cell_type_list_LakeDFC,'RowNames',{'Gradeint1_r'});
rhos_G2 = cell2table(rho_gradient_2_schaefer_400, "VariableNames", cell_type_list_LakeDFC,'RowNames',{'Gradeint2_r'});
rho_pvalues_G1_G2 = [rhos_G1; pvalues_G1; rhos_G2; pvalues_G2];
writetable(rho_pvalues_G1_G2,[output_data_dir,'/CellType_Gradients_corr_HCP_LakeDFC_NormZscore0.3.csv'],'WriteRowNames',true);

% VIS
% gradient 1
for ii=1:18
    [rho_gradient_1_schaefer_400{ii}, pval_gradient_1_schaefer_400{ii}]=pvalvsNull(gradient_1_schaefer_400_fsa5_l_dir,...
        gradient_1_schaefer_400_fsa5_r_dir,...
        schaefer_400_LakeVIS_l_dir,...
        schaefer_400_LakeVIS_r_dir,...
        ii,1000,...
        gradient_1_schaefer_400_null_dir, ...
        schaefer_400_VIS_mask_l, schaefer_400_VIS_mask_r);
end

% gradient 2
for ii=1:18
    [rho_gradient_2_schaefer_400{ii}, pval_gradient_2_schaefer_400{ii}]=pvalvsNull(gradient_2_schaefer_400_fsa5_l_dir,...
        gradient_2_schaefer_400_fsa5_r_dir,...
        schaefer_400_LakeVIS_l_dir,...
        schaefer_400_LakeVIS_r_dir,...
        ii,1000,...
        gradient_2_schaefer_400_null_dir, ...
        schaefer_400_VIS_mask_l, schaefer_400_VIS_mask_r);
end

pvalues_G1 = cell2table(pval_gradient_1_schaefer_400, "VariableNames", cell_type_list_LakeVIS,'RowNames',{'Gradeint1_pvalue'});
pvalues_G2 = cell2table(pval_gradient_2_schaefer_400, "VariableNames", cell_type_list_LakeVIS,'RowNames',{'Gradeint2_pvalue'});
rhos_G1 = cell2table(rho_gradient_1_schaefer_400, "VariableNames", cell_type_list_LakeVIS,'RowNames',{'Gradeint1_r'});
rhos_G2 = cell2table(rho_gradient_2_schaefer_400, "VariableNames", cell_type_list_LakeVIS,'RowNames',{'Gradeint2_r'});
rho_pvalues_G1_G2 = [rhos_G1; pvalues_G1; rhos_G2; pvalues_G2];
writetable(rho_pvalues_G1_G2,[output_data_dir,'/CellType_Gradients_corr_HCP_LakeVIS_NormZscore0.3.csv'],'WriteRowNames',true);

%%
read_annotation(fullfile('/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Parcellations/FreeSurfer5.3/fsaverage6/label/lh.Schaefer2018_400Parcels_17Networks_order.annot'));
read_annotation(fullfile(fshome,'/subjects/fsaverage5/label/lh.aparc.a2009s.annot'));