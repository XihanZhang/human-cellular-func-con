
clear variables
close all

%% set up dirs
base_dir = '/Volumes/GoogleDrive/My Drive/Gradient_Shift_Cellular_Basis/Milgram';
data_dir = [base_dir, '/data/Gradients_Margulies2016'];
out_dir = [base_dir, '/data/permcca'];
figure_dir = [base_dir,'/figures/permcca'];

path(path,[base_dir,'/software/PermCCA-master']);

%% load the data
% 1. gradients
gradient_1 = csvread([data_dir,'/gradient_1_schaefer_400_fsa5.csv']);
gradient_1 = gradient_1';

gradient_2 = csvread([data_dir,'/gradient_2_schaefer_400_fsa5.csv']);
gradient_2 = gradient_2';

% 2. permutated gradients
load([data_dir,'/gradient_1_schaefer_400_null.mat']);
gradient_1_null = [bigrotl bigrotr];
gradient_1_null = gradient_1_null';

load([data_dir,'/gradient_2_schaefer_400_null.mat']);
gradient_2_null = [bigrotl bigrotr];
gradient_2_null = gradient_2_null';
clear bigrotl bigrotr

% 3. cell type abundance
LakeDFC = readtable([data_dir,'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3.csv'],'ReadRowNames',true);
LakeDFC_mat = table2array(LakeDFC)';
LakeDFC_cell_type_list = LakeDFC.Properties.RowNames;
clear LakeDFC

LakeVIS = readtable([data_dir,'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3.csv'],'ReadRowNames',true);
LakeVIS_mat = table2array(LakeVIS)';
LakeVIS_cell_type_list = LakeVIS.Properties.RowNames;
clear LakeVIS

% Removing In2 and Ex2, leaving cells common between Lake_DFC and
LakeDFC_mat(:,[2]) = []; % Ex2
LakeDFC_cell_type_list(2) = [];
LakeVIS_mat(:,[8]) = []; % In2
LakeVIS_cell_type_list(8) = [];

cell_num = size(LakeDFC_mat,2);

for i=1:cell_num
    LakeDFC_cell_type_list{i}=strrep(LakeDFC_cell_type_list{i},'DFC_','');
    LakeVIS_cell_type_list{i}=strrep(LakeVIS_cell_type_list{i},'VIS_','');
end

% 4. drop NaN parcels
nan_rows = find(isnan(LakeDFC_mat(:,1)));
LakeDFC_mat(nan_rows,:) = [];
LakeVIS_mat(nan_rows,:) = [];
gradient_1(nan_rows,:) = [];
gradient_1_null(nan_rows,:) = [];
gradient_2(nan_rows,:) = [];
gradient_2_null(nan_rows,:) = [];

%% Gradually remove significant univariant cell-types from PermCCA
nP = 1000; % permutation number
%% Gradient 1 Lake_DFC
% transmodal-distributed cell types: In1, Ex1
% unimodal-distributed cell types: End, PVALB

% 0. Remove In1, Ex1, and PVALB
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[1,7,10]) = [];
% permcca
[p_G1_DFC_In1Ex1PVALB,r_G1_DFC_In1Ex1PVALB,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 1. Remove Ex1 (Index=1)
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[1]) = [];
% permcca
[p_G1_DFC_Ex1,r_G1_DFC_Ex1,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 2. Remove In1 (Index=7)
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[7]) = [];
% permcca
[p_G1_DFC_In1,r_G1_DFC_In1,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 3. Remove End (Index=12)
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[12]) = [];
% permcca
[p_G1_DFC_End,r_G1_DFC_End,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 4. Remove PVALB (Index=10)
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[10]) = [];
% permcca
[p_G1_DFC_PVALB,r_G1_DFC_PVALB,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 5. Remove Ex1 and End
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[1,12]) = [];
% permcca
[p_G1_DFC_Ex1End,r_G1_DFC_Ex1End,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 6. Remove In1 and PVALB
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[7,10]) = [];
% permcca
[p_G1_DFC_In1PVALB,r_G1_DFC_In1PVALB,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 7. Remove Ex1 and PVALB
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[1,10]) = [];
% permcca
[p_G1_DFC_Ex1PVALB,r_G1_DFC_Ex1PVALB,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 8. Remove In1 and Ex1
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[1,7]) = [];
% permcca
[p_G1_DFC_In1Ex1,r_G1_DFC_In1Ex1,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 9. Remove End and PVALB
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[12,10]) = [];
% permcca
[p_G1_DFC_EndPVALB,r_G1_DFC_EndPVALB,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 10. Remove In1, Ex1, End, and PVALB
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[1,7,12,10]) = [];
% permcca
[p_G1_DFC_In1Ex1EndPVALB,r_G1_DFC_In1Ex1EndPVALB,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_1_null);

% Output table for p-value and r for CCA
cca_r_p = [r_G1_DFC_In1Ex1PVALB,p_G1_DFC_In1Ex1PVALB;...
    r_G1_DFC_Ex1,p_G1_DFC_Ex1;...
    r_G1_DFC_In1,p_G1_DFC_In1;...
    r_G1_DFC_End,p_G1_DFC_End;...
    r_G1_DFC_PVALB,p_G1_DFC_PVALB;...
    r_G1_DFC_Ex1End,p_G1_DFC_Ex1End;...
    r_G1_DFC_In1PVALB,p_G1_DFC_In1PVALB;...
    r_G1_DFC_Ex1PVALB,p_G1_DFC_Ex1PVALB;...
    r_G1_DFC_In1Ex1,p_G1_DFC_In1Ex1;...
    r_G1_DFC_EndPVALB,p_G1_DFC_EndPVALB;...
    r_G1_DFC_In1Ex1EndPVALB,p_G1_DFC_In1Ex1EndPVALB];

RowNames = {'In1Ex1PVALB','Ex1','In1','End','PVALB','Ex1End','In1PVALB',...
    'Ex1PVALB','In1Ex1','EndPVALB','In1Ex1EndPVALB'};
ColNames = {'r','p'}
cca_p_r_table = array2table(cca_r_p,'RowNames',RowNames,'VariableNames',ColNames);
writetable(cca_p_r_table,[out_dir,'/ccaGradualRemove_G1_DFC_r_p.csv'],'WriteRowNames',true);

%% Gradient 1 Lake_VIS
% transmodal-distributed cell types: In1, Ex1, Ex5
% unimodal-distributed cell types: PVALB

% 0. Remove In1,Ex1,PVALB
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[1,7,10]) = [];
% permcca
[p_G1_In1Ex1PVALB,r_G1_In1Ex1PVALB,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 1. Remove Ex1 (Index=1)
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[1]) = [];
% permcca
[p_G1_Ex1,r_G1_Ex1,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 2. Remove In1 (Index=7)
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[7]) = [];
% permcca
[p_G1_In1,r_G1_In1,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 3. Remove Ex5 (Index=4)
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[4]) = [];
% permcca
[p_G1_Ex5,r_G1_Ex5,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 4. Remove PVALB (Index=10)
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[10]) = [];
% permcca
[p_G1_PVALB,r_G1_PVALB,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 5. Remove Ex5, PVALB
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[4,10]) = [];
% permcca
[p_G1_Ex5PVALB,r_G1_Ex5PVALB,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 6. Remove Ex1, PVALB
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[1,10]) = [];
% permcca
[p_G1_Ex1PVALB,r_G1_Ex1PVALB,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 7. Remove In1, PVALB
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[7,10]) = [];
% permcca
[p_G1_In1PVALB,r_G1_In1PVALB,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 8. Remove Ex5, Ex1, In1
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[1,4,7]) = [];
% permcca
[p_G1_Ex5Ex1In1,r_G1_Ex5Ex1In1,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

% 9. Remove Ex5, Ex1, In1, PVALB
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[1,4,7,10,12]) = [];
% permcca
[p_G1_Ex5Ex1In1PVALB,r_G1_Ex5Ex1In1PVALB,A_G1,B_G1,U_G1,V_G1] = permcca(gradient_1,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_1_null);

cca_r_p = [r_G1_In1Ex1PVALB,p_G1_In1Ex1PVALB;...
    r_G1_Ex1,p_G1_Ex1;...
    r_G1_In1,p_G1_In1;...
    r_G1_Ex5,p_G1_Ex5;...
    r_G1_PVALB,p_G1_PVALB;...
    r_G1_Ex5PVALB,p_G1_Ex5PVALB;...
    r_G1_Ex1PVALB,p_G1_Ex1PVALB;...
    r_G1_In1PVALB,p_G1_In1PVALB;...
    r_G1_Ex5Ex1In1,p_G1_Ex5Ex1In1;...
    r_G1_Ex5Ex1In1PVALB,p_G1_Ex5Ex1In1PVALB];

RowNames = {'In1Ex1PVALB','Ex1','In1','Ex5','PVALB','Ex5PVALB','Ex1PVALB','In1PVALB',...
    'Ex5Ex1In1','Ex5Ex1In1PVALB'};
ColNames = {'r','p'}
cca_p_r_table = array2table(cca_r_p,'RowNames',RowNames,'VariableNames',ColNames);
writetable(cca_p_r_table,[out_dir,'/ccaGradualRemove_G1_VIS_r_p.csv'],'WriteRowNames',true);

%% Gradient 2 Lake_DFC
% vision-distributed cell types: Ex4, End

% 0. Remove Ex4, End
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[3,12]) = [];
% permcca
[p_G2_Ex4End,r_G2_Ex4End,A_G2,B_G2,U_G2,V_G2] = permcca(gradient_2,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_2_null);

% 1. Remove Ex4
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[3]) = [];
% permcca
[p_G2_Ex4,r_G2_Ex4,A_G2,B_G1,U_G1,V_G1] = permcca(gradient_2,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_2_null);

% 2. Remove End
LakeDFC_mat_reduced = LakeDFC_mat;
LakeDFC_mat_reduced(:,[12]) = [];
% permcca
[p_G2_End,r_G2_End,A_G2,B_G2,U_G2,V_G2] = permcca(gradient_2,LakeDFC_mat_reduced,nP,[],[],[],[],gradient_2_null);

cca_r_p = [r_G2_Ex4End,p_G2_Ex4End;...
    r_G2_Ex4,p_G2_Ex4;...
    r_G2_End,p_G2_End];

RowNames = {'Ex4End','Ex4','End'};
ColNames = {'r','p'};
cca_p_r_table = array2table(cca_r_p,'RowNames',RowNames,'VariableNames',ColNames);
writetable(cca_p_r_table,[out_dir,'/ccaGradualRemove_G2_DFC_r_p.csv'],'WriteRowNames',true);

%% Gradient 2 Lake_VIS
% vision-distributed cell types: Ex4

% 1. Remove Ex4
LakeVIS_mat_reduced = LakeVIS_mat;
LakeVIS_mat_reduced(:,[3]) = [];
% permcca
[p_G2_Ex4,r_G2_Ex4,A_G2,B_G1,U_G1,V_G1] = permcca(gradient_2,LakeVIS_mat_reduced,nP,[],[],[],[],gradient_2_null);

cca_r_p = [r_G2_Ex4,p_G2_Ex4];

RowNames = {'Ex4'};
ColNames = {'r','p'};
cca_p_r_table = array2table(cca_r_p,'RowNames',RowNames,'VariableNames',ColNames);
writetable(cca_p_r_table,[out_dir,'/ccaGradualRemove_G2_VIS_r_p.csv'],'WriteRowNames',true);
