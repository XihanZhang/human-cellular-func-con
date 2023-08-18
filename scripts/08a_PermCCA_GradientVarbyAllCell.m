
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
    LakeDFC_cell_type_list{i}=strrep(LakeDFC_cell_type_list{i},'DFC_','')
    LakeVIS_cell_type_list{i}=strrep(LakeVIS_cell_type_list{i},'VIS_','')
end

% 4. drop NaN parcels
nan_rows = find(isnan(LakeDFC_mat(:,1)));
LakeDFC_mat(nan_rows,:) = [];
LakeVIS_mat(nan_rows,:) = [];
gradient_1(nan_rows,:) = [];
gradient_1_null(nan_rows,:) = [];
gradient_2(nan_rows,:) = [];
gradient_2_null(nan_rows,:) = [];

%% permulation CCA with all 17 cells
nP = 1000;

[p_G1_DFC,r_G1_DFC,A_G1_DFC,B_G1_DFC,U_G1_DFC,V_G1_DFC] = permcca(gradient_1,LakeDFC_mat,nP,[],[],[],[],gradient_1_null);
[p_G1_VIS,r_G1_VIS,A_G1_VIS,B_G1_VIS,U_G1_VIS,V_G1_VIS] = permcca(gradient_1,LakeVIS_mat,nP,[],[],[],[],gradient_1_null);

loading_G1_DFC=zeros(cell_num,2);
loading_G1_VIS=zeros(cell_num,2);
for i=1:cell_num
    [loading_G1_DFC(i,1),loading_G1_DFC(i,2)]=corr(LakeDFC_mat(:,i),U_G1_DFC);
    [loading_G1_VIS(i,1),loading_G1_VIS(i,2)]=corr(LakeVIS_mat(:,i),U_G1_VIS);
end

mdl = fitlm(U_G1_DFC,V_G1_DFC);
figure('Position', [10 10 1300 400]);
subplot(1,2,1)
bar(loading_G1_DFC(:,1));
set(gca, 'XTick', 1:cell_num, 'XTickLabels', LakeDFC_cell_type_list);
for ii=1:cell_num
    if loading_G1_DFC(ii,2)<(0.05/cell_num)
        if loading_G1_DFC(ii,1)>0
            text(ii, loading_G1_DFC(ii,1)+0.015, '*');
        else
            text(ii, loading_G1_DFC(ii,1)-0.015, '*');
        end
    end
end

subplot(1,2,2)
scatter(U_G1_DFC,V_G1_DFC);
plot(mdl);
xlabel('U');
ylabel('V');
title(['Gradient1 r=',num2str(r_G1_DFC),', p=',num2str(p_G1_DFC),' (Lake DFC)']);
saveas(gcf,[figure_dir,'/LakeDFC_G1.png'])

%--
mdl = fitlm(U_G1_VIS,V_G1_VIS);
figure('Position', [10 10 1300 400]);
subplot(1,2,1)
bar(loading_G1_VIS(:,1));
set(gca, 'XTick', 1:cell_num, 'XTickLabels', LakeVIS_cell_type_list);
for ii=1:cell_num
    if loading_G1_VIS(ii,2)<(0.05/cell_num)
        if loading_G1_VIS(ii,1)>0
            text(ii, loading_G1_VIS(ii,1)+0.015, '*');
        else
            text(ii, loading_G1_VIS(ii,1)-0.015, '*');
        end
    end
end

subplot(1,2,2)
scatter(U_G1_VIS,V_G1_VIS);
plot(mdl);
xlabel('U');
ylabel('V');
title(['Gradient1 r=',num2str(r_G1_VIS),', p=',num2str(p_G1_VIS),' (Lake VIS)']);
saveas(gcf,[figure_dir,'/LakeVIS_G1.png'])

%% Gradient 2
[p_G2_DFC,r_G2_DFC,A_G2_DFC,B_G2_DFC,U_G2_DFC,V_G2_DFC] = permcca(gradient_2,LakeDFC_mat,nP,[],[],[],[],gradient_2_null);
[p_G2_VIS,r_G2_VIS,A_G2_VIS,B_G2_VIS,U_G2_VIS,V_G2_VIS] = permcca(gradient_2,LakeVIS_mat,nP,[],[],[],[],gradient_2_null);

loading_G2_DFC=zeros(cell_num,2);
loading_G2_VIS=zeros(cell_num,2);
for i=1:cell_num
    [loading_G2_DFC(i,1),loading_G2_DFC(i,2)]=corr(LakeDFC_mat(:,i),U_G2_DFC);
    [loading_G2_VIS(i,1),loading_G2_VIS(i,2)]=corr(LakeVIS_mat(:,i),U_G2_VIS);
end

mdl = fitlm(U_G2_DFC,V_G2_DFC);
figure('Position', [10 10 1300 400]);
subplot(1,2,1)
bar(loading_G2_DFC(:,1));
set(gca, 'XTick', 1:cell_num, 'XTickLabels', LakeDFC_cell_type_list);
for ii=1:cell_num
    if loading_G2_DFC(ii,2)<(0.05/cell_num)
        if loading_G2_DFC(ii,1)>0
            text(ii, loading_G2_DFC(ii,1)+0.015, '*');
        else
            text(ii, loading_G2_DFC(ii,1)-0.015, '*');
        end
    end
end

subplot(1,2,2)
scatter(U_G2_DFC,V_G2_DFC);
plot(mdl);
xlabel('U');
ylabel('V');
title(['Gradient2 r=',num2str(r_G2_DFC),', p=',num2str(p_G2_DFC),' (Lake DFC)']);
saveas(gcf,[figure_dir,'/LakeDFC_G2.png'])

%--
mdl = fitlm(U_G2_VIS,V_G2_VIS);
figure('Position', [10 10 1300 400]);
subplot(1,2,1)
bar(loading_G2_VIS(:,1));
set(gca, 'XTick', 1:cell_num, 'XTickLabels', LakeVIS_cell_type_list);
for ii=1:cell_num
    if loading_G2_VIS(ii,2)<(0.05/cell_num)
        if loading_G2_VIS(ii,1)>0
            text(ii, loading_G2_VIS(ii,1)+0.015, '*');
        else
            text(ii, loading_G2_VIS(ii,1)-0.015, '*');
        end
    end
end

subplot(1,2,2)
scatter(U_G2_VIS,V_G2_VIS);
plot(mdl);
xlabel('U');
ylabel('V');
title(['Gradient2 r=',num2str(r_G2_VIS),', p=',num2str(p_G2_VIS),' (Lake VIS)']);
saveas(gcf,[figure_dir,'/LakeVIS_G2.png'])

%% Saving
% 1. output everything
save([out_dir,'/permcca_Gradients_Celltypes.mat'],...,
    'p_G1_DFC','r_G1_DFC','A_G1_DFC','B_G1_DFC','U_G1_DFC','V_G1_DFC',...
    'p_G1_VIS','r_G1_VIS','A_G1_VIS','B_G1_VIS','U_G1_VIS','V_G1_VIS',...
    'p_G2_DFC','r_G2_DFC','A_G2_DFC','B_G2_DFC','U_G2_DFC','V_G2_DFC',...
    'p_G2_VIS','r_G2_VIS','A_G2_VIS','B_G2_VIS','U_G2_VIS','V_G2_VIS',...
    'loading_G1_DFC','loading_G1_VIS','loading_G2_DFC','loading_G2_VIS',...
    'LakeDFC_cell_type_list','LakeVIS_cell_type_list');

load([out_dir,'/permcca_Gradients_Celltypes.mat']);

% 2. make table for loadings and p-values
loadings_DFC = [loading_G1_DFC,loading_G2_DFC];
ColNames = {'loading_G1';'loading_G1_p';...
    'loading_G2';'loading_G2_p'}
loading_DFC_table = array2table(loadings_DFC,'RowNames',LakeDFC_cell_type_list,'VariableNames',ColNames);
writetable(loading_DFC_table,[out_dir,'/loading_DFC_table.csv'],'WriteRowNames',true);

loadings_VIS = [loading_G1_VIS,loading_G2_VIS];
ColNames = {'loading_G1';'loading_G1_p';...
    'loading_G2';'loading_G2_p'}
loading_VIS_table = array2table(loadings_VIS,'RowNames',LakeVIS_cell_type_list,'VariableNames',ColNames);
writetable(loading_VIS_table,[out_dir,'/loading_VIS_table.csv'],'WriteRowNames',true);

% 3. make table for U and V
U_V_all = [U_G1_DFC,V_G1_DFC,U_G1_VIS,V_G1_VIS,...
    U_G2_DFC,V_G2_DFC,U_G2_VIS,V_G2_VIS];
ColNames = {'U_G1_DFC','V_G1_DFC','U_G1_VIS','V_G1_VIS',...
    'U_G2_DFC','V_G2_DFC','U_G2_VIS','V_G2_VIS'};
U_V_table = array2table(U_V_all,'VariableNames',ColNames);
writetable(U_V_table,[out_dir,'/U_V_table.csv']);

% 4. make table for p-value and r for CCA
cca_p_r_all = [p_G1_DFC,p_G1_VIS,p_G2_DFC,p_G2_VIS;...
    r_G1_DFC,r_G1_VIS,r_G2_DFC,r_G2_VIS];
ColNames = {'G1_DFC','G1_VIS','G2_DFC','G2_VIS'};
RowNames = {'p','r'}
cca_p_r_table = array2table(cca_p_r_all,'RowNames',RowNames,'VariableNames',ColNames);
writetable(cca_p_r_table,[out_dir,'/cca_p_r_table.csv']);
