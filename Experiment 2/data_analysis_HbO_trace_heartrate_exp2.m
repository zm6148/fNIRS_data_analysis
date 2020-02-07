%% load data folder
% difference between ITD and noITD trace plot
clear all;
%close all;

% path = uigetdir(pwd,'choose reuslts folder');
% filepattern=fullfile(path,'/*.mat');
% files=dir(filepattern);

path = [pwd,'\Data\Processed'];
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);

% % based on experiment type choose which ones to load
% splited = strsplit(path,'\');
% experiment_type = splited{8};
% 
% if strcmp (experiment_type, 'ITD')
%     % load the selected channels
%     selected_channel_folder_name = ['C:\Users\mz86\Desktop\fNIRS analysis\experiment 2\Data\processed\selected channels\', 'selected_channels.mat'];
%     load(selected_channel_folder_name);
% else
%     % load the selected channels
%     selected_channel_folder_name = ['C:\Users\mz86\Desktop\fNIRS analysis\experiment 2\Data\processed\selected channels\', 'selected_channels.mat'];
%     load(selected_channel_folder_name);
% end


%% define basic stuff
% Difine window
window_b=[-5,40];
window_t=[-5,40];

%% load individual data

%(C, A, D, B, CA, DB, CD, AB) * (17 subjects)

breath_blockAvg_all = [];

task_blockAvg_all_ITD= [];

heartrate_blockAvg_all_ITD = [];

task_blockAvg_all_noITD= [];

heartrate_blockAvg_all_noITD = [];

dcAvg_all_ITD_noITD_diff= [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change here for Hb type #####################
Hb_type = 'HbO';
sys_or_fun = 'fun';
window_size = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type = Hb_type; %HbR
switch type
    case 'HbO'
        number = 1;
    case 'HbR'
        number = 2;
    case 'HbT'
        number = 3;
end

for index=1:length(files)
    base_name = files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [path,'\',name,'.mat'];
    load(new_name);
    disp(base_name);
    
    % build Refrence channel regressor
    % for source C RC: 18; channels: 13, 16, 15, 17
    % for source D RC: 19; channels: 20, 23, 21, 24
    % for source A RC: 4; channels: 2, 6, 3, 1
    % for source B RC: 10; channels: 12, 9, 11, 7
    % use regressor closest to the source
    
    % implement HMS before doing any block average
    new_dc = [];
    for ii = 1:size(dc,3)
        [dc_HbO_f, dc_HbR_f] = HMS(dc(:,1, ii), dc(:,2, ii),sys_or_fun);
        temp = cat(2, dc_HbO_f', dc_HbR_f', dc_HbO_f'+ dc_HbR_f');
        new_dc = cat(3, new_dc, temp);
    end
    % assign it back to old dc
    % dc = new_dc;
    
    RC_C = 18;
    channels_C = [13, 16, 15, 17];%, 14];
    % breath task (datapoints * 3 * 4)
    breath_blockAvg_subject_C = hmrBlockAvg(dc(:,:,channels_C),s_b,t,window_t);
    % ITD task
    task_blockAvg_subject_C_ITD= hmrBlockAvg(dc_task(:,:,channels_C),s_ITD_task,t,window_t);
    task_blockAvg_subject_C_ITD = method_normalized_by_b( breath_blockAvg_subject_C, task_blockAvg_subject_C_ITD);
    % noITD task
    task_blockAvg_subject_C_noITD = hmrBlockAvg(dc_task(:,:,channels_C),s_noITD_task,t,window_t);
    task_blockAvg_subject_C_noITD = method_normalized_by_b( breath_blockAvg_subject_C, task_blockAvg_subject_C_noITD);
    % mean heart rate
    mean_heartrate_at_C_ITD= method_heart_rate(dod_raw, s_ITD_task, t, window_t, first_task_index, last_task_index, window_size);
    mean_heartrate_at_C_noITD= method_heart_rate(dod_raw, s_noITD_task, t, window_t, first_task_index, last_task_index, window_size);

    RC_A = 4;
    channels_A = [2, 6, 3, 1];%, 5];
    % breath task
    breath_blockAvg_subject_A = hmrBlockAvg(dc(:,:,channels_A),s_b,t,window_t);
    % ITD task
    task_blockAvg_subject_A_ITD= hmrBlockAvg(dc_task(:,:,channels_A),s_ITD_task,t,window_t);
    task_blockAvg_subject_A_ITD = method_normalized_by_b( breath_blockAvg_subject_A, task_blockAvg_subject_A_ITD);
    % noITD task
    task_blockAvg_subject_A_noITD = hmrBlockAvg(dc_task(:,:,channels_A),s_noITD_task,t,window_t);
    task_blockAvg_subject_A_noITD = method_normalized_by_b( breath_blockAvg_subject_A, task_blockAvg_subject_A_noITD);
    % mean heart rate
    mean_heartrate_at_A_ITD= method_heart_rate(dod_raw, s_ITD_task, t, window_t, first_task_index, last_task_index, window_size);
    mean_heartrate_at_A_noITD= method_heart_rate(dod_raw, s_noITD_task, t, window_t, first_task_index, last_task_index, window_size);

    RC_D = 19;
    channels_D = [20, 23, 21, 24];%, 22];
    % breath task
    breath_blockAvg_subject_D = hmrBlockAvg(dc(:,:,channels_D),s_b,t,window_t);
    % ITD task
    task_blockAvg_subject_D_ITD = hmrBlockAvg(dc_task(:,:,channels_D),s_ITD_task,t,window_t);
    task_blockAvg_subject_D_ITD = method_normalized_by_b( breath_blockAvg_subject_D, task_blockAvg_subject_D_ITD);
    % noITD task
    task_blockAvg_subject_D_noITD = hmrBlockAvg(dc_task(:,:,channels_D),s_noITD_task,t,window_t);
    task_blockAvg_subject_D_noITD = method_normalized_by_b( breath_blockAvg_subject_D, task_blockAvg_subject_D_noITD);
    % mean heart rate
    mean_heartrate_at_D_ITD= method_heart_rate(dod_raw, s_ITD_task, t, window_t, first_task_index, last_task_index, window_size);
    mean_heartrate_at_D_noITD= method_heart_rate(dod_raw, s_noITD_task, t, window_t, first_task_index, last_task_index, window_size);

    RC_B = 10;
    channels_B = [12, 9, 11, 7];%, 8];
    % breath task
    breath_blockAvg_subject_B = hmrBlockAvg(dc(:,:,channels_B),s_b,t,window_t);
    % ITD task
    task_blockAvg_subject_B_ITD= hmrBlockAvg(dc_task(:,:,channels_B),s_ITD_task,t,window_t);
    task_blockAvg_subject_B_ITD = method_normalized_by_b( breath_blockAvg_subject_B, task_blockAvg_subject_B_ITD);
    % noITD task
    task_blockAvg_subject_B_noITD = hmrBlockAvg(dc_task(:,:,channels_B),s_noITD_task,t,window_t);
    task_blockAvg_subject_B_noITD = method_normalized_by_b( breath_blockAvg_subject_B, task_blockAvg_subject_B_noITD);
    % mean heart rate
    mean_heartrate_at_B_ITD= method_heart_rate(dod_raw, s_ITD_task, t, window_t, first_task_index, last_task_index, window_size);
    mean_heartrate_at_B_noITD= method_heart_rate(dod_raw, s_noITD_task, t, window_t, first_task_index, last_task_index, window_size);

%     % combine CA and DB
%     RC_CA = [18,4];
%     channels_CA = [13, 16, 15, 17, 2, 6, 3, 1];
%     % breath task
%     breath_blockAvg_subject_CA = hmrBlockAvg(dc(:,:,channels_CA),s_b,t,window_t);
%     % ITD task
%     task_blockAvg_subject_CA_ITD= hmrBlockAvg(dc_task(:,:,channels_CA),s_ITD_task,t,window_t);
%     % noITD task
%     task_blockAvg_subject_CA_noITD = hmrBlockAvg(dc_task(:,:,channels_CA),s_noITD_task,t,window_t);
%     % ITD no_ITD difference
%     dcAvg_subject_ITD_noITD_diff_CA = method_substraction_by_b(breath_blockAvg_subject_CA, task_blockAvg_subject_CA_ITD, task_blockAvg_subject_CA_noITD);
%     
%     RC_DB = [19, 10];
%     channels_DB = [12, 9, 11, 7, 20, 23, 21, 24];
%     % breath task
%     breath_blockAvg_subject_DB = hmrBlockAvg(dc(:,:,channels_DB),s_b,t,window_t);
%     % ITD task
%     task_blockAvg_subject_DB_ITD= hmrBlockAvg(dc_task(:,:,channels_DB),s_ITD_task,t,window_t);
%     % noITD task
%     task_blockAvg_subject_DB_noITD = hmrBlockAvg(dc_task(:,:,channels_DB),s_noITD_task,t,window_t);
%     % ITD no_ITD difference
%     dcAvg_subject_ITD_noITD_diff_DB = method_substraction_by_b(breath_blockAvg_subject_DB, task_blockAvg_subject_DB_ITD, task_blockAvg_subject_DB_noITD);
%     
%     % combine CD (left) and AB (right)
%     RC_CD = [18, 19];
%     channels_CD = [13, 16, 15, 17, 20, 23, 21, 24];
%     % breath task
%     breath_blockAvg_subject_CD = hmrBlockAvg(dc(:,:,channels_CD),s_b,t,window_t);
%     % ITD task
%     task_blockAvg_subject_CD_ITD= hmrBlockAvg(dc_task(:,:,channels_CD),s_ITD_task,t,window_t);
%     % noITD task
%     task_blockAvg_subject_CD_noITD = hmrBlockAvg(dc_task(:,:,channels_CD),s_noITD_task,t,window_t);
%     % ITD no_ITD difference
%     dcAvg_subject_ITD_noITD_diff_CD = method_substraction_by_b(breath_blockAvg_subject_CD, task_blockAvg_subject_CD_ITD, task_blockAvg_subject_CD_noITD);
%     
%     RC_AB = [4, 10];
%     channels_AB = [2, 6, 3, 1, 12, 9, 11, 7];
%     % breath task
%     breath_blockAvg_subject_AB = hmrBlockAvg(dc(:,:,channels_AB),s_b,t,window_t);
%     % ITD task
%     task_blockAvg_subject_AB_ITD= hmrBlockAvg(dc_task(:,:,channels_AB),s_ITD_task,t,window_t);
%     % noITD task
%     task_blockAvg_subject_AB_noITD = hmrBlockAvg(dc_task(:,:,channels_AB),s_noITD_task,t,window_t);
%     % ITD no_ITD difference
%     dcAvg_subject_ITD_noITD_diff_AB = method_substraction_by_b(breath_blockAvg_subject_AB, task_blockAvg_subject_AB_ITD, task_blockAvg_subject_AB_noITD);
    
    
    breath_blockAvg_subject = cat (4, mean(breath_blockAvg_subject_C,3), ...
        mean(breath_blockAvg_subject_A,3), ...
        mean(breath_blockAvg_subject_D,3), ...
        mean(breath_blockAvg_subject_B,3));
%         mean(breath_blockAvg_subject_CA,3), ...
%         mean(breath_blockAvg_subject_DB,3), ...
%         mean(breath_blockAvg_subject_CD,3), ...
%         mean(breath_blockAvg_subject_AB,3));
    
    task_blockAvg_subject_ITD = cat (4, mean(task_blockAvg_subject_C_ITD,3), ...
        mean(task_blockAvg_subject_A_ITD,3), ...
        mean(task_blockAvg_subject_D_ITD,3), ...
        mean(task_blockAvg_subject_B_ITD,3));
%         mean(task_blockAvg_subject_CA_ITD,3), ...
%         mean(task_blockAvg_subject_DB_ITD,3), ...
%         mean(task_blockAvg_subject_CD_ITD,3), ...
%         mean(task_blockAvg_subject_AB_ITD,3));

    heartrate_blockAvg_subject_ITD = cat (4, mean_heartrate_at_C_ITD, ...
        mean_heartrate_at_A_ITD, ...
        mean_heartrate_at_D_ITD, ...
        mean_heartrate_at_B_ITD);

    task_blockAvg_subject_noITD = cat (4, mean(task_blockAvg_subject_C_noITD,3), ...
        mean(task_blockAvg_subject_A_noITD,3), ...
        mean(task_blockAvg_subject_D_noITD,3), ...
        mean(task_blockAvg_subject_B_noITD,3));
%         mean(task_blockAvg_subject_CA_noITD,3), ...
%         mean(task_blockAvg_subject_DB_noITD,3), ...
%         mean(task_blockAvg_subject_CD_noITD,3), ...
%         mean(task_blockAvg_subject_AB_noITD,3));

    heartrate_blockAvg_subject_noITD = cat (4, mean_heartrate_at_C_noITD, ...
        mean_heartrate_at_A_noITD, ...
        mean_heartrate_at_D_noITD, ...
        mean_heartrate_at_B_noITD);
    
%    dcAvg_subject_ITD_noITD_diff = cat (3, dcAvg_subject_ITD_noITD_diff_C, ...
%         dcAvg_subject_ITD_noITD_diff_A, ...
%         dcAvg_subject_ITD_noITD_diff_D, ...
%         dcAvg_subject_ITD_noITD_diff_B);
%         dcAvg_subject_ITD_noITD_diff_CA, ...
%         dcAvg_subject_ITD_noITD_diff_DB, ...
%         dcAvg_subject_ITD_noITD_diff_CD, ...
%         dcAvg_subject_ITD_noITD_diff_AB);
    
    breath_blockAvg_all = cat(5, breath_blockAvg_all, breath_blockAvg_subject);
    
    task_blockAvg_all_ITD = cat(5, task_blockAvg_all_ITD, task_blockAvg_subject_ITD);

    heartrate_blockAvg_all_ITD = cat(5, heartrate_blockAvg_all_ITD, heartrate_blockAvg_subject_ITD);
    
    task_blockAvg_all_noITD = cat(5, task_blockAvg_all_noITD, task_blockAvg_subject_noITD);

    heartrate_blockAvg_all_noITD = cat(5, heartrate_blockAvg_all_noITD, heartrate_blockAvg_subject_noITD);
    
%    dcAvg_all_ITD_noITD_diff = cat(4, dcAvg_all_ITD_noITD_diff, dcAvg_subject_ITD_noITD_diff);
    
end

%% plot and see

% breath_blockAvg_all
% 2251           3           1           8          17
% data           Hb type     channel     sites      subject

% task_blockAvg_all_ITD; task_blockAvg_all_noITD
% 2251           3           1           8          17
% data           Hb type     channel     sites      subject

% dcAvg_all_ITD_noITD_diff
% 2251           3           8          17
% data           Hb type     sites      subject

% plot by the site of recording
% site order
% C, A, D, B, CA, DB, CD, AB

a = figure;
t = (window_t(1):1/50:window_t(2))';
region = ['C', 'A', 'D', 'B'];
all_block_p = [];
for ii = 1 : 4 % only plot C, A, D, B size(breath_blockAvg_all,4)
    
    % each region has different channels selected
    subject_selected = 1:14;
    
    % HbO 1
    breath_site_matrix_HbO = reshape(breath_blockAvg_all(:, number, 1, ii,subject_selected), [size(breath_blockAvg_all,1),length(subject_selected)]);
    ITD_site_matrix_HbO = reshape(task_blockAvg_all_ITD(:, number, 1, ii,subject_selected), [size(task_blockAvg_all_ITD,1),length(subject_selected)]);
    noITD_site_matrix_HbO = reshape(task_blockAvg_all_noITD(:, number, 1, ii,subject_selected), [size(task_blockAvg_all_noITD,1),length(subject_selected)]);
    %diff_site_matrix_HbO = reshape(dcAvg_all_ITD_noITD_diff(:, 1, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
    
    % HbR 2
    breath_site_matrix_HbR = reshape(breath_blockAvg_all(:, 2, 1, ii,subject_selected), [size(breath_blockAvg_all,1),length(subject_selected)]);
    ITD_site_matrix_HbR = reshape(task_blockAvg_all_ITD(:, 2, 1, ii,subject_selected), [size(task_blockAvg_all_ITD,1),length(subject_selected)]);
    noITD_site_matrix_HbR = reshape(task_blockAvg_all_noITD(:, 2, 1, ii,subject_selected), [size(task_blockAvg_all_noITD,1),length(subject_selected)]);
    %diff_site_matrix_HbR = reshape(dcAvg_all_ITD_noITD_diff(:, 2, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
    
    % HbT 3
    breath_site_matrix_HbT = reshape(breath_blockAvg_all(:, 3, 1, ii,subject_selected), [size(breath_blockAvg_all,1),length(subject_selected)]);
    ITD_site_matrix_HbT = reshape(task_blockAvg_all_ITD(:, 3, 1, ii,subject_selected), [size(task_blockAvg_all_ITD,1),length(subject_selected)]);
    noITD_site_matrix_HbT = reshape(task_blockAvg_all_noITD(:, 3, 1, ii,subject_selected), [size(task_blockAvg_all_noITD,1),length(subject_selected)]);
    %diff_site_matrix_HbT = reshape(dcAvg_all_ITD_noITD_diff(:, 3, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
    
%     % use HMS to get functional 
%     [ITD_site_matrix_HbO_f, ITD_site_matrix_HbR_f] = HMS_blvg(ITD_site_matrix_HbO, ITD_site_matrix_HbR);
%     [noITD_site_matrix_HbO_f, noITD_site_matrix_HbR_f] = HMS_blvg(noITD_site_matrix_HbO, noITD_site_matrix_HbR);
%     [breath_site_matrix_HbO_f, breath_site_matrix_HbR_f] = HMS_blvg(breath_site_matrix_HbO, breath_site_matrix_HbR);
    
    % before plotting do t-test on normalized concetration on 3 sections
    % 1: 5s before onset; 2: 5~20s while task in on; 3: 20~25s after the
    % task
    block_1_index = 1:5*50; %5s and 50 Hz sampling frequency
    block_2_index = (5*50+1) : 35*50;
    block_3_index = (35*50+1): 45*50;
    % select ITD noITD data based on block index
    % ITD
    ITD_site_matrix_HbO_block_1 = ITD_site_matrix_HbO(block_1_index,:);
    ITD_site_matrix_HbO_block_2 = ITD_site_matrix_HbO(block_2_index,:);
    ITD_site_matrix_HbO_block_3 = ITD_site_matrix_HbO(block_3_index,:);
    % noITD
    noITD_site_matrix_HbO_block_1 = noITD_site_matrix_HbO(block_1_index,:);
    noITD_site_matrix_HbO_block_2 = noITD_site_matrix_HbO(block_2_index,:);
    noITD_site_matrix_HbO_block_3 = noITD_site_matrix_HbO(block_3_index,:);
    % perform t-test between ITD and noITD data on all 3 sections of
    % recording
    % ITD
    ITD_site_matrix_HbO_block_1_average = mean(ITD_site_matrix_HbO_block_1,2);
    ITD_site_matrix_HbO_block_2_average = mean(ITD_site_matrix_HbO_block_2,2);
    ITD_site_matrix_HbO_block_3_average = mean(ITD_site_matrix_HbO_block_3,2);
    % noITD
    noITD_site_matrix_HbO_block_1_average = mean(noITD_site_matrix_HbO_block_1,2);
    noITD_site_matrix_HbO_block_2_average = mean(noITD_site_matrix_HbO_block_2,2);
    noITD_site_matrix_HbO_block_3_average = mean(noITD_site_matrix_HbO_block_3,2);
    % paired ttest
    [h,p1,ci,stats] = ttest2(ITD_site_matrix_HbO_block_1_average, noITD_site_matrix_HbO_block_1_average);
    [h,p2,ci,stats] = ttest2(ITD_site_matrix_HbO_block_2_average, noITD_site_matrix_HbO_block_2_average);
    [h,p3,ci,stats] = ttest2(ITD_site_matrix_HbO_block_3_average, noITD_site_matrix_HbO_block_3_average);
    all_block_p = cat(1, all_block_p, [p1,p2,p3]);
    
    % start plotting
    % breath hold trace
%     subplot(4,2,2*ii-1);
%     hold on;
%     all_average = mean(breath_site_matrix_HbO,2);
%     SE = std(breath_site_matrix_HbO,0,2);
%     plot(t,all_average,'b','LineWidth',2);
%     %calcualte SE and plot error bar
%     % ITD with err bar
%     err = SE/sqrt(size(breath_site_matrix_HbO,2));
%     h = fill([t;flipud(t)],[all_average-err;flipud(all_average+err)],[0 0 1],'linestyle','none');
%     set(h,'facealpha',.3);
%     % add lines
%     xlim([-5,40])
%     ylim([-1*10^(-6) 1*10^(-6)]);
%     y_limit=get(gca,'ylim');
%     x_limit=get(gca,'xlim');
%     line(x_limit,[0 0],'Color','k')
%     plot([0 0],y_limit,'r')
%     plot([15 15],y_limit,'r')
%     set(gca,'TickDir','out');
%     title(['breath bolding at ', region(ii)]);
%     
%     % ITD and no_ITD trace
%     subplot(4,2,2*ii);
%     hold on;
% 
%     % noITD
%     noITD_average = mean(noITD_site_matrix_HbO,2);
%     noITD_SE = std(noITD_site_matrix_HbO,0,2);
%     % ITD
%     ITD_average = mean(ITD_site_matrix_HbO,2);
%     ITD_SE = std(ITD_site_matrix_HbO,0,2);
%     
%     plot(t,noITD_average,'b','LineWidth',2);
%     plot(t,ITD_average,'r','LineWidth',2);
%     legend('pitch only','pitch and space','Location', 'best');
%     
%     %calcualte SE and plot error bar
%     % ITD with err bar
%     noITD_err = noITD_SE/sqrt(size(noITD_site_matrix_HbO,2));
%     h = fill([t;flipud(t)],[noITD_average-noITD_err;flipud(noITD_average+noITD_err)],[0 0 1],'linestyle','none');
%     set(h,'facealpha',.3);
%  
%     %calcualte SE and plot error bar
%     % ITD with err bar
%     ITD_err = ITD_SE/sqrt(size(ITD_site_matrix_HbO,2));
%     h = fill([t;flipud(t)],[ITD_average-ITD_err;flipud(ITD_average+ITD_err)],[1 0 0],'linestyle','none');
%     set(h,'facealpha',.3);
% 
%     % add lines
%     xlim([-5,40]);
%     ylim([-0.8,0.8]);
%     y_limit=get(gca,'ylim');
%     x_limit=get(gca,'xlim');
%     line(x_limit,[0 0],'Color','k')
%     plot([0 0],y_limit,'r')
%     plot([15 15],y_limit,'r')
%     set(gca,'TickDir','out');
%     title(['task at ', region(ii)]);
    

    % NEW HMS start plotting
    % breath hold trace
    subplot(4,2,2*ii-1);
    hold on;
    all_average = mean(breath_site_matrix_HbO,2);
    SE = std(breath_site_matrix_HbO,0,2);
    plot(t,all_average,'b','LineWidth',2);
    %calcualte SE and plot error bar
    % ITD with err bar
    err = SE/sqrt(size(breath_site_matrix_HbO,2));
    h = fill([t;flipud(t)],[all_average-err;flipud(all_average+err)],[0 0 1],'linestyle','none');
    set(h,'facealpha',.3);
    % add lines
    xlim([-5,40])
    ylim([-1*10^(-6) 1*10^(-6)]);
    y_limit=get(gca,'ylim');
    x_limit=get(gca,'xlim');
    line(x_limit,[0 0],'Color','k')
    plot([0 0],y_limit,'r')
    plot([15 15],y_limit,'r')
    set(gca,'TickDir','out');
    title(['breath bolding at ', region(ii)]);
    
    % ITD and no_ITD trace
    subplot(4,2,2*ii);
    hold on;
    % noITD
    noITD_average = mean(noITD_site_matrix_HbO,2);
    noITD_SE = std(noITD_site_matrix_HbO,0,2);
    % ITD
    ITD_average = mean(ITD_site_matrix_HbO,2);
    ITD_SE = std(ITD_site_matrix_HbO,0,2);
    
    plot(t,noITD_average,'b','LineWidth',2);
    plot(t,ITD_average,'r','LineWidth',2);
    legend('pitch only','pitch and space','Location', 'best');
    
    %calcualte SE and plot error bar
    % ITD with err bar
    noITD_err = noITD_SE/sqrt(size(noITD_site_matrix_HbO,2));
    h = fill([t;flipud(t)],[noITD_average-noITD_err;flipud(noITD_average+noITD_err)],[0 0 1],'linestyle','none');
    set(h,'facealpha',.3);
 
    %calcualte SE and plot error bar
    % ITD with err bar
    ITD_err = ITD_SE/sqrt(size(ITD_site_matrix_HbO,2));
    h = fill([t;flipud(t)],[ITD_average-ITD_err;flipud(ITD_average+ITD_err)],[1 0 0],'linestyle','none');
    set(h,'facealpha',.3);

    % add lines
    xlim([-5,40]);
    %ylim([-0.8,0.8]);
    y_limit=get(gca,'ylim');
    x_limit=get(gca,'xlim');
    line(x_limit,[0 0],'Color','k')
    plot([0 0],y_limit,'r')
    plot([15 15],y_limit,'r')
    set(gca,'TickDir','out');
    title(['task at ', region(ii)]);
end
%print(a,'-fillpage','substraction_exp2','-dpdf','-r0')


b = figure;
t = (window_t(1):1/50:window_t(2))';
region = ['C', 'A', 'D', 'B'];
all_block_p = [];
for ii = 1 : 4 % only plot C, A, D, B size(breath_blockAvg_all,4)
    
    % each region has different channels selected
    subject_selected = 1:14;
    
    % HbO 1
    ITD_site_matrix_HbO = reshape(heartrate_blockAvg_all_ITD(:, 1, 1, ii,subject_selected), [size(heartrate_blockAvg_all_ITD,1),length(subject_selected)]);
    noITD_site_matrix_HbO = reshape(heartrate_blockAvg_all_noITD(:, 1, 1, ii,subject_selected), [size(heartrate_blockAvg_all_ITD,1),length(subject_selected)]);
    %diff_site_matrix_HbO = reshape(dcAvg_all_ITD_noITD_diff(:, 1, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
     
    % ITD and no_ITD trace
    subplot(4,1,ii);
    hold on;
    % noITD
    noITD_average = mean(noITD_site_matrix_HbO,2);
    noITD_SE = std(noITD_site_matrix_HbO,0,2);
    % ITD
    ITD_average = mean(ITD_site_matrix_HbO,2);
    ITD_SE = std(ITD_site_matrix_HbO,0,2);
    
    plot(t,noITD_average,'b','LineWidth',2);
    plot(t,ITD_average,'r','LineWidth',2);
    legend('pitch only','pitch and space','Location', 'best');
    
    %calcualte SE and plot error bar
    % ITD with err bar
    noITD_err = noITD_SE/sqrt(size(noITD_site_matrix_HbO,2));
    h = fill([t;flipud(t)],[noITD_average-noITD_err;flipud(noITD_average+noITD_err)],[0 0 1],'linestyle','none');
    set(h,'facealpha',.3);
 
    %calcualte SE and plot error bar
    % ITD with err bar
    ITD_err = ITD_SE/sqrt(size(ITD_site_matrix_HbO,2));
    h = fill([t;flipud(t)],[ITD_average-ITD_err;flipud(ITD_average+ITD_err)],[1 0 0],'linestyle','none');
    set(h,'facealpha',.3);

    % add lines
    xlim([-5,40]);
    %ylim([-0.8,0.8]);
    y_limit=get(gca,'ylim');
    x_limit=get(gca,'xlim');
    line(x_limit,[0 0],'Color','k')
    plot([0 0],y_limit,'r')
    plot([15 15],y_limit,'r')
    set(gca,'TickDir','out');
    title(['task at ', region(ii)]);
end