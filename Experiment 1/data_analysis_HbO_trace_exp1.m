%% load data folder
% analysis using block average trace during task
clear;
%close all;

% path = uigetdir(pwd,'choose reuslts folder');
% filepattern=fullfile(path,'/*.mat');
% files=dir(filepattern);

% path = [pwd,'\Data\Speech\Processed'];
% filepattern=fullfile(path,'/*.mat');
% files=dir(filepattern);

path = [pwd,'\Data\SSN\Processed'];
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);

% based on experiment type choose which ones to load
splited = strsplit(path,'\');
experiment_type = splited{8};

%%  use selected subject only

% intersect of SSN and ITD selected Channels
% selected_subject = ['C:\Users\mz86\Desktop\fNIRS analysis\experiment 1\new Control data\', 'selected_channels_inter'];
% selected_subject = ['C:\Users\mz86\Desktop\fNIRS analysis\experiment 1\data\', 'selected_channels_inter'];
% load(selected_subject);
% both_selected = selected_subject_inter;

%% for only breath holding
%selected_channel_folder_name = 'C:\Users\mz86\Desktop\fNIRS analysis\experiment 1\new Control data\all\selected channels\selected_channels_SSN.mat';
%load(selected_channel_folder_name);
% load('C:\Users\mz86\Desktop\fNIRS analysis\experiment 1\Data\All Subjects Breath Only\selected channels\selected_channels_breathonly.mat');


%% define basic stuff
% Difine window
window_b=[-5,40];
window_t=[-5,40];

%% load individual data

%(C, A, D, B, CA, DB, CD, AB) * (17 subjects)

breath_blockAvg_all = [];

task_blockAvg_all= [];

% task_blockAvg_all_noITD= [];
% 
% dcAvg_all_ITD_noITD_diff= [];

Hb_type = 'HbO';
sys_or_fun = 'fun';

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
        [dc_HbO_f, dc_HbR_f] = HMS(dc(:,1, ii), dc(:,2, ii), sys_or_fun);
        temp = cat(2, dc_HbO_f', dc_HbR_f', dc_HbO_f'+ dc_HbR_f');
        new_dc = cat(3, new_dc, temp);
    end
    % assign it back to old dc
    % dc = new_dc;

    
    s_task = s_ITD_task + s_noITD_task;
    
    RC_C = 18;
    channels_C = [13, 16, 15, 17];%, 14];
    % breath task (datapoints * 3 * 4)
    breath_blockAvg_subject_C = hmrBlockAvg(dc(:,:,channels_C),s_b,t,window_t);
    % task
    task_blockAvg_subject_C= hmrBlockAvg(dc_task(:,:,channels_C),s_task,t,window_t);
    % mean 
    mean_at_C = method_normalized_by_b(breath_blockAvg_subject_C, task_blockAvg_subject_C);
   
    
    RC_A = 4;
    channels_A = [2, 6, 3, 1];%, 5];
    % breath task
    breath_blockAvg_subject_A = hmrBlockAvg(dc(:,:,channels_A),s_b,t,window_t);
    % ITD task
    task_blockAvg_subject_A= hmrBlockAvg(dc_task(:,:,channels_A),s_task,t,window_t);
    % mean 
    mean_at_A = method_normalized_by_b(breath_blockAvg_subject_A, task_blockAvg_subject_A);

    
    RC_D = 19;
    channels_D = [20, 23, 21, 24];%, 22];
    % breath task
    breath_blockAvg_subject_D = hmrBlockAvg(dc(:,:,channels_D),s_b,t,window_t);
    % ITD task
    task_blockAvg_subject_D= hmrBlockAvg(dc_task(:,:,channels_D),s_task,t,window_t);
    % mean 
    mean_at_D = method_normalized_by_b(breath_blockAvg_subject_D, task_blockAvg_subject_D);

    RC_B = 10;
    channels_B = [12, 9, 11, 7];%, 8];
    % breath task
    breath_blockAvg_subject_B = hmrBlockAvg(dc(:,:,channels_B),s_b,t,window_t);
    % ITD task
    task_blockAvg_subject_B= hmrBlockAvg(dc_task(:,:,channels_B),s_task,t,window_t);
    % mean 
    mean_at_B = method_normalized_by_b(breath_blockAvg_subject_B, task_blockAvg_subject_B);

%     
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
%     
    
    breath_blockAvg_subject = cat (4, mean(breath_blockAvg_subject_C,3), ...
        mean(breath_blockAvg_subject_A,3), ...
        mean(breath_blockAvg_subject_D,3), ...
        mean(breath_blockAvg_subject_B,3));
%         mean(breath_blockAvg_subject_CA,3), ...
%         mean(breath_blockAvg_subject_DB,3), ...
%         mean(breath_blockAvg_subject_CD,3), ...
%         mean(breath_blockAvg_subject_AB,3));
    
    task_blockAvg_subject = cat (4, mean_at_C, ...
        mean_at_A, ...
        mean_at_D, ...
        mean_at_B);
%         mean(task_blockAvg_subject_CA_ITD,3), ...
%         mean(task_blockAvg_subject_DB_ITD,3), ...
%         mean(task_blockAvg_subject_CD_ITD,3), ...
%         mean(task_blockAvg_subject_AB_ITD,3));
    
%     task_blockAvg_subject_noITD = cat (4, mean(task_blockAvg_subject_C_noITD,3), ...
%         mean(task_blockAvg_subject_A_noITD,3), ...
%         mean(task_blockAvg_subject_D_noITD,3), ...
%         mean(task_blockAvg_subject_B_noITD,3), ...
%         mean(task_blockAvg_subject_CA_noITD,3), ...
%         mean(task_blockAvg_subject_DB_noITD,3), ...
%         mean(task_blockAvg_subject_CD_noITD,3), ...
%         mean(task_blockAvg_subject_AB_noITD,3));
%     
%     dcAvg_subject_ITD_noITD_diff = cat (3, dcAvg_subject_ITD_noITD_diff_C, ...
%         dcAvg_subject_ITD_noITD_diff_A, ...
%         dcAvg_subject_ITD_noITD_diff_D, ...
%         dcAvg_subject_ITD_noITD_diff_B, ...
%         dcAvg_subject_ITD_noITD_diff_CA, ...
%         dcAvg_subject_ITD_noITD_diff_DB, ...
%         dcAvg_subject_ITD_noITD_diff_CD, ...
%         dcAvg_subject_ITD_noITD_diff_AB);
    
    breath_blockAvg_all = cat(5, breath_blockAvg_all, breath_blockAvg_subject);
    
    task_blockAvg_all = cat(5, task_blockAvg_all, task_blockAvg_subject);
    
%     task_blockAvg_all_noITD = cat(5, task_blockAvg_all_noITD, task_blockAvg_subject_noITD);
%     
%     dcAvg_all_ITD_noITD_diff = cat(4, dcAvg_all_ITD_noITD_diff, dcAvg_subject_ITD_noITD_diff);
    
end

% save('Speech_breath_blockavg_data.mat', 'breath_blockAvg_all');
% save('Speech_task_blockavg_data.mat', 'task_blockAvg_all');

save('SSN_breath_blockavg_data.mat', 'breath_blockAvg_all');
save('SSN_task_blockavg_data.mat', 'task_blockAvg_all');


% %% plot and see

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

for ii = 1 : size(breath_blockAvg_all,4)
    
    % each region has different channels selected
    subject_selected = [1:14]; %both_selected{ii};
  
    % HbO 1
    breath_site_matrix_HbO = reshape(breath_blockAvg_all(:, number, 1, ii,subject_selected), [size(breath_blockAvg_all,1),length(subject_selected)]);
    task_site_matrix_HbO = reshape(task_blockAvg_all(:, number, 1, ii,subject_selected), [size(task_blockAvg_all,1),length(subject_selected)]);
%     noITD_site_matrix_HbO = reshape(task_blockAvg_all_noITD(:, 1, 1, ii,:), [size(task_blockAvg_all_noITD,1),size(task_blockAvg_all_noITD,5)]);
%     diff_site_matrix_HbO = reshape(dcAvg_all_ITD_noITD_diff(:, 1, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
%     
    % HbR 2
    breath_site_matrix_HbR = reshape(breath_blockAvg_all(:, 2, 1, ii,subject_selected), [size(breath_blockAvg_all,1),length(subject_selected)]);
    task_site_matrix_HbR = reshape(task_blockAvg_all(:, 2, 1, ii,subject_selected), [size(task_blockAvg_all,1),length(subject_selected)]);
%     noITD_site_matrix_HbR = reshape(task_blockAvg_all_noITD(:, 2, 1, ii,:), [size(task_blockAvg_all_noITD,1),size(task_blockAvg_all_noITD,5)]);
%     diff_site_matrix_HbR = reshape(dcAvg_all_ITD_noITD_diff(:, 2, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
%     
    % HbT 3
    breath_site_matrix_HbT = reshape(breath_blockAvg_all(:, 3, 1, ii,subject_selected), [size(breath_blockAvg_all,1),length(subject_selected)]);
    task_site_matrix_HbT = reshape(task_blockAvg_all(:, 3, 1, ii,subject_selected), [size(task_blockAvg_all,1),length(subject_selected)]);
%     noITD_site_matrix_HbT = reshape(task_blockAvg_all_noITD(:, 3, 1, ii,:), [size(task_blockAvg_all_noITD,1),size(task_blockAvg_all_noITD,5)]);
%     diff_site_matrix_HbT = reshape(dcAvg_all_ITD_noITD_diff(:, 3, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
%     
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
    
    % task trace
    subplot(4,2,2*ii);
    hold on;
    all_average = mean(task_site_matrix_HbO,2);
    SE = std(task_site_matrix_HbO,0,2);
    plot(t,all_average,'r','LineWidth',2);
    %calcualte SE and plot error bar
    % ITD with err bar
    err = SE/sqrt(size(task_site_matrix_HbO,2));
    h = fill([t;flipud(t)],[all_average-err;flipud(all_average+err)],[1 0 0],'linestyle','none');
    set(h,'facealpha',.3);
    % add lines
    xlim([-5,40])
    ylim([-0.8,0.8]);
    y_limit=get(gca,'ylim');
    x_limit=get(gca,'xlim');
    line(x_limit,[0 0],'Color','k')
    plot([-2.8 -2.8],y_limit,'r')
    plot([15 15],y_limit,'r')
    set(gca,'TickDir','out');
    title(['Task at ', region(ii)]);
    
%     % difference trace
%     subplot(8,3,3*ii);
%     hold on;
%     all_average = mean(diff_site_matrix_HbO,2);
%     SE = std(diff_site_matrix_HbO,0,2);
%     plot(t,all_average,'b','LineWidth',2);
%     %calcualte SE and plot error bar
%     % ITD - noITD
%     err = SE/sqrt(size(diff_site_matrix_HbR,2));
%     h = fill([t;flipud(t)],[all_average-err;flipud(all_average+err)],[0 0 1],'linestyle','none');
%     set(h,'facealpha',.3);
%     % add lines
%     xlim([-5,40])
%     y_limit=get(gca,'ylim');
%     x_limit=get(gca,'xlim');
%     line(x_limit,[0 0],'Color','k')
%     plot([0 0],y_limit,'r')
%     plot([15 15],y_limit,'r')
end
