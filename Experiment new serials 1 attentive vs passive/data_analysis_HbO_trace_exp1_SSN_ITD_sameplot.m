clear;
close;

%% load data folder
load('Speech_breath_blockavg_data.mat');
Speech_breath_blockAvg_all = breath_blockAvg_all;
load('Speech_task_blockavg_data.mat');
Speech_task_blockAvg_all = task_blockAvg_all;

load('SSN_breath_blockavg_data.mat');
SSN_breath_blockAvg_all = breath_blockAvg_all;
load('SSN_task_blockavg_data.mat');
SSN_task_blockAvg_all = task_blockAvg_all;

% breath_blockAvg_all = Speech_breath_blockAvg_all;
task_blockAvg_all_ITD = Speech_task_blockAvg_all;
task_blockAvg_all_noITD = SSN_task_blockAvg_all;

%% define basic stuff
% Difine window
window_b=[-5,40];
window_t=[-5,40];

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

figure;
t = (window_t(1):1/50:window_t(2))';
region = {'right cIFS', 'left cIFS', 'right STG', 'left STG'};
all_block_p = [];

Hb_type = 'HbT';
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

data_for_R = [];
for ii = 1 : 4 % only plot C, A, D, B size(breath_blockAvg_all,4)
    
    % each region has different channels selected
    subject_selected = [1:14];
    
    % HbO 1
    breath_site_matrix_HbO_speech = reshape(Speech_breath_blockAvg_all(:, number, 1, ii,subject_selected), [size(Speech_breath_blockAvg_all,1),length(subject_selected)]);
    breath_site_matrix_HbO_ssn = reshape(SSN_breath_blockAvg_all(:, number, 1, ii,subject_selected), [size(SSN_breath_blockAvg_all,1),length(subject_selected)]);
    breath_site_matrix_HbO = (breath_site_matrix_HbO_speech + breath_site_matrix_HbO_ssn)/2;

    ITD_site_matrix_HbO = reshape(task_blockAvg_all_ITD(:, number, 1, ii,subject_selected), [size(task_blockAvg_all_ITD,1),length(subject_selected)]);
    noITD_site_matrix_HbO = reshape(task_blockAvg_all_noITD(:, number, 1, ii,subject_selected), [size(task_blockAvg_all_noITD,1),length(subject_selected)]);
    %diff_site_matrix_HbO = reshape(dcAvg_all_ITD_noITD_diff(:, 1, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
    
%     % HbR 2
%     breath_site_matrix_HbR = reshape(breath_blockAvg_all(:, 2, 1, ii,subject_selected), [size(breath_blockAvg_all,1),length(subject_selected)]);
%     ITD_site_matrix_HbR = reshape(task_blockAvg_all_ITD(:, 2, 1, ii,subject_selected), [size(task_blockAvg_all_ITD,1),length(subject_selected)]);
%     noITD_site_matrix_HbR = reshape(task_blockAvg_all_noITD(:, 2, 1, ii,subject_selected), [size(task_blockAvg_all_noITD,1),length(subject_selected)]);
%     %diff_site_matrix_HbR = reshape(dcAvg_all_ITD_noITD_diff(:, 2, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
%     
%     % HbT 3
%     breath_site_matrix_HbT = reshape(breath_blockAvg_all(:, 3, 1, ii,subject_selected), [size(breath_blockAvg_all,1),length(subject_selected)]);
%     ITD_site_matrix_HbT = reshape(task_blockAvg_all_ITD(:, 3, 1, ii,subject_selected), [size(task_blockAvg_all_ITD,1),length(subject_selected)]);
%     noITD_site_matrix_HbT = reshape(task_blockAvg_all_noITD(:, 3, 1, ii,subject_selected), [size(task_blockAvg_all_noITD,1),length(subject_selected)]);
%     %diff_site_matrix_HbT = reshape(dcAvg_all_ITD_noITD_diff(:, 3, ii,:), [size(dcAvg_all_ITD_noITD_diff,1),size(dcAvg_all_ITD_noITD_diff,4)]);
    
    % use HMS to get functional 
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
    title(['breath bolding at ', region{ii}]);
    
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
    legend('passive','active','Location', 'best');
    
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
    ylim([-1.5,1.5]);
    y_limit=get(gca,'ylim');
    x_limit=get(gca,'xlim');
    line(x_limit,[0 0],'Color','k')
    plot([0 0],y_limit,'r')
    plot([15 15],y_limit,'r')
    set(gca,'TickDir','out');
    title(['task at ', region{ii}]);
    
    
%     % plot noITD - ITD
%     subplot(4,3,3*ii);
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
%     plot([-2.8 -2.8],y_limit,'r')
%     plot([15 15],y_limit,'r')
%     title(['difference at ', region(ii)]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % separate code to save as the format R needed for analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data format:
    % participant hearing resolution Group context(location) time pupil puuil_l
    % ITD_site_matrix_HbO     speech condition
    % noITD_site_matrix_HbO   noise condiiton
    
    
    
    for jj = 1 : size(ITD_site_matrix_HbO,2)
        
        fnirs_data_speech = ITD_site_matrix_HbO(:,jj);
        fnirs_data_noise = noITD_site_matrix_HbO(:,jj);
        
        sid = ones(size(ITD_site_matrix_HbO,1),1) * jj;
        
        % ii 1:E; 2:A; 3:F; 4:B
        % further divide into left/right, and location (STG, ciFS/tgPCs)
        %  left/right: 1:left; 2:right
        %  STG : 1; cIFS: 2
        switch ii
            case 1 % E
                hemsphere = ones(size(ITD_site_matrix_HbO,1),1) * 2; %right
                location = ones(size(ITD_site_matrix_HbO,1),1) * 1; % STG
            case 2 % A
                hemsphere = ones(size(ITD_site_matrix_HbO,1),1) * 1; %left
                location = ones(size(ITD_site_matrix_HbO,1),1) * 1; % STG
            case 3 % F
                hemsphere = ones(size(ITD_site_matrix_HbO,1),1) * 2; %right
                location = ones(size(ITD_site_matrix_HbO,1),1) * 2; % cIFS
            case 4 % B
                hemsphere = ones(size(ITD_site_matrix_HbO,1),1) * 1; %left
                location = ones(size(ITD_site_matrix_HbO,1),1) * 2; % cIFS
        end
        
        %roi_code E,A,F,B: 1,2,3,4
        roi_code = ones(size(ITD_site_matrix_HbO,1),1) * ii;
                 
        % task_type 1:speech; 2: noise
        task_type_speech  = ones(size(ITD_site_matrix_HbO,1),1);
        task_type_noise  = ones(size(ITD_site_matrix_HbO,1),1) * 2;
        
        % time 
        time = -5000 : 20 : 40000;
        
        dummy_speech = cat(2, sid, hemsphere, location, roi_code, task_type_speech, time', fnirs_data_speech);
        dummy_noise = cat(2, sid, hemsphere, location, roi_code, task_type_noise, time', fnirs_data_noise);
        
        data_for_R = cat(1, data_for_R, dummy_speech);
        data_for_R = cat(1, data_for_R, dummy_noise);
    end
    
    
    
end
% save as cvs
csvwrite('R_data_active_passive_HbO.csv',data_for_R)
