clear;
close;

%% load data folder experiment 1
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


%% load data folder experiment 3
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


%% select which type of blood
Hb_type = 'HbO';

type = Hb_type; %HbR
switch type
    case 'HbO'
        number = 1;
    case 'HbR'
        number = 2;
    case 'HbT'
        number = 3;
end

%% Difine window
window_b=[-5,40];
window_t=[-5,40];

%% for each time point do t-test on it
t_each_ROI_time = [];
p_each_ROI_time = [];

for time = 1 : size(task_blockAvg_all_ITD, 1)
    % for each region of ROI
    t_each_ROI = [];
    p_each_ROI = [];
    
    for ROI = 1:4
        Speech_task_at_time = Speech_task_blockAvg_all(time, number,1, ROI,:);
        SSN_task_at_time = SSN_task_blockAvg_all(time, number,1, ROI,:);
        
        A = reshape(Speech_task_at_time, [1,14]);
        B = reshape(SSN_task_at_time, [1,14]);
        [t, p] = t_test(A, B);
        
        t_each_ROI(ROI) = t;
        p_each_ROI(ROI) = p;
    end
    
    t_each_ROI_time = cat(1, t_each_ROI_time, t_each_ROI);
    p_each_ROI_time = cat(1, p_each_ROI_time, p_each_ROI);
end


%% plot
region = {'right cIFS', 'left cIFS', 'right STG', 'left STG'};
for ROI = 1 : 4 % plot E,A,F,B
    
    time = (window_t(1):1/50:window_t(2))';
    subplot(2,2,ROI);
    hold on;
    
    plot(time, t_each_ROI_time(:,ROI));
    
    % find index where p < 0.05
    significant_index = find(p_each_ROI_time(:,ROI) < 0.05);
    
    y_limit=get(gca,'ylim');
    start_index = time(significant_index(1));
    end_index = time(significant_index(length(significant_index)));
    plot([start_index start_index],y_limit,'r')
    plot([end_index end_index],y_limit,'r')

    xlim([-5,40]);
    title(['t-test at ', region{ROI}]);
    
end

