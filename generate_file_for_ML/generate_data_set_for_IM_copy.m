% generate data set for ML on colab
clear;
close all;
%% define basic stuff
% Difine window
window_b=[-5,40];
window_t=[-5,40];
Fs = 50;

type = 'HbO'; %HbR, HbT
switch type
    case 'HbO'
        number = 1;
    case 'HbR'
        number = 2;
    case 'HbT'
        number = 3;
end

ROI_list = [];
task_type_list = [];
blvg_all_subject = [];
block_left_cIFS_HbO = [];

%% active data marked as 1
% for exp1 data only take region cIFS
path = [pwd,'\data\attentive_1'];
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);

for index=1:length(files)
    % load file and display name
    base_name = files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [path,'\',name,'.mat'];
    load(new_name);
    disp(base_name);
    splited = strsplit(base_name,'_');
    experiment_type = splited{2};
    
    % task marks
    s_task = s_ITD_task + s_noITD_task;
    
    % which experiment
    if strcmp(experiment_type,'ITD')
        % location code 1
        left_cIFS = [20, 23, 21, 24];
        % location code 2
        right_cIFS = [12, 9, 11, 7];
        % location code 3
        left_STG = [];
        % location code 4
        right_STG = [];
    else
        % location code 1
        left_cIFS = [1,2,3,5];
        % location code 2
        right_cIFS = [11,13,15,12];
        % location code 3
        left_STG = [8,9,10,6];
        % location code 4
        right_STG = [16,19,20,18];
    end
    
    % normalized by controlled breathing task
    % controlled breathing max
    breath_blockAvg_left_cIFS = hmrBlockAvg(dc(:,:,right_cIFS),s_b,t,window_t);
    % normalized dc_task
    dc_task_left_cIFS = method_normalized_by_b(breath_blockAvg_left_cIFS, dc_task(:,:,right_cIFS));
    % select HbO only this time
    dc_task_left_cIFS_HbO = dc_task_left_cIFS(:,number);
    
    % find all task start index from s_task and save all block into one
    % matirix
    task_start_index = find(s_task == 1);
    
    for i = 1 : length(task_start_index)
        start_index = task_start_index + window_t(1)*Fs;
        end_index = task_start_index + window_t(2)*Fs;
        dummy = dc_task_left_cIFS_HbO(start_index : end_index);
        block_left_cIFS_HbO = cat(2, block_left_cIFS_HbO, dummy);
    end
    
    
    
    
end
%% passive data markerd as 0
% for exp1 data only take region cIFS
path = [pwd,'\data\passive_0'];
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);


%% plot and check
t = (window_t(1):1/50:window_t(2))';
% task trace
figure;
hold on;
all_average = mean(block_left_cIFS_HbO,2);
SE = std(block_left_cIFS_HbO,0,2);
plot(t,all_average,'r','LineWidth',2);
%calcualte SE and plot error bar
% ITD with err bar
err = SE/sqrt(size(block_left_cIFS_HbO,2));
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


