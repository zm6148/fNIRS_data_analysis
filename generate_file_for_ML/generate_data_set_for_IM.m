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

data = [];
ROI = [];
SID = [];
label = [];

%% active data marked as 1
% for exp1 data only take region cIFS
path = [pwd,'\data\attentive_1'];
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);
ROI_list_all_subject = [];
blockAvg_all_subject = [];
subject_all_code = [];

for index=1:length(files)
    % load file and display name
    base_name = files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [path,'\',name,'.mat'];
    load(new_name);
    disp(base_name);
    splited = strsplit(base_name,'_');
    experiment_type = splited{2};
    
    % initiate things we want to save
    ROI_list = [];
    blockAvg_subject = [];
    subject_code = [];
    % task marks
    s_task = s_ITD_task + s_noITD_task;
    
    % which experiment
    if strcmp(experiment_type,'ITD')
        % no stg channel
        channels = [[20, 23, 21, 24];...% ROI code 1, left_cIFS
                    [12, 9, 11, 7]];    % ROI code 2, right_cIFS
    else
        channels = [[1,2,3,5]; ...    % ROI code 1, left_cIFS
                    [11,13,15,12];... % ROI code 2, right_cIFS
                    [8,9,10,6]; ...   % ROI code 3, left_STG
                    [16,19,20,18]];   % ROI code 4, right_STG
    end
    
    for ii = 1 : size(channels,1) 
        % disp(ii);
        % ROI code
        ROI_code = ii;
        % channels for this ROI
        channels_ROI = channels(ii,:);
        % normalized by controlled breathing task
        % controlled breathing max
        breath_blockAvg_ROI = hmrBlockAvg(dc(:,:,channels_ROI), s_b, t, window_t);
        % normalized dc_task
        dc_task_ROI = method_normalized_by_b(breath_blockAvg_ROI, dc_task(:,:,channels_ROI));
        % select HbO only this time
        dc_task_left_ROI = dc_task_ROI(:,number);
        
        % find all task start index from s_task and save all block into one
        % matirix
        task_start_index = find(s_task == 1);
        for jj = 1 : length(task_start_index)
            start_index = task_start_index(jj) + window_t(1)*Fs;
            end_index = task_start_index(jj) + window_t(2)*Fs;
            dummy = dc_task_left_ROI(start_index : end_index) - dc_task_left_ROI(start_index);
            blockAvg_subject = cat(2, blockAvg_subject, dummy);
            ROI_list = [ROI_list; ROI_code];
            subject_code = [subject_code; index];
        end
    end
    
    ROI_list_all_subject = cat(1, ROI_list_all_subject, ROI_list);
    subject_all_code = cat(1, subject_all_code, subject_code);
    blockAvg_all_subject = cat(1, blockAvg_all_subject, blockAvg_subject'); 
end

%% add to data, ROI, SID, label for active condition
data = cat(1, data, blockAvg_all_subject);
ROI = cat(1, ROI, ROI_list_all_subject);
SID = cat(1, SID, subject_all_code);
% active condition
label_active = ones(size(blockAvg_all_subject, 1),1);
label = cat(1, label, label_active);

%% passive data markerd as 0
% for exp1 data only take region cIFS
path = [pwd,'\data\passive_0'];
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);

ROI_list_all_subject = [];
blockAvg_all_subject = [];
subject_all_code = [];

for index=1:length(files)
    % load file and display name
    base_name = files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [path,'\',name,'.mat'];
    load(new_name);
    disp(base_name);
    splited = strsplit(base_name,'_');
    experiment_type = splited{2};
    
    % initiate things we want to save
    ROI_list = [];
    blockAvg_subject = [];
    subject_code = [];
    % task marks
    s_task = s_ITD_task + s_noITD_task;
    
    % which experiment
    if strcmp(experiment_type,'SSN')
        % no stg channel
        channels = [[20, 23, 21, 24];...% ROI code 1, left_cIFS
                    [12, 9, 11, 7]];    % ROI code 2, right_cIFS
    else
        channels = [[1,2,3,5]; ...    % ROI code 1, left_cIFS
                    [11,13,15,12];... % ROI code 2, right_cIFS
                    [8,9,10,6]; ...   % ROI code 3, left_STG
                    [16,19,20,18]];   % ROI code 4, right_STG
    end
    
    for ii = 1 : size(channels,1) 
        % disp(ii);
        % ROI code
        ROI_code = ii;
        % channels for this ROI
        channels_ROI = channels(ii,:);
        % normalized by controlled breathing task
        % controlled breathing max
        breath_blockAvg_ROI = hmrBlockAvg(dc(:,:,channels_ROI), s_b, t, window_t);
        % normalized dc_task
        dc_task_ROI = method_normalized_by_b(breath_blockAvg_ROI, dc_task(:,:,channels_ROI));
        % select HbO only this time
        dc_task_left_ROI = dc_task_ROI(:,number);
        
        % find all task start index from s_task and save all block into one
        % matirix
        task_start_index = find(s_task == 1);
        for jj = 1 : length(task_start_index)
            start_index = task_start_index(jj) + window_t(1)*Fs;
            end_index = task_start_index(jj) + window_t(2)*Fs;
            dummy = dc_task_left_ROI(start_index : end_index) - dc_task_left_ROI(start_index);
            blockAvg_subject = cat(2, blockAvg_subject, dummy);
            ROI_list = [ROI_list; ROI_code];
            subject_code = [subject_code; index];
        end
    end
    
    ROI_list_all_subject = cat(1, ROI_list_all_subject, ROI_list);
    subject_all_code = cat(1, subject_all_code, subject_code);
    blockAvg_all_subject = cat(1, blockAvg_all_subject, blockAvg_subject'); 
end

%% add to data, ROI, SID, label for passive condition
data = cat(1, data, blockAvg_all_subject);
ROI = cat(1, ROI, ROI_list_all_subject);
SID = cat(1, SID, subject_all_code);
% active condition
label_active = zeros(size(blockAvg_all_subject, 1),1);
label = cat(1, label, label_active);

% save data, ROI, SID, label to .mat file
save('ML_dataset.mat','data','ROI', 'SID', 'label');

%% plot and check
window_t=[-5,40];
t = (window_t(1):1/50:window_t(2))';
% see each region
index = find(label == 0);
data_to_plot = data(index, :)';
% task trace
figure(1);
hold on;
all_average = mean(data_to_plot,2);
SE = std(data_to_plot,0,2);
plot(t,all_average,'r','LineWidth',2);
%calcualte SE and plot error bar
% ITD with err bar
err = SE/sqrt(size(data_to_plot,2));
h = fill([t;flipud(t)],[all_average-err;flipud(all_average+err)],[1 0 0],'linestyle','none');
set(h,'facealpha',.3);



