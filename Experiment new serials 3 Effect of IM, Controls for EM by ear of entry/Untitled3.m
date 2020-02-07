close all
%load('C:\Users\mz86\Desktop\Fnirs data and analysis\Experiment new serials 3 Effect of IM, Controls for EM by ear of entry\TOH_2_20191022_1227_06_processed_data.mat');
% channel
RC = 14;
channels = [11,13,15,12];

% Difine window
window_b=[-5,40];
window_t=[-5,40];

type = 'HbO'; %HbR
switch type
    case 'HbO'
        number = 1;
        onset_time = 10;
    case 'HbR'
        number = 2;
        onset_time = 15;
    case 'HbT'
        number = 3;
end



[target_hrf_ITD, target_hrf_noITD] = target_HRF_onset_time(s_ITD_task,s_noITD_task, type, onset_time);

% each channel normalized by max_b
dc_HbO_RC_sum = 0;
for ii = 1:length(RC)
    % calcualte the max value of breathholding block average
    breath_blockAvg_subject = hmrBlockAvg(dc(:,number,RC), s_b, t, window_b);
    max_b = max(abs(breath_blockAvg_subject));
    dc_HbO_RC_sum = dc_HbO_RC_sum + dc_task(:,number,RC)/max_b;
end
dc_HbO_RC = dc_HbO_RC_sum/length(RC);

% HbO data before RCS
% dc_HbO_task_channels=dc_task(:,1,channels);
% each channel normalized by max_b
dc_HbO_task_channels_sum = 0;
for ii = 1:length(channels)
    % calcualte the max value of breathholding block average
    breath_blockAvg_subject = hmrBlockAvg(dc(:,number,channels(ii)), s_b, t, window_b);
    max_b = max(abs(breath_blockAvg_subject));
    dc_HbO_task_channels_sum = dc_HbO_task_channels_sum + dc_task(:,number,channels(ii))/max_b;
end
dc_HbO_task_channels_mean = dc_HbO_task_channels_sum/length(channels);


% build matrix for glm
% b = [ITD, noITD, combined, RC]
%                       x1 b(2)         x2 b(3)           x3 b(4)
%disp(size(target_hrf_ITD));
if strcmp(type,'HbR')
    glm_regressors = cat(2, target_hrf_ITD, target_hrf_noITD);%, dc_HbO_RC); %[0;diff(target_hrf_ITD)], [0;diff(target_hrf_noITD)]);
else
    glm_regressors = cat(2, target_hrf_ITD, target_hrf_noITD, dc_HbO_RC); %[0;diff(target_hrf_ITD)], [0;diff(target_hrf_noITD)]);
end

%                                        regressor       observed
mdl = fitlm(glm_regressors, dc_HbO_task_channels_mean');

% output results betas and r value
b = mdl.Coefficients.Estimate; %beta values
r = mdl.Rsquared.Ordinary; %values
p = mdl.Coefficients.pValue; % p value

beta_ITD = b(2);
beta_noITD = b(3);


blkavg_ITD = hmrBlockAvg(dc_HbO_task_channels_mean - dc_HbO_RC * b(1), s_ITD_task,t,window_t);
blkavg_noITD = hmrBlockAvg(dc_HbO_task_channels_mean - dc_HbO_RC * b(1), s_noITD_task,t,window_t);

blkavg_fit_ITD = hmrBlockAvg(target_hrf_ITD, s_ITD_task,t,window_t)* beta_ITD;
blkavg_fit_noITD = hmrBlockAvg(target_hrf_noITD, s_noITD_task,t,window_t)* beta_noITD;

%% plot
figure (1)
hold on;
t = (0:1/50:(length(blkavg_fit_ITD)-1)/50)-5;
plot(t, blkavg_ITD,'g');
plot(t, blkavg_fit_ITD,'g--');
plot(t, blkavg_noITD,'b');
plot(t, blkavg_fit_noITD,'b--');
%% plot
figure(2)
t = 0: 1/50: (length(target_hrf_ITD)-1)/50;
hold on;
plot(t, dc_HbO_RC * b(1), 'k');
plot(t, dc_HbO_task_channels_mean,...
    'r',...
    'LineWidth',2);
plot(t, target_hrf_ITD, 'g--');
plot(t, target_hrf_noITD, 'b--');

y_limit=get(gca,'ylim');
x_limit=get(gca,'xlim');
ones_index = find(s_ITD_task==1);
for ii = 1 : length(ones_index)
%     plot([ones_index/50, ones_index/50],y_limit,'r')
%     plot([ones_index/50 + 15, ones_index/50 + 15],y_limit,'r')
   
    % Add lines
    %h1 = line([ones_index/50, ones_index/50],y_limit);
    %h2 = line([ones_index/50 + 15, ones_index/50 + 15],y_limit);
    % Set properties of lines
    %set([h1 h2],'Color','k','LineWidth',2)
    % Add a patch
    h = patch([ones_index(ii)/50, ones_index(ii)/50 + 15, ones_index(ii)/50 + 15, ones_index(ii)/50],[y_limit(1) y_limit(1) y_limit(2) y_limit(2)],'k');
    set(h,'facealpha',.3);
end

ones_index = find(s_noITD_task==1);
for ii = 1 : length(ones_index)
    h = patch([ones_index(ii)/50, ones_index(ii)/50 + 15, ones_index(ii)/50 + 15, ones_index(ii)/50],[y_limit(1) y_limit(1) y_limit(2) y_limit(2)],'k');
    set(h,'facealpha',.3);
end
plot(t, s_ITD_task, 'k');
plot(t, s_noITD_task, 'k');



