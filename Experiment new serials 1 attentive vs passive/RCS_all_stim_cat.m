function [ b, r, p ] = RCS_all_stim_cat(dc_task_both, s_task_speech_both, s_task_SSN_both, RC, channels, breath_blockAvg_speech, breath_blockAvg_SSN,Hb_type)
% refrence channel as regressor
% GLM analysis on on HbO only
% dc_task(:,1,RC(ii)) middle number is for 1:HbO, 2: HbR, 3: HbT

% % find breath holding value at 20 second mark for speech and SSN
% speech_avg = mean(breath_blockAvg_speech(:,1,channels),3);
% SSN_avg = mean(breath_blockAvg_SSN(:,1,channels),3);
%
% % value at 20s mark
% speech_max = speech_avg(20*50);
% SSN_max = SSN_avg(20*50);
%
% both_max = (speech_max + SSN_max)/2;
%
% if both_max < 0
%     sign = 1;
% else
%     sign = -1;
% end

type = Hb_type; %HbR
switch type
    case 'HbO'
        number = 1;
    case 'HbR'
        number = 2;
    case 'HbT'
        number = 3;
end

% creat HRF function
[target_hrf_speech, target_hrf_SSN] = target_HRF_cat(s_task_speech_both, s_task_SSN_both, type);

first_deriv_speech = [0; diff(s_task_speech_both)];
first_deriv_SSN = [0; diff(s_task_SSN_both)];

% %dc_task;
% RC = 18;
% channels = [13, 16, 15, 17];

% RC channel
dc_HbO_RC_sum = 0;
for ii = 1:length(RC)
    
    % calcualte the max value of breathholding block average
    max_b_speech = max(abs(breath_blockAvg_speech(:,number,RC(ii))));
    max_b_SSN = max(abs(breath_blockAvg_SSN(:,number,RC(ii))));
    
    max_b = (max_b_speech + max_b_SSN) / 2;
    
    dc_HbO_RC_sum = dc_HbO_RC_sum + dc_task_both(:,number,RC(ii))/max_b;
    
end
dc_HbO_RC = dc_HbO_RC_sum/length(RC); %/(sign);

% HbO data before RCS
% dc_HbO_task_channels=dc_task(:,1,channels);
% each channel normalized by its max
dc_HbO_task_channels_sum = 0;

for ii = 1:length(channels)
    
    % calcualte the max value of breathholding block average
    max_b_speech = max(abs(breath_blockAvg_speech(:,number,channels(ii))));
    max_b_SSN = max(abs(breath_blockAvg_SSN(:,number,channels(ii))));
    
    max_b = (max_b_speech + max_b_SSN) / 2;
    
    dc_HbO_task_channels_sum = dc_HbO_task_channels_sum + dc_task_both(:,number,channels(ii))/max_b;
end

dc_HbO_task_channels_mean = dc_HbO_task_channels_sum/length(channels);


% build matrix for glm
% b = [ITD, noITD, combined, RC]
%                       x1 b(2)           x2 b(3)           x3 b(4)
% glm_regressors = cat(2, target_hrf_speech, target_hrf_SSN,  dc_HbO_RC);

if strcmp(type,'HbR')
    glm_regressors = cat(2, target_hrf_speech, target_hrf_SSN, first_deriv_speech, first_deriv_SSN); %, dc_HbO_RC);
else
    glm_regressors = cat(2, target_hrf_speech, target_hrf_SSN, first_deriv_speech, first_deriv_SSN); % dc_HbO_RC);
end

%           regressor       observed
mdl = fitlm(glm_regressors, dc_HbO_task_channels_mean);

% output results betas and r value
b = mdl.Coefficients.Estimate; %beta values
r = mdl.Rsquared.Ordinary; %values
p = mdl.Coefficients.pValue; % p value
end

