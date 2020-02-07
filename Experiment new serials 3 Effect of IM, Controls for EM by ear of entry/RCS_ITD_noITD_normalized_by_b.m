function [ b, r, p ] = RCS_ITD_noITD_normalized_by_b(dc, s_b, dc_task, t, window_b, s_ITD_task, s_noITD_task, RC, channels, Hb_type)

% refrence channel as regressor
% GLM analysis on on HbO only
% dc_task(:,1,RC(ii)) middle number is for 1:HbO, 2: HbR, 3: HbT


type = Hb_type; %HbR
switch type
    case 'HbO'
        number = 1;
    case 'HbR'
        number = 2;
    case 'HbT'
        number = 3;
end


[target_hrf_ITD, target_hrf_noITD] = target_HRF(s_ITD_task,s_noITD_task, type);

first_deriv_ITD= [0; diff(target_hrf_ITD)];
first_deriv_noITD= [0; diff(target_hrf_noITD)];
% second_deriv = [0; diff(first_deriv)];
% third_deriv = [0; diff(second_deriv)];

% %dc_task;
% RC = 18;
% channels = [13, 16, 15, 17];

% each channel normalized by max_b
dc_HbO_RC_sum = 0;
for ii = 1:length(RC)
    % calcualte the max value of breathholding block average
    breath_blockAvg_subject = hmrBlockAvg(dc(:,number,RC(ii)), s_b, t, window_b);
    max_b = max(abs(breath_blockAvg_subject));
    dc_HbO_RC_sum = dc_HbO_RC_sum + dc_task(:,number,RC(ii))/max_b;
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
    glm_regressors = cat(2, target_hrf_ITD, target_hrf_noITD, first_deriv_ITD, first_deriv_noITD);%, dc_HbO_RC); %[0;diff(target_hrf_ITD)], [0;diff(target_hrf_noITD)]);
else
    glm_regressors = cat(2, target_hrf_ITD, target_hrf_noITD, first_deriv_ITD, first_deriv_noITD);%, dc_HbO_RC); %[0;diff(target_hrf_ITD)], [0;diff(target_hrf_noITD)]);
end

%                                        regressor       observed
mdl = fitlm(glm_regressors, dc_HbO_task_channels_mean');

% output results betas and r value
b = mdl.Coefficients.Estimate; %beta values
r = mdl.Rsquared.Ordinary; %values
p = mdl.Coefficients.pValue; % p value
end

