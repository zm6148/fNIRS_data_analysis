function [ b, r, p ] = RCS_breath_catselection( dc, s_b, RC, channels, Hb_type)
%refrence channel as regressor

% before anything trim the data
% 5 secs before the first stim mark and 45 secs after the last stim mark
% first_index = min(find(s_b==1));
% last_index = max(find(s_b==1));
% 
% start_index = first_index-50*5;
% end_index = last_index+50*45;
% 
% dc = dc(start_index:end_index,:,:);
% s_b = s_b(start_index:end_index);

type = Hb_type; %HbR
switch type
    case 'HbO'
        number = 1;
    case 'HbR'
        number = 2;
    case 'HbT'
        number = 3;
end

% HRF model
% 1: gamma function
% 2: difference between 2 gaussiz
breath_hrf = breath_HRF(s_b, Hb_type);

% each channel normalized by its variance
dc_HbO_RC_sum = 0;
for ii = 1:length(RC)
    dc_HbO_RC_sum = dc_HbO_RC_sum + dc(:,number,RC(ii));%/std(dc(:,number,RC(ii)));
end
dc_HbO_RC = dc_HbO_RC_sum/length(RC);

% HbO data before RCS
dc_HbO_breath_channels_sum = 0;
for ii = 1:length(channels)
    dc_HbO_breath_channels_sum = dc_HbO_breath_channels_sum + dc(:,number,channels(ii)); %/std(dc(:,number,channels(ii)));
end
dc_HbO_breath_channels_mean = dc_HbO_breath_channels_sum/length(channels);

% build matrix for glm
% b = [ITD, noITD, combined, RC]
%                       x1 b(2)         x2 b(3)           x3 b(4)    
glm_regressors = cat(2, breath_hrf, dc_HbO_RC);
%           regressor       observed
mdl = fitlm(glm_regressors, dc_HbO_breath_channels_mean);

% output results betas and r value 
b = mdl.Coefficients.Estimate; %beta values
r = mdl.Rsquared.Ordinary; %values
p = mdl.Coefficients.pValue; % p value
end

