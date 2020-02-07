function mean_diff = method_normalized_by_b(breath_blockAvg_subject, task_blockAvg_subject_ITD)

% normalize channel data based on peak value during breath holding task

diff = task_blockAvg_subject_ITD;
% max_breath = max(abs(breath_blockAvg_subject));

channels = size(breath_blockAvg_subject,3);
hb_type = size(breath_blockAvg_subject,2);

for ii = 1 : channels
    
    for jj = 1 : hb_type
        
        breath_data = breath_blockAvg_subject(:,jj,ii);
        max_breath = max(abs(breath_data));
        
        max_breath_type_channel = max_breath; 
        diff(:,jj, ii) = diff(:,jj, ii)/max_breath_type_channel;
        
    end
end

mean_diff = mean(diff,3);

end