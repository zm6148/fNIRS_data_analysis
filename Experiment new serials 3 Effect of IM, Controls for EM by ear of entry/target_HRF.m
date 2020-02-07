function [target_hrf_ITD,target_hrf_noITD] = target_HRF( s_ITD_task,s_noITD_task, type)

% task block train as regressor
% GLM analysis on on HbO only
% dc_task(:,1,RC(ii)) middle number is for 1:HbO, 2: HbR, 3: HbT


data_length = length(s_ITD_task);

% t=0:1/50:50;
% h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);
% h = h/max(h);

switch type
    case 'HbO'
        RT = 0.02;
        T = 45;
        
        p = [6 16 1 1 6 0 45];
        [h_ITD,p] = spm_hrf(RT,p,T);
        
        p = [6 16 1 1 6 0 45];
        [h_noITD,p] = spm_hrf(RT,p,T);
        
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);
            
    case 'HbR'
        RT = 0.02;
        T = 45;
        
        p = [16 16 1 1 6 0 45];
        [h_ITD,p] = spm_hrf(RT,p,T);
        
        p = [16 16 1 1 6 0 45];
        [h_noITD,p] = spm_hrf(RT,p,T);
        

        %                 t=0:1/50:50;
        %                 h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);
        
    case 'HbT'
        RT = 0.02;
        T = 45;
        
        p = [6 16 1 1 6 0 45];
        [h_ITD,p] = spm_hrf(RT,p,T);
        
        p = [6 16 1 1 6 0 45];
        [h_noITD,p] = spm_hrf(RT,p,T);
        
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

end         


tick = find(s_ITD_task==1);
s_ITD_task_square=s_ITD_task;
for idx = 1:numel(tick)
    element = tick(idx);
    s_ITD_task_square(element:element+50*15)=1;
end

tick = find(s_noITD_task==1);
s_noITD_task_square=s_noITD_task;
for idx = 1:numel(tick)
    element = tick(idx);
    s_noITD_task_square(element:element+50*15)=1;
end


target_hrf_ITD = conv(h_ITD,s_ITD_task_square);
target_hrf_ITD = target_hrf_ITD(1:data_length)/max(target_hrf_ITD);

target_hrf_noITD = conv(h_noITD,s_noITD_task_square);
target_hrf_noITD = target_hrf_noITD(1:data_length)/max(target_hrf_noITD);


end

