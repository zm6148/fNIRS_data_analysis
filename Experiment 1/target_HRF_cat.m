function [target_hrf_speech,target_hrf_SSN] = target_HRF_cat( s_speech_task,s_SSN_task, type)

%function to build HRF for task 

data_length = length(s_speech_task);

switch type
    case 'HbO'
        RT = 0.02;
        p = [15 25 1 1 100 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);

        
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

        
    case 'HbR'
        RT = 0.02;
        p = [15 35 1 1 100 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);

%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

    case 'HbT'
        RT = 0.02;
        p = [3.33 25 1 1 80 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);

        
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

end

tick = find(s_speech_task==1);
s_ITD_task_square=s_speech_task;
for idx = 1:numel(tick)
    element = tick(idx);
    %disp(element);
    s_ITD_task_square(element:element+15*50)=1;
    
end

tick = find(s_SSN_task==1);
s_noITD_task_square=s_SSN_task;
for idx = 1:numel(tick)
    element = tick(idx);
    s_noITD_task_square(element:element+15*50)=1;
end


target_hrf_speech = conv(h,s_ITD_task_square);
target_hrf_speech=target_hrf_speech(1:data_length)/max(abs(target_hrf_speech));

target_hrf_SSN = conv(h,s_noITD_task_square);
target_hrf_SSN=target_hrf_SSN(1:data_length)/max(abs(target_hrf_SSN));


end

