clear
load bla
close all
RT = 0.02;
T = 45;
data_length = length(s_ITD_task);

onset_time = 1;
for onset_time = 5
    %     p = [6 16 1 1 6 onset_time 30];
    %
    %     [h_ITD,p] = spm_hrf(RT,p,T);
    %     Fs = 50;
    %     h_ITD = [h_ITD; zeros(round(45*Fs)-length(h_ITD)+1,1)];
    %     ht = 0 : 1/Fs : (length(h_ITD)-1)/Fs;
    %
    %
    %
    %     %h_ITD = h_ITD./max(abs(h_ITD(:)));
    %     figure(2);plot(ht,h_ITD)
    %
    %     one_block = 0.*ht(ht<50);
    %     one_block(5<ht&ht<20) = 1;
    %     %t_hrt = real(ifft(fft(h_ITD).*fft(one_block')));
    %     t_hrt = conv(h_ITD,one_block');
    %     t_hrt = t_hrt(1:length(one_block));
    %     figure(3);hold on;plot(t_hrt)
    
    RT = 0.02;
    T = 45;
    
    p = [6 35 1 1 6 0 45];
    [h_ITD,p] = spm_hrf(RT,p,T);
    
    p = [6 35 1 1 6 0 45];
    [h_noITD,p] = spm_hrf(RT,p,T);
    
    
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
    
    
    target_hrf_ITD = conv(s_ITD_task_square,h_ITD,'same');
    %target_hrf_ITD = target_hrf_ITD(1:data_length)/max(target_hrf_ITD);
    
    t = 0: 1/50: (length(target_hrf_ITD)-1)/50;
    figure(1)
    hold on;
    
    plot(t, target_hrf_ITD)
end