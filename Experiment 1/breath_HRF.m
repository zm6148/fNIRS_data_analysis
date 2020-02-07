function breath_hrf = breath_HRF(s_b, hrf_case)

% HRF function of breath holding task
% GLM analysis on on HbO only 
% dc_task(:,1,RC(ii)) middle number is for 1:HbO, 2: HbR, 3: HbT

data_length = length(s_b);

% switch which HRF to model breath holding response
switch hrf_case
    case 'HbO'
        RT = 0.02;
        p = [5 25 1 1 6 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);
              
    case 'HbR'
        RT = 0.02;
        p = [25 45 1 1 6 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);

    case 'HbT'
        RT = 0.02;
        p = [5 25 1 1 6 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);

        
    case 3
        t=0:1/50:50;
        h = 0.6.*t.^2.1.*exp(-t/1.6)-0.0023.*t.^3.54.*exp(-t/4.25);
        
end


tick = find(s_b==1);
s_b_square=s_b;

for idx = 1:numel(tick)
    element = tick(idx);
    s_b_square(element:element+15*50)=1;
end

breath_hrf = conv(h,s_b_square);
breath_hrf = breath_hrf(1:data_length)/max(abs(breath_hrf));
%breath_hrf = flip(breath_hrf);
end

