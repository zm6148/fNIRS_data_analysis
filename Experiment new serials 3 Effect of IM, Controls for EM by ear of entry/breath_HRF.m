function breath_hrf = breath_HRF(s_b, hrf_case)
% HRF function for breath task

data_length = length(s_b);

% switch which HRF to model breath holding response
% gamma is used
switch hrf_case
    case 'HbO'
        
        RT = 0.02;
        p = [6 16 1 1 6 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);
        
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);
%         h = h/max(h);       
    case 'HbR'
        
        RT = 0.02;
        p = [16 26 1 1 6 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);
        
%         t=0:1/50:50;
%         h = 0.6.*t.^2.1.*exp(-t/1.6)-0.0023.*t.^3.54.*exp(-t/4.25);
%         h = h/max(h);     

%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);
%         h = h/max(h);       

    case 'HbT'
        
        RT = 0.02;
        p = [6 16 1 1 6 0 45];
        T = 45;
        [h,p] = spm_hrf(RT,p,T);
        
%         t=0:1/50:50;
%         h = (1/gamma(6)).*t.^5.*exp(-t)-(1/(6*gamma(16))).*t.^15.*exp(-t);
%         h = h/max(h);       

end


tick = find(s_b==1);
s_b_square=s_b;

for idx = 1:numel(tick)
    element = tick(idx);
    s_b_square(element:element+15*50)=1;
end

breath_hrf = conv(h,s_b_square);
breath_hrf=breath_hrf(1:data_length)/max(abs(breath_hrf));

end

