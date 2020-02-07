function [ HbO_f, HbR_f] = HMS_blvg(dc_HbO, dc_HbR)

k_s = 0.4;
k_f = -0.6;

matrix = [-k_s, 1; -k_s * k_f, k_f];

HbO_f = [];
HbR_f = [];
for subject = 1:size(dc_HbO,2)
    f_results = 1/(k_f - k_s) * matrix * [dc_HbO(:,subject)'; dc_HbR(:,subject)'];
    HbO_f = [HbO_f; f_results(1,:)];
    HbR_f = [HbR_f; f_results(2,:)];
    
end

HbO_f = HbO_f';
HbR_f = HbR_f';
end