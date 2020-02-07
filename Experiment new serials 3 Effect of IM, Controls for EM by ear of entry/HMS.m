function [ HbO_f, HbR_f] = HMS(dc_HbO, dc_HbR, type)

k_s = 0.4;
k_f = -0.6;

switch type
    
    case 'sys'
        matrix = [k_s, -1; k_f * k_s, -k_s];
        
    case 'fun'
        matrix = [-k_s, 1; -k_s * k_f, k_f];
        
end


f_results = 1/(k_f - k_s) * matrix * [dc_HbO'; dc_HbR'];

HbO_f = f_results(1,:);
HbR_f = f_results(2,:);

end