function [b_E, r_E, p_E,...
          b_A, r_A, p_A,...
          b_F, r_F, p_F,...
          b_B, r_B, p_B] = select_best_HRF(t_HRF_allresults)


r_all_onset_time_E = [];
r_all_onset_time_A = [];
r_all_onset_time_F = [];
r_all_onset_time_B = [];

beta_onset_time_E_ITD = [];
beta_onset_time_A_ITD = [];
beta_onset_time_F_ITD = [];
beta_onset_time_B_ITD = [];

beta_onset_time_E_noITD = [];
beta_onset_time_A_noITD = [];
beta_onset_time_F_noITD = [];
beta_onset_time_B_noITD = [];


for ii = 1:length(t_HRF_allresults)
    %     E = t_HRF_allresults{ii}{1};
    %     A = t_HRF_allresults{ii}{2};
    %     F = t_HRF_allresults{ii}{3};
    %     B = t_HRF_allresults{ii}{4};
    
    r_all_onset_time_E = [r_all_onset_time_E, t_HRF_allresults{ii}{1}{2}];
    r_all_onset_time_A = [r_all_onset_time_A, t_HRF_allresults{ii}{2}{2}];
    r_all_onset_time_F = [r_all_onset_time_F, t_HRF_allresults{ii}{3}{2}];
    r_all_onset_time_B = [r_all_onset_time_B, t_HRF_allresults{ii}{4}{2}];
    
    beta_onset_time_E_ITD = [beta_onset_time_E_ITD, t_HRF_allresults{ii}{1}{1}(2)];
    beta_onset_time_A_ITD = [beta_onset_time_A_ITD, t_HRF_allresults{ii}{2}{1}(2)];
    beta_onset_time_F_ITD = [beta_onset_time_F_ITD, t_HRF_allresults{ii}{3}{1}(2)];
    beta_onset_time_B_ITD = [beta_onset_time_B_ITD, t_HRF_allresults{ii}{4}{1}(2)];
    
    beta_onset_time_E_noITD = [beta_onset_time_E_noITD, t_HRF_allresults{ii}{1}{1}(3)];
    beta_onset_time_A_noITD = [beta_onset_time_A_noITD, t_HRF_allresults{ii}{2}{1}(3)];
    beta_onset_time_F_noITD = [beta_onset_time_F_noITD, t_HRF_allresults{ii}{3}{1}(3)];
    beta_onset_time_B_noITD = [beta_onset_time_B_noITD, t_HRF_allresults{ii}{4}{1}(3)];
    
end



% [M_E,I_E] = max(r_all_onset_time_E);
% [M_A,I_A] = max(r_all_onset_time_A);
% [M_F,I_F] = max(r_all_onset_time_F);
% [M_B,I_B] = max(r_all_onset_time_B);

% diff_r_E = [eps diff(r_all_onset_time_E)];
% diff_r_A = [eps diff(r_all_onset_time_A)];
% diff_r_F = [eps diff(r_all_onset_time_F)];
% diff_r_B = [eps diff(r_all_onset_time_B)];

[pks,I_E] = findpeaks(r_all_onset_time_E);
[pks,I_A] = findpeaks(r_all_onset_time_A);
[pks,I_F] = findpeaks(r_all_onset_time_F);
[pks,I_B] = findpeaks(r_all_onset_time_B);

% disp(I_E);
% disp(I_A);
% disp(I_F);
% disp(I_B);

b_E = t_HRF_allresults{I_E(1)}{1}{1};
r_E = t_HRF_allresults{I_E(1)}{1}{2};
p_E = t_HRF_allresults{I_E(1)}{1}{3};

b_A = t_HRF_allresults{I_A(1)}{2}{1};
r_A = t_HRF_allresults{I_A(1)}{2}{2};
p_A = t_HRF_allresults{I_A(1)}{2}{3};

b_F = t_HRF_allresults{I_F(1)}{3}{1};
r_F = t_HRF_allresults{I_F(1)}{3}{2};
p_F = t_HRF_allresults{I_F(1)}{3}{3};

b_B = t_HRF_allresults{I_B(1)}{4}{1};
r_B = t_HRF_allresults{I_B(1)}{4}{2};
p_B = t_HRF_allresults{I_B(1)}{4}{3};
end