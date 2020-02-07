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
    
    
    
    E = t_HRF_allresults{ii}{1};
    A = t_HRF_allresults{ii}{2};
    F = t_HRF_allresults{ii}{3};
    B = t_HRF_allresults{ii}{4};
    
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


figure(1)
hold on;
plot(beta_onset_time_E_ITD,'r');
plot(beta_onset_time_A_ITD,'g');
plot(beta_onset_time_F_ITD,'b');
plot(beta_onset_time_B_ITD,'k');

plot(beta_onset_time_E_noITD,'r--');
plot(beta_onset_time_A_noITD,'g--');
plot(beta_onset_time_F_noITD,'b--');
plot(beta_onset_time_B_noITD,'k--');

figure(2)
hold on;
plot(r_all_onset_time_E,'r');
plot(r_all_onset_time_A,'g');
plot(r_all_onset_time_F,'b');
plot(r_all_onset_time_B,'k');