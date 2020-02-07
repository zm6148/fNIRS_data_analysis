%% data analysisi on ITD and noITD condition seperatly
%% load data folder
clear;
close all;

path = uigetdir(pwd,'choose reuslts folder');
filepattern=fullfile(path,'/*.mat');
files=dir(filepattern);

% path = [pwd,'\Data\Processed\trace_selected'];
% filepattern=fullfile(path,'/*.mat');
% files=dir(filepattern);

splited = strsplit(path,'\');
experiment_type = splited{8};

%% define basic stuff
% Difine window
window_b=[-5,40];
window_t=[-5,40];

%% load individual data

%(C, A, D, B, CA, DB, CD, AB) * (17 subjects)

subject_beta = [];
all_glm_beta = [];
subject_r2 = [];
all_glm_r2 = [];
subject_p = [];
all_glm_p = [];

subject_beta_breath = [];
all_glm_beta_breath = [];
subject_r2_breath = [];
all_glm_r2_breath = [];
subject_p_breath = [];
all_glm_p_breath = [];

Hb_type = 'HbO';
sys_or_fun = 'fun';

type = Hb_type; %HbR
switch type
    case 'HbO'
        number = 1;
    case 'HbR'
        number = 2;
    case 'HbT'
        number = 3;
end

for index=1:length(files)
    base_name = files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [path,'\',name,'.mat'];
    load(new_name);
    disp(base_name);
    
    % implement HMS before doing any block average
    dc_b = dc;
    new_dc = [];
    for ii = 1:size(dc,3)
        [dc_HbO_f, dc_HbR_f] = HMS(dc(:,1, ii), dc(:,2, ii), sys_or_fun);
        temp = cat(2, dc_HbO_f', dc_HbR_f', dc_HbO_f'+ dc_HbR_f');
        new_dc = cat(3, new_dc, temp);
    end
    % assign it back to old dc
    %dc = new_dc;
    
    % build Refrence channel regressor
    % for source C RC: 18; channels: 13, 16, 15, 17
    % for source D RC: 19; channels: 20, 23, 21, 24
    % for source A RC: 4; channels: 2, 6, 3, 1
    % for source B RC: 10; channels: 12, 9, 11, 7
    % use regressor closest to the source
    
    RC_A = 4;
    channels_A = [1,2,3,5];
    
    RC_B = 7;
    channels_B = [8,9,10,6];%[5,7,8,9];
    
    RC_E = 14;
    channels_E = [11,13,15,12];
    
    RC_F = 17;
    channels_F = [16,19,20,18];
    
    % combine CA and DB
    %     RC_CA = [18,4];
    %     channels_CA = [13, 16, 15, 17, 2, 6, 3, 1];
    %
    %     RC_DB = [19, 10];
    %     channels_DB = [12, 9, 11, 7, 20, 23, 21, 24];
    
    % combine CD (left) and AB (right)
    %     RC_CD = [18, 19];
    %     channels_CD = [13, 16, 15, 17, 20, 23, 21, 24];
    %
    %     RC_AB = [4, 10];
    %     channels_AB = [2, 6, 3, 1, 12, 9, 11, 7];
    % GLM analysis on HbO of the breath holding data
    [ b_E_breath, r_E_breath, p_E_breath] = RCS_breath( dc_b, s_b, RC_E, channels_E, Hb_type);
    [ b_A_breath, r_A_breath, p_A_breath] = RCS_breath( dc_b, s_b, RC_A, channels_A, Hb_type);
    [ b_F_breath, r_F_breath, p_F_breath] = RCS_breath( dc_b, s_b, RC_F, channels_F, Hb_type);
    [ b_B_breath, r_B_breath, p_B_breath] = RCS_breath( dc_b, s_b, RC_B, channels_B, Hb_type);
    %     [ b_CA_breath, r_CA_breath, p_CA_breath] = RCS_breath( dc, s_b, RC_CA, channels_CA  );
    %     [ b_DB_breath, r_DB_breath, p_DB_breath] = RCS_breath( dc, s_b, RC_DB, channels_DB  );
    %     [ b_CD_breath, r_CD_breath, p_CD_breath] = RCS_breath( dc, s_b, RC_CD, channels_CD  );
    %     [ b_AB_breath, r_AB_breath, p_AB_breath] = RCS_breath( dc, s_b, RC_AB, channels_AB  );
    
    % for all on set time loop through all possible
    % 1:45
    glm_beta_onset_time = [];
    for onset_time = 0:5:25
        
        
        % GLM analysis on HbO of the ITD and no_ITD
        [ b_E, r_E, p_E] = RCS_ITD_noITD_normalized_by_b_onset_time( dc, s_b, dc_task, t, window_b, s_ITD_task, s_noITD_task, RC_E, channels_E, Hb_type, onset_time);
        [ b_A, r_A, p_A] = RCS_ITD_noITD_normalized_by_b_onset_time( dc, s_b, dc_task, t, window_b, s_ITD_task, s_noITD_task, RC_A, channels_A, Hb_type, onset_time);
        [ b_F, r_F, p_F] = RCS_ITD_noITD_normalized_by_b_onset_time( dc, s_b, dc_task, t, window_b, s_ITD_task, s_noITD_task, RC_F, channels_F, Hb_type, onset_time);
        [ b_B, r_B, p_B] = RCS_ITD_noITD_normalized_by_b_onset_time( dc, s_b, dc_task, t, window_b, s_ITD_task, s_noITD_task, RC_B, channels_B, Hb_type, onset_time);
        
        
        %     [ b_CA, r_CA, p_CA] = RCS( dc_task, s_ITD_task,s_noITD_task, RC_CA, channels_CA  );
        %     [ b_DB, r_DB, p_DB] = RCS( dc_task, s_ITD_task,s_noITD_task, RC_DB, channels_DB  );
        %     [ b_CD, r_CD, p_CD] = RCS( dc_task, s_ITD_task,s_noITD_task, RC_CD, channels_CD  );
        %     [ b_AB, r_AB, p_AB] = RCS( dc_task, s_ITD_task,s_noITD_task, RC_AB, channels_AB  );
        
        % save data
        % build a matrix that has all data for all subjects
        % ITD b(2); noITD b(3); RC b(4)
        subject_beta = cat(2, b_E, b_A, b_F, b_B);%, b_CA, b_DB, b_CD, b_AB);
        glm_beta_onset_time = cat(3, glm_beta_onset_time, subject_beta);
        subject_r2 = cat(1, r_E, r_A, r_F, r_B); %, r_CA, r_DB, r_CD, r_AB);
        all_glm_r2 = cat(2, all_glm_r2, subject_r2);
        subject_p = cat(2, p_E, p_A, p_F, p_B); %, p_CA, p_DB, p_CD, p_AB);
        all_glm_p = cat(3, all_glm_p, subject_p);
        
        
        
        
    end
    all_glm_beta = cat(4, all_glm_beta, glm_beta_onset_time);
    
    % save data breath holding datas
    % build a matrix that has all data for all subjects
    % ITD b(2); noITD b(3); RC b(4)
    subject_beta_breath = cat(2, b_E_breath, b_A_breath, b_F_breath, b_B_breath); %, b_CA_breath, b_DB_breath, b_CD_breath, b_AB_breath);
    all_glm_beta_breath = cat(3, all_glm_beta_breath, subject_beta_breath);
    subject_r2_breath = cat(1, r_E_breath, r_A_breath, r_F_breath, r_B_breath); %, r_CA_breath, r_DB_breath, r_CD_breath, r_AB_breath);
    all_glm_r2_breath = cat(2, all_glm_r2_breath, subject_r2_breath);
    subject_p_breath = cat(2, p_E_breath, p_A_breath, p_F_breath, p_B_breath);%, p_CA_breath, p_DB_breath, p_CD_breath, p_AB_breath);
    all_glm_p_breath = cat(3, all_glm_p_breath, subject_p_breath);
end

%% GLM data analysis

% check breath holding GLM model and only select the ones with good GLM
% model to proceed with ITD noITD analysis

% all subject with breath holding r2 bigger than 0.5
% all_good_fit_subject_index 8 * 17 cell
% 8: C, A, D, B, CA, DB, CD, AB
% 17: 17 subjects in total

all_good_fit_subject_index = cell(size(all_glm_r2_breath,1),1);
for ii = 1:size(all_glm_r2_breath,1)
    % r2 value for breath
    good_r2_index = find(all_glm_r2_breath(ii,:) > 0);% 0.5
    % p vlaue for breath
    good_pvalue_index = find(all_glm_p_breath(2,ii,:) < 50);% 0.05
    % beta value for breath
    good_bets_value_index = find(all_glm_beta_breath(2,ii,:) < 50);
    % select only good p and good r2
    good_fit_index_1 = intersect(good_r2_index, good_pvalue_index);
    % select only beta smaller than 0
    good_fit_index = intersect(good_bets_value_index, good_fit_index_1);
    
    all_good_fit_subject_index{ii} = good_fit_index;
    
    
    n = length(good_fit_index);
    
    subject_r2_breath = all_glm_r2_breath(ii,good_fit_index);
    
    subject_pvalue_breath = all_glm_p_breath(2,ii,good_fit_index);
    subject_pvalue_breath = reshape(subject_pvalue_breath,[1,n]);
    
    subject_beta_breath = all_glm_beta_breath(2,ii,good_fit_index);
    subject_beta_breath = reshape(subject_beta_breath,[1,n]);
    
    
    subject_breath_table{ii} = cat(2,good_fit_index, subject_beta_breath', subject_pvalue_breath', subject_r2_breath');
    
end

%

% at each region select only good subjects
% table row : E A F B
% table column : subject_index, subject_beta_ITD', subject_beta_noITD', subject_pvalue_ITD', subject_pvalue_noITD', subject_r2'
for ii = 1:length(all_good_fit_subject_index)
    
    subject_index = all_good_fit_subject_index{ii};
    
    subject_r2 = all_glm_r2(ii,subject_index);
    
    n = length(subject_index);
    
    subject_pvalue_ITD = all_glm_p(2,ii,subject_index);
    subject_pvalue_ITD = reshape(subject_pvalue_ITD,[1,n]);
    subject_pvalue_noITD = all_glm_p(3,ii,subject_index);
    subject_pvalue_noITD = reshape(subject_pvalue_noITD,[1,n]);
    
    subject_beta_ITD = all_glm_beta(2,ii,subject_index);
    subject_beta_ITD = reshape(subject_beta_ITD,[1,n]);
    subject_beta_noITD = all_glm_beta(3,ii,subject_index);
    subject_beta_noITD = reshape(subject_beta_noITD,[1,n]);
    
    subject_selection_by_breath_table{ii} = cat(2,subject_index, subject_beta_ITD', subject_beta_noITD', subject_pvalue_ITD', subject_pvalue_noITD', subject_r2');
end

% %% further the selection based on task portion model
%
% for ii = 1:length(all_good_fit_subject_index)
%
%     results_region = subject_selection_by_breath_table{ii};
%     all_index = 1:(size(results_region,1));
%
%     % exclude p value not smaller than 0.05
%     % for ITD and no_ITD both
%     r2_exclude = find(results_region(:,6) < 0.5); % exclude r2
%     p_ITD_exclude = find(results_region(:,4) > 0.05); % exclude p ITD
%     p_noITD_exclude = find(results_region(:,5) > 0.05); % exclude p noITD
%
%     exclude_index = union(p_ITD_exclude, p_noITD_exclude);
%     exclude_index = union(exclude_index, r2_exclude);
%     all_index(exclude_index)=[];
%
%
%     subject_selection_by_breath_table_further{ii} = results_region(all_index,:);
% end

%% save the selected subject code as in a file

% region E
channels_at_E = subject_selection_by_breath_table{1}(:,1);
% region A
channels_at_A = subject_selection_by_breath_table{2}(:,1);
% region F
channels_at_F = subject_selection_by_breath_table{3}(:,1);
% region B
channels_at_B = subject_selection_by_breath_table{4}(:,1);

selected_channels = {channels_at_E; channels_at_A; channels_at_F; channels_at_B};

% % if ITD or SSN saved in different name
% if strcmp(experiment_type,'ITD')
%     subject_data_name = [path, '/selected channels/', 'selected_channels.mat'];
%     save(subject_data_name,'selected_channels');
% else
%     subject_data_name = [path, '/selected channels/', 'selected_channels.mat'];
%     save(subject_data_name,'selected_channels');
% end



%% plot and see breath portion
% different region may have different subjects
% 2 cases:1, only select by breath. 2, further selection by ITD noITD
figure(1);
title('beta values at recording site');
names = {'E (right STG) ';...
    'A (left STG)';...
    'F (right cIFS/tgPCS)';...
    'B (left cIFS/tgPCS)'
    };

breath_plot_data = [];
for ii = 1:length(all_good_fit_subject_index)
    
    subplot(2,2,ii);
    
    results_region = subject_breath_table{ii};
    
    
    b_beta = results_region(:,2);
    
    % t-test
    [h,p,ci,stats] = ttest(b_beta);
    
    % signed rank test
    %[p,h] = signrank(b_beta);
    
    beta_mean = [mean(b_beta)];
    beta_SEM = [std(b_beta)/sqrt(length(b_beta))];
    
    breath_plot_data = [breath_plot_data; [beta_mean, beta_SEM, p]];
    
    hold on;
    bar(2, beta_mean(1),'r');
    errorbar(2, beta_mean,beta_SEM,'.');
    set(gca,'XTick', [0 1 2]);
    set(gca,'XTickLabel',{'', 'Breath Holding', ''})
    title(['breath ', names{ii}, ';p = ', num2str(p), ';tstats=', num2str(stats.tstat),';df=',num2str(stats.df)]);
    ylim([-6*10^(-7) 6*10^(-7)]);
end



%% plot and see task portion
% different region may have different subjects
% 2 cases:1, only select by breath. 2, further selection by ITD noITD
which = 1;
%figure;
%title('beta values at recording site');
names = {'E (right STG) ';...
    'A (left STG)';...
    'F (right cIFS/tgPCS)';...
    'B (left cIFS/tgPCS)'
    };


task_ITD_plot_data = [];
task_noITD_plot_data = [];
for ii = 1:length(all_good_fit_subject_index)
    
    %subplot(2,2,ii);
    switch which
        case 1
            results_region = subject_selection_by_breath_table{ii};
        case 2
            results_region = subject_selection_by_breath_table_further{ii};
    end
    
    ITD_beta = results_region(:,2);
    noITD_beta = results_region(:,3);
    
    % signed rank test agansit meidan 0
    %     [p_ITD_beta,h_ITD_beta,rank_ITD] = signrank(ITD_beta);
    %     [p_noITD_beta,h_noITD_beta,rank_noITD] = signrank(noITD_beta);
    %     [p_c_s,h_c_s,rank_both] = signrank(ITD_beta,noITD_beta);
    
    [h_ITD, p_ITD, ci_ITD, stats_ITD] = ttest(ITD_beta);
    [h_noITD, p_noITD, ci_noITD, stats_noITD] = ttest(noITD_beta);
    [h_c, p_c, ci_c, stats_c] = ttest(ITD_beta,noITD_beta);
    
    % t-test
    %     [h_ITD_beta,p_ITD_beta] = ttest(ITD_beta);
    %     [h_noITD_beta,p_noITD_beta] = ttest(noITD_beta);
    %     [h_c_s,p_c_s] = ttest(ITD_beta,noITD_beta);
    
    % t-test
    %[h,p,ci,stats] = ttest2(ITD_beta,noITD_beta);
    
    beta_mean = [mean(ITD_beta), mean(noITD_beta)];
    beta_SEM = [std(ITD_beta)/sqrt(length(ITD_beta)), std(noITD_beta)/sqrt(length(noITD_beta))];
    
    % add adjustement for rANOVA
    diff_from_mean_ITD = kron(ones(size(ITD_beta,1),1),mean(ITD_beta,1)) - kron(ones(size(ITD_beta)),mean(mean(ITD_beta,1),2));
    adjusted_data_ITD = ITD_beta - diff_from_mean_ITD;
    adjusted_std_error_ITD = std(adjusted_data_ITD)./sqrt(length(ITD_beta));
    
    diff_from_mean_noITD = kron(ones(size(noITD_beta,1),1),mean(noITD_beta,1)) - kron(ones(size(noITD_beta)),mean(mean(noITD_beta,1),2));
    adjusted_data_noITD = noITD_beta - diff_from_mean_noITD;
    adjusted_std_error_noITD = std(adjusted_data_noITD)./sqrt(length(noITD_beta));
    % end of adjustment
    
    % save for plot later
    ITD_data = [mean(ITD_beta), std(ITD_beta)/sqrt(length(ITD_beta)), p_ITD ];
    noITD_data = [mean(noITD_beta), std(noITD_beta)/sqrt(length(noITD_beta)), p_noITD];
    
    %     ITD_data = [mean(adjusted_data_ITD),  adjusted_std_error_ITD, p_ITD ];
    %     noITD_data = [mean(adjusted_data_noITD), adjusted_std_error_noITD, p_noITD];
    
    task_ITD_plot_data = [task_ITD_plot_data; ITD_data];
    task_noITD_plot_data = [task_noITD_plot_data; noITD_data];
    
    %     hold on
    %     bar(1,beta_mean(1),'r');
    %     bar(2,beta_mean(2),'y');
    %     errorbar(1:2,beta_mean,beta_SEM,'.');
    %     set(gca,'XTick', [0 1 2 3]);
    %     set(gca,'XTickLabel',{'', 'ITD', 'no ITD', ''})
    %     title(['Task ', names{ii}, ' p=', num2str(p_c), ' p\_ITD =', num2str(p_ITD), ' p\_noITD=',num2str(p_noITD)]);
    
end

% save the subject plot
% set(gcf, 'PaperUnits', 'inches');
% x_width=8 ;y_width=8;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% switch which
%     case 1
%         saveas(gcf, [path, '/GLM/', 'beta_selection_by_breath_ITD_noITD'], 'pdf') %Save figure
%     case 2
%         saveas(gcf, [path, '/GLM/', 'beta_selection_by_breath_further'], 'pdf') %Save figure
% end

%% plot all 4 region in one bar plot
figure(2);
hold on;
x = [1:0.5:2.5]';
y = cat(2, task_ITD_plot_data(:,1), task_noITD_plot_data(:,1));
bar(x, y, 'BarWidth',0.9);
%ylim([-0.05 0.45]);

% for i1=1:4
%     x_location = [1-0.075, 1.5-0.075, 2-0.075, 2.5-0.075];
%     text(x(i1), y(i1), num2str(task_ITD_plot_data(i1,3),'%10.3e'),...
%                'HorizontalAlignment','center',...
%                'VerticalAlignment','bottom');
% end
%
% for i1=1:4
%     x_location = [1+0.075, 1.5+0.075, 2+0.075, 2.5+0.075];
%     text(x(i1), y(i1), num2str(task_noITD_plot_data(i1,3),'%10.3e'),...
%                'HorizontalAlignment','center',...
%                'VerticalAlignment','bottom');
% end

errorbar([1-0.075, 1.5-0.075, 2-0.075, 2.5-0.075], task_ITD_plot_data(:,1),task_ITD_plot_data(:,2),'.');
errorbar([1+0.075, 1.5+0.075, 2+0.075, 2.5+0.075], task_noITD_plot_data(:,1),task_noITD_plot_data(:,2),'.');

set(gca,'XTick', [0.5 1 1.5 2 2.5 3]);
set(gca,'XTickLabel',{'', 'E (right cIFS)', 'A (left cIFS)', 'F (right STG)', 'B (left STG)'});
set(gca,'TickDir','out');
%ylim([-1*10^(-7) 3*10^(-7)]);
title('perceptual task beta value');
legend('ITD', 'Opposite Ear');

%% plot and see
% % different region have same number of  subjects
% % 2 cases:1, only select by breath. 2, further selection by ITD noITD
% which = 1;
% figure;
% title('beta values at recording site');
% names = {'C (left tgPCS) ';...
%     'A (right tgPCS)';...
%     'D (left cIFS)';...
%     'B (right cIFS)';...
%     'CA (tgPCS)';...
%     'DB (cIFS)';...
%     'CD (left)';...
%     'AB (right)'
%     };
%
%
% for ii = 1:length(all_good_fit_subject_index)
%
%     subplot(4,2,ii);
%
%     results_region = subject_selection_by_breath_table_further{ii};
%
%     index = [];
%     for jj = 1:length(common)
%         dummy = find(results_region(:,1)==common(jj));
%         index = [index, dummy];
%     end
%
%     results_region = results_region(index,:);
%
%     ITD_beta = results_region(:,2);
%     noITD_beta = results_region(:,3);
%
%     beta_mean = [mean(ITD_beta), mean(noITD_beta)];
%     beta_SEM = [std(ITD_beta)/sqrt(length(ITD_beta)), std(noITD_beta)/sqrt(length(noITD_beta))];
%
%     hold on
%     bar(1,beta_mean(1),'r');
%     bar(2,beta_mean(2),'y');
%     errorbar(1:2,beta_mean,beta_SEM,'.');
%     set(gca,'XTick', [0 1 2 3]);
%     set(gca,'XTickLabel',{'', 'ITD', 'no ITD', ''})
%     title(names{ii});
% end
%
% % save the subject plot
% set(gcf, 'PaperUnits', 'inches');
% x_width=8 ;y_width=8;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% switch which
%     case 1
%         saveas(gcf, [path, '/GLM/', 'beta_selection_by_breath_same'], 'pdf') %Save figure
%     case 2
%         saveas(gcf, [path, '/GLM/', 'beta_selection_by_breath_further_same'], 'pdf') %Save figure
% end

%%
% save all subject GLM data
subject_data_name = [pwd, '\GLM_data_exp2.mat'];
save(subject_data_name,'subject_selection_by_breath_table');

%% plot beta over time

all_glm_beta_E_ITD = [];
all_glm_beta_E_noITD = [];
all_glm_beta_A_ITD = [];
all_glm_beta_A_noITD = [];
all_glm_beta_F_ITD = [];
all_glm_beta_F_noITD = [];
all_glm_beta_B_ITD = [];
all_glm_beta_B_noITD = [];


for jj = 1:size(all_glm_beta,4)
    
    glm_beta_E_ITD = [];
    glm_beta_E_noITD = [];
    glm_beta_A_ITD = [];
    glm_beta_A_noITD = [];
    glm_beta_F_ITD = [];
    glm_beta_F_noITD = [];
    glm_beta_B_ITD = [];
    glm_beta_B_noITD = [];
    
    for ii = 1:size(all_glm_beta,3)
        
        glm_beta = all_glm_beta(:,:,ii,jj);
        
        glm_beta_E_ITD = [glm_beta_E_ITD, glm_beta(1,1)];
        glm_beta_E_noITD =  [glm_beta_E_noITD glm_beta(2,1)];
        glm_beta_A_ITD =  [glm_beta_A_ITD glm_beta(1,2)];
        glm_beta_A_noITD =  [glm_beta_A_noITD glm_beta(2,2)];
        glm_beta_F_ITD =  [glm_beta_F_ITD glm_beta(1,3)];
        glm_beta_F_noITD = [glm_beta_F_noITD glm_beta(2,3)];
        glm_beta_B_ITD = [glm_beta_B_ITD  glm_beta(1,4)];
        glm_beta_B_noITD = [glm_beta_B_noITD  glm_beta(2,4)];
    end
    
    all_glm_beta_E_ITD = [all_glm_beta_E_ITD; glm_beta_E_ITD];
    all_glm_beta_E_noITD = [all_glm_beta_E_noITD; glm_beta_E_noITD];
    all_glm_beta_A_ITD = [all_glm_beta_A_ITD; glm_beta_A_ITD];
    all_glm_beta_A_noITD = [all_glm_beta_A_noITD; glm_beta_A_noITD];
    all_glm_beta_F_ITD = [all_glm_beta_F_ITD; glm_beta_F_ITD];
    all_glm_beta_F_noITD = [all_glm_beta_F_noITD; glm_beta_F_noITD];
    all_glm_beta_B_ITD = [all_glm_beta_B_ITD; glm_beta_B_ITD];
    all_glm_beta_B_noITD = [all_glm_beta_B_noITD; glm_beta_B_noITD];
    
end

all_data = {all_glm_beta_E_ITD; all_glm_beta_E_noITD; ...
    all_glm_beta_A_ITD; all_glm_beta_A_noITD; ...
    all_glm_beta_F_ITD; all_glm_beta_F_noITD; ...
    all_glm_beta_B_ITD; all_glm_beta_B_noITD};



t = (0:5:25)';
plot_num = 1;
region = {'E right cIFS'; 'A left cIFS'; 'F right STG'; 'B left STG'};
for kk = 1:2:8
    subplot(2,2,plot_num);
    hold on;
    all_glm_beta_E_ITD = all_data{kk};
    all_glm_beta_E_noITD = all_data {kk+1};
    % ITD
    ITD_average = mean(all_glm_beta_E_ITD,1)';
    ITD_SE = std(all_glm_beta_E_ITD,0,1)';
    % noITD
    noITD_average = mean(all_glm_beta_E_noITD,1)';
    noITD_SE = std(all_glm_beta_E_noITD,0,1)';
    
    plot(t,noITD_average,'b','LineWidth',2);
    plot(t,ITD_average,'r','LineWidth',2);
    legend('Opposite ear','ITD','Location', 'best');
    
    %calcualte SE and plot error bar
    % ITD with err bar
    noITD_err = noITD_SE/sqrt(size(all_glm_beta_E_noITD,1));
    h = fill([t;flipud(t)],[noITD_average-noITD_err; flipud(noITD_average+noITD_err)],[0 0 1],'linestyle','none');
    set(h,'facealpha',.3);
    
    %calcualte SE and plot error bar
    % ITD with err bar
    ITD_err = ITD_SE/sqrt(size(all_glm_beta_E_noITD,1));
    h = fill([t;flipud(t)],[ITD_average-ITD_err; flipud(ITD_average+ITD_err)],[1 0 0],'linestyle','none');
    set(h,'facealpha',.3);
    title(['task at ', region{plot_num}]);
    plot_num = plot_num +1;
    
    
end




