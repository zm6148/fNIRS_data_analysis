%% load data folder
% data analysis using GLM model
clear;
%close all;

% speech_path = uigetdir(pwd,'choose speech reuslts folder');
% speech_filepattern=fullfile(speech_path,'/*.mat');
% speech_files=dir(speech_filepattern);
% 
% SSN_path = uigetdir(pwd,'choose SSN reuslts folder');
% SSN_filepattern=fullfile(SSN_path,'/*.mat');
% SSN_files=dir(SSN_filepattern);

speech_path = [pwd,'\Data\Speech\Processed'];
speech_filepattern=fullfile(speech_path,'/*.mat');
speech_files=dir(speech_filepattern);

SSN_path = [pwd,'\Data\SSN\Processed'];
SSN_filepattern=fullfile(SSN_path,'/*.mat');
SSN_files=dir(SSN_filepattern);



%% define basic stuff
% Difine window
window_b=[-5,40];
window_t=[-5,40];

% build Refrence channel regressor
% for source C RC: 18; channels: 13, 16, 15, 17
% for source D RC: 19; channels: 20, 23, 21, 24
% for source A RC: 4; channels: 2, 6, 3, 1
% for source B RC: 10; channels: 12, 9, 11, 7
% use regressor closest to the source

% RC_C = 18;
% channels_C = [13, 16, 15, 17];%, 14];
% 
% RC_A = 4;
% channels_A = [2, 6, 3, 1];%,, 5];
% 
% RC_D = 19;
% channels_D = [20, 23, 21, 24];%,, 22];
% 
% RC_B = 10;
% channels_B = [12, 9, 7, 11];%,, 8];

RC_C = 18;
RC_A = 4;
RC_D = 19;
RC_B = 10;

channels_C = 16; %[13, 15, 17, 16];
channels_A = 1;  %[3,  2,  6,  1];
channels_D = 23; %[20, 21, 24, 23];
channels_B = 11; %[7,  12, 9,  11]];

% combine CA and DB
RC_CA = [18,4];
channels_CA = [13, 16, 15, 17, 2, 6, 3, 1];

RC_DB = [19, 10];
channels_DB = [12, 9, 11, 7, 20, 23, 21, 24];

% combine CD (left) and AB (right)
RC_CD = [18, 19];
channels_CD = [13, 16, 15, 17, 20, 23, 21, 24];

RC_AB = [4, 10];
channels_AB = [2, 6, 3, 1, 12, 9, 11, 7];
%% load individual data

%(C, A, D, B, CA, DB, CD, AB) * (17 subjects)

%%
subject_beta_speech = [];
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

speech_blavg_max = [];
SSN_blavg_max = [];

Hb_type = 'HbO'; %'HbR'
sys_or_fun = 'fun'; %'sys'

for index=1:length(speech_files)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cat breath holding part from speech and SSN recordings
    
    % first load the speech part of breathholding
    base_name = speech_files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [speech_path,'\',name,'.mat'];
    load(new_name);
    disp(base_name);
    dc_speech_b = dc;
    % implement HMS before doing any block average
    new_dc = [];
    for ii = 1:size(dc,3)
        [dc_HbO_f, dc_HbR_f] = HMS(dc(:,1, ii), dc(:,2, ii), sys_or_fun);
        temp = cat(2, dc_HbO_f', dc_HbR_f', dc_HbO_f'+ dc_HbR_f');
        new_dc = cat(3, new_dc, temp);
    end
    % assign it back to old dc
    % dc = new_dc;
    
    aux_speech = aux;
    dc_speech = dc;
    dc_task_speech = dc_task;
    first_task_index_speech = first_task_index;
    last_task_index_speech = last_task_index;
    s_b_speech = s_b;
    s_ITD_task_speech = s_ITD_task;
    s_noITD_task_speech = s_noITD_task;
    s_task_speech = s_ITD_task_speech + s_noITD_task_speech;
    t_speech = t;
    
    % then load the SSN part of the breathholding
    base_name = SSN_files(index).name;
    [folder, name, extension] = fileparts(base_name);
    new_name = [SSN_path,'\',name,'.mat'];
    load(new_name);
    disp(base_name); 
    dc_SSN_b = dc;
    % implement HMS before doing any block average
    new_dc = [];
    for ii = 1:size(dc,3)
        [dc_HbO_f, dc_HbR_f] = HMS(dc(:,1, ii), dc(:,2, ii),sys_or_fun);
        temp = cat(2, dc_HbO_f', dc_HbR_f', dc_HbO_f'+ dc_HbR_f');
        new_dc = cat(3, new_dc, temp);
    end
    % assign it back to old dc
    % dc = new_dc;
    
    aux_SSN = aux;
    dc_SSN = dc;
    dc_task_SSN = dc_task;
    first_task_index_SSN = first_task_index;
    last_task_index_SSN = last_task_index;
    s_b_SSN = s_b;
    s_ITD_task_SSN = s_ITD_task;
    s_noITD_task_SSN = s_noITD_task;
    s_task_SSN = s_ITD_task_SSN + s_noITD_task_SSN;
    t_SSN = t;
    
    % block average at breahhodling part
    breath_blockAvg_speech = hmrBlockAvg(dc_speech,s_b_speech,t_speech,window_b);
    breath_blockAvg_SSN = hmrBlockAvg(dc_SSN,s_b_SSN,t_SSN,window_b);
        
    % concatanate breath
    first_index_speech = min(find(s_b_speech==1));
    last_index_speech = max(find(s_b_speech==1));
    
    first_index_SSN = min(find(s_b_SSN==1));
    last_index_SSN = max(find(s_b_SSN==1));
    
    dc_b_speech = dc_speech_b(first_index_speech : last_index_speech+45*50,:,:);
    s_b_speech = s_b_speech(first_index_speech : last_index_speech+45*50);
    
    dc_b_SSN = dc_SSN_b(first_index_SSN : last_index_SSN+45*50,:,:);
    s_b_SSN = s_b_SSN(first_index_SSN : last_index_SSN+45*50);
    
    dc_b_both = cat(1, dc_b_speech, dc_b_SSN);
    s_b_both = cat(1, s_b_speech, s_b_SSN);
    
%     dc_b_both = dc_b_speech;
%     s_b_both =  s_b_speech;
    
    % concatanate task part
    % speech first then SSN
    % speech 0 before, 45 after
    % SSN 0 before, 45 after
    
    first_index_speech_task = min(find(s_task_speech==1));
    last_index_speech_task = max(find(s_task_speech==1));
    
    first_index_SSN_task = min(find(s_task_SSN==1));
    last_index_SSN_task = max(find(s_task_SSN==1));
    
    dc_task_speech = dc_speech(first_task_index_speech : last_task_index_speech+45*50, :, :);
    dc_task_SSN = dc_SSN(first_task_index_SSN : last_task_index_SSN+45*50, :, :);
    
    s_task_speech = s_task_speech(first_index_speech_task : last_index_speech_task+45*50);
    s_task_SSN = s_task_SSN(first_index_SSN_task : last_index_SSN_task+45*50);
    
    dc_task_both = cat(1, dc_task_speech, dc_task_SSN);
    s_task_speech_both = cat(1, s_task_speech, 0.*s_task_SSN);
    s_task_SSN_both = cat(1, 0.*s_task_speech, s_task_SSN);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % GLM analysis on HbO of the breath holding data
    [ b_C_breath, r_C_breath, p_C_breath] = RCS_breath_catselection( dc_b_both, s_b_both, RC_C, channels_C, Hb_type);
    [ b_A_breath, r_A_breath, p_A_breath] = RCS_breath_catselection( dc_b_both, s_b_both, RC_A, channels_A, Hb_type);
    [ b_D_breath, r_D_breath, p_D_breath] = RCS_breath_catselection( dc_b_both, s_b_both, RC_D, channels_D, Hb_type);
    [ b_B_breath, r_B_breath, p_B_breath] = RCS_breath_catselection( dc_b_both, s_b_both, RC_B, channels_B, Hb_type);
    
    
    % GLM analysis on HbO of the ITD and no_ITD
    [ b_C, r_C, p_C] = RCS_all_stim_cat(  dc_task_both, s_task_speech_both, s_task_SSN_both, RC_C, channels_C, breath_blockAvg_speech, breath_blockAvg_SSN, Hb_type);
    [ b_A, r_A, p_A] = RCS_all_stim_cat(  dc_task_both, s_task_speech_both, s_task_SSN_both, RC_A, channels_A, breath_blockAvg_speech, breath_blockAvg_SSN, Hb_type);
    [ b_D, r_D, p_D] = RCS_all_stim_cat(  dc_task_both, s_task_speech_both, s_task_SSN_both, RC_D, channels_D, breath_blockAvg_speech, breath_blockAvg_SSN, Hb_type);
    [ b_B, r_B, p_B] = RCS_all_stim_cat(  dc_task_both, s_task_speech_both, s_task_SSN_both, RC_B, channels_B, breath_blockAvg_speech, breath_blockAvg_SSN, Hb_type);
    
    
    % save data
    % build a matrix that has all data for all subjects
    % speech b(2); SSN b(3); RC b(4)
    subject_beta_speech = cat(2, b_C, b_A, b_D, b_B);%, b_CA, b_DB, b_CD, b_AB);
    all_glm_beta = cat(3, all_glm_beta, subject_beta_speech);
    subject_r2 = cat(1, r_C, r_A, r_D, r_B); %, r_CA, r_DB, r_CD, r_AB);
    all_glm_r2 = cat(2, all_glm_r2, subject_r2);
    subject_p = cat(2, p_C, p_A, p_D, p_B); %, p_CA, p_DB, p_CD, p_AB);
    all_glm_p = cat(3, all_glm_p, subject_p);
    
    % save data breath holding datas
    % build a matrix that has all data for all subjects
    % speech b(2); SSN b(3); RC b(4)
    subject_beta_breath = cat(2, b_C_breath, b_A_breath, b_D_breath, b_B_breath); %, b_CA_breath, b_DB_breath, b_CD_breath, b_AB_breath);
    all_glm_beta_breath = cat(3, all_glm_beta_breath, subject_beta_breath);
    subject_r2_breath = cat(1, r_C_breath, r_A_breath, r_D_breath, r_B_breath); %, r_CA_breath, r_DB_breath, r_CD_breath, r_AB_breath);
    all_glm_r2_breath = cat(2, all_glm_r2_breath, subject_r2_breath);
    subject_p_breath = cat(2, p_C_breath, p_A_breath, p_D_breath, p_B_breath);%, p_CA_breath, p_DB_breath, p_CD_breath, p_AB_breath);
    all_glm_p_breath = cat(3, all_glm_p_breath, subject_p_breath);
end

% save all subject GLM data
% subject_data_name = [path, '/GLM/', 'GLM_data.mat'];
% save(subject_data_name,'all_glm_beta','all_glm_p','all_glm_r2', ...
%     'all_glm_beta_breath','all_glm_p_breath','all_glm_r2_breath');

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
    good_r2_index = find(all_glm_r2_breath(ii,:) > 0); %0.5
    % p vlaue for breath
    good_pvalue_index = find(all_glm_p_breath(2,ii,:) < 50);%0.05
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

% at each region select only good subjects based on breath holding
% table row : C A D B CA DB CD AB
% table column : subject_index, subject_beta_ITD', subject_beta_noITD', subject_pvalue_ITD', subject_pvalue_noITD', subject_r2'
for ii = 1:length(all_good_fit_subject_index)
    
    subject_index = all_good_fit_subject_index{ii};
    
    subject_r2 = all_glm_r2(ii,subject_index);
    
    n = length(subject_index);
    
    %speech
    subject_pvalue_speech = all_glm_p(2,ii,subject_index);
    subject_pvalue_speech = reshape(subject_pvalue_speech,[1,n]);
    
    subject_beta_speech = all_glm_beta(2,ii,subject_index);
    subject_beta_speech = reshape(subject_beta_speech,[1,n]);
    
    %SSN
    subject_pvalue_SSN = all_glm_p(3,ii,subject_index);
    subject_pvalue_SSN = reshape(subject_pvalue_SSN,[1,n]);
    
    subject_beta_SSN = all_glm_beta(3,ii,subject_index);
    subject_beta_SSN = reshape(subject_beta_SSN,[1,n]);
    
    subject_selection_by_breath_table{ii} = cat(2,subject_index, subject_beta_speech', subject_beta_SSN', subject_pvalue_speech', subject_pvalue_SSN', subject_r2');
end

%% further the selection based on task portion model

for ii = 1:length(all_good_fit_subject_index)
    
    results_region = subject_selection_by_breath_table{ii};
    all_index = 1:(size(results_region,1));
    
    % exclude p value not smaller than 0.05
    % for ITD and no_ITD both
    % r2_exclude = find(results_region(:,6) < 0); % exclude r2
    p_ITD_exclude = find(results_region(:,4) > 0.05); % exclude p ITD
    p_noITD_exclude = find(results_region(:,5) > 0.05); % exclude p noITD
    
    exclude_index = union(p_ITD_exclude, p_noITD_exclude);
    % exclude_index = union(exclude_index, r2_exclude);
    all_index(exclude_index)=[];
       
 
   subject_selection_by_breath_table_further{ii} = results_region(all_index,:);
end

%% save the selected subject code as in a file

% region C
channels_at_C = subject_selection_by_breath_table{1}(:,1);
% region A
channels_at_A = subject_selection_by_breath_table{2}(:,1);
% region D
channels_at_D = subject_selection_by_breath_table{3}(:,1);
% region B
channels_at_B = subject_selection_by_breath_table{4}(:,1);

selected_channels = {channels_at_C; channels_at_A; channels_at_D; channels_at_B};

% save selected channels
% subject_data_name = ['C:\Users\mz86\Desktop\fNIRS analysis\experiment 1\new Control data\', 'selected_channels_cat.mat'];
% save(subject_data_name,'selected_channels');
%% plot and see breath portion
% different region may have different subjects
% 2 cases:1, only select by breath. 2, further selection by ITD noITD
figure(1);
title('beta values at recording site');
names = {'C (left tgPCS) ';...
    'A (right tgPCS)';...
    'D (left cIFS)';...
    'B (right cIFS)';...
    'CA (tgPCS)';...
    'DB (cIFS)';...
    'CD (left)';...
    'AB (right)'
    };

breath_plot_data = [];
for ii = 1:length(all_good_fit_subject_index)
    
    subplot(2,2,ii);
    
    results_region = subject_breath_table{ii};
    
    
    b_beta = results_region(:,2);
    
    % t-test
    [h,p,ci,stats] = ttest(b_beta);
    
    % signed rank test
    %[p,h,stats_b] = signrank(b_beta);
    
    
    beta_mean = [mean(b_beta)];
    beta_SEM = [std(b_beta)/sqrt(length(b_beta))];
    
    breath_plot_data = [breath_plot_data; [beta_mean, beta_SEM, p]];
    
    hold on
    bar(beta_mean, 'r');
    errorbar(beta_mean,beta_SEM,'.');
    set(gca,'XTick', [0 1 2]);
    set(gca,'XTickLabel',{'', 'Breath Holding', ''})
    title(['breath ', names{ii}, ';p= ', num2str(p), ';tstats=', num2str(stats.tstat),';df=',num2str(stats.df)]);
    ylim([-6*10^(-7) 6*10^(-7)]);
end
% set(gcf, 'PaperUnits', 'inches');
% x_width=8 ;y_width=8;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); 

%% plot in graph all breath beta value plot
% figure;
% hold on;
% x = [1:1:4]';
% y = breath_plot_data(:,1);
% bar(x, y,'g');
% for i1=1:numel(y)
%     text(x(i1), y(i1), num2str(breath_plot_data(i1,3),'%10.3e'),...
%         'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom');
% end
% errorbar(breath_plot_data(:,1),breath_plot_data(:,2),'.');
% set(gca,'XTick', [0 1 2 3 4]);
% set(gca,'XTickLabel',{'', 'C', 'A', 'D', 'B'});
% set(gca,'TickDir','out');
% ylim([-6*10^(-7) 1*10^(-7)]);
% title('breath holding task beta value');



%% plot and see task portion
% different region may have different subjects
% 2 cases:1, only select by breath. 2, further selection by ITD noITD
which = 1;
%title('beta values at recording site');
names = {'C (left tgPCS) ';...
    'A (right tgPCS)';...
    'D (left cIFS)';...
    'B (right cIFS)';...
    'CA (tgPCS)';...
    'DB (cIFS)';...
    'CD (left)';...
    'AB (right)'
    };
task_speech_plot_data = [];
task_SSN_plot_data = [];
for ii = 1:length(all_good_fit_subject_index)
    
    %subplot(2,2,ii);
    switch which
        case 1
            results_region = subject_selection_by_breath_table{ii};
        case 2
            results_region = subject_selection_by_breath_table_further{ii};
    end
    
    speech_beta = results_region(:,2);
    SSN_beta = results_region(:,3);
    
    % signed rank test agansit meidan 0
%     [p_speech_beta,h_speech_beta,stats_speech] = signrank(speech_beta);%,0,'method','approximate');
%     [p_SSN_beta,h_SSN_beta,stats_SSN] = signrank(SSN_beta);%,0,'method','approximate');
%     [p_c_s,h_c_s, stats_2] = signrank(speech_beta,SSN_beta);%,'method','approximate');
    
    % t-test agansit meidan 0
    [h_speech, p_speech, ci_speech, stats_speech ] = ttest(speech_beta);
    [h_SSN ,p_SSN, ci_SSN, stats_SSN] = ttest(SSN_beta);
    [h_c_s, p_c_s, ci_c_s, stats_c_s] = ttest(speech_beta,SSN_beta);
    
    % t-test
    %[h,p,ci,stats] = ttest2(ITD_beta,noITD_beta);
    
    beta_mean = [mean(speech_beta), mean(SSN_beta)];
    beta_SEM = [std(speech_beta)/sqrt(length(speech_beta)), std(SSN_beta)/sqrt(length(SSN_beta))];

    % add adjustement for rANOVA
    diff_from_mean_speech = kron(ones(size(speech_beta,1),1),mean(speech_beta,1)) - kron(ones(size(speech_beta)),mean(mean(speech_beta,1),2));
    adjusted_data_speech = speech_beta - diff_from_mean_speech;
    adjusted_std_error_speech = std(adjusted_data_speech)./sqrt(length(speech_beta));
    
    diff_from_mean_SSN = kron(ones(size(SSN_beta,1),1),mean(SSN_beta,1)) - kron(ones(size(SSN_beta)),mean(mean(SSN_beta,1),2));
    adjusted_data_SSN = SSN_beta - diff_from_mean_SSN;
    adjusted_std_error_SSN = std(adjusted_data_SSN)./sqrt(length(SSN_beta));
    % end of adjustment

    % save for plot later
%     speech_data = [mean(speech_beta), std(speech_beta)/sqrt(length(speech_beta)), p_speech ];
%     SSN_data = [mean(SSN_beta), std(SSN_beta)/sqrt(length(SSN_beta)), p_SSN];

    speech_data = [mean(adjusted_data_speech),  adjusted_std_error_speech, p_speech ];
    SSN_data = [mean(adjusted_data_SSN), adjusted_std_error_SSN, p_SSN];
    
    task_speech_plot_data = [task_speech_plot_data; speech_data];
    task_SSN_plot_data = [task_SSN_plot_data; SSN_data];
    
%     hold on
%     bar(1,beta_mean(1),'r');
%     bar(2,beta_mean(2),'y');
%     errorbar(1:2,beta_mean,beta_SEM,'.');
%     set(gca,'XTick', [0 1 2 3]);
%     set(gca,'XTickLabel',{'', 'speech', 'SSN', ''})
%     title(['Task ', names{ii}, ' p = ', num2str(p_c_s), 'p speech = ', num2str(p_speech), 'p SSN = ',num2str(p_SSN)]);
%     ylim([-0.05 0.45]);
end

%% plot all 4 region in one bar plot
figure(2);
hold on;
x = [1:0.5:2.5]';
y = cat(2, task_speech_plot_data(:,1), task_SSN_plot_data(:,1));
bar(x, y, 'BarWidth',0.9);
%ylim([-0.05 0.45]);

% for i1=1:4
%     x_location = [1-0.075, 1.5-0.075, 2-0.075, 2.5-0.075];
%     text(x_location(i1), y(i1), num2str(task_speech_plot_data(i1,3),'%10.3e'),...
%         'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom');
% end
% 
% for i1=1:4
%     x_location = [1+0.075, 1.5+0.075, 2+0.075, 2.5+0.075];
%     text(x_location(i1), y(i1), num2str(task_SSN_plot_data(i1,3),'%10.3e'),...
%         'HorizontalAlignment','center',...
%         'VerticalAlignment','bottom');
% end

% adjusted error bar
% diff_from_mean_speech = kron(ones(size(task_speech_plot_data(:,1),1),1),mean(task_speech_plot_data(:,1),1)) - kron(ones(size(task_speech_plot_data(:,1))),mean(mean(task_speech_plot_data(:,1),1),2));
% adjusted_data_speech = task_speech_plot_data(:,1) - diff_from_mean_speech;
% adjusted_std_error_speech = std(adjusted_data_speech,0,2)./sqrt(length(speech_beta));
% 
% diff_from_mean_SSN = kron(ones(size(task_SSN_plot_data(:,1),1),1),mean(task_SSN_plot_data(:,1),1)) - kron(ones(size(task_SSN_plot_data(:,1))),mean(mean(task_SSN_plot_data(:,1),1),2));
% adjusted_data_SSN = task_SSN_plot_data(:,1) - diff_from_mean_SSN;
% adjusted_std_error_SSN = std(adjusted_data_SSN,0,2)./sqrt(length(SSN_beta));
% end of adjusted error bar

errorbar([1-0.075, 1.5-0.075, 2-0.075, 2.5-0.075], task_speech_plot_data(:,1),task_speech_plot_data(:,2),'.');
errorbar([1+0.075, 1.5+0.075, 2+0.075, 2.5+0.075],task_SSN_plot_data(:,1),task_SSN_plot_data(:,2),'.');

set(gca,'XTick', [0.5 1 1.5 2 2.5 3]);
set(gca,'XTickLabel',{'', 'C', 'A', 'D', 'B',''});
set(gca,'TickDir','out');
title('perceptual task beta value');
legend('Speech', 'SSN');


%% new plot with beta value for each area
% C_b_beta = [-3.67285622292623e-07;-2.64825439692080e-07;-1.57645458533783e-07;-3.89084146976586e-07;-8.47443878232698e-08;-1.01240101932271e-07;-7.17024287946947e-08;-2.00629289280951e-07;-1.69572354683513e-07;-1.70142580523898e-07;-1.67460324921895e-07;-3.29071437156442e-07;-2.46165725321329e-07];
% A_b_beta = [-6.09965770146579e-07;-5.10776918073187e-07;-3.77678682148899e-07;-4.53236587563945e-07;-1.30314727080577e-07;-1.41556047549662e-07;-6.22056242740678e-08;-2.01245855285846e-07;-2.75050826254769e-07;-3.60401894206441e-07;-1.93216847855377e-07;-4.80556297287418e-07;-3.32101138311541e-07];
% D_b_beta = [-4.91789423556818e-07;-3.36196167612743e-07;-2.85931137226238e-07;-7.81962899815182e-08;4.73441081576513e-08;-1.78952036347796e-07;-1.34982388550137e-07;-4.55943262997277e-07;-7.91838562375741e-08;-3.27612860990100e-07;-1.47919022472556e-07;-3.05064570346488e-07;-2.87824827665189e-07];
% B_b_beta = [-5.88455900804928e-07;-5.28203497270147e-07;-2.88615923796249e-07;5.65583717790143e-08;-6.81666593301597e-09;-1.36538764496955e-07;1.75134897407519e-08;-8.67924360246317e-07;-1.00639683935361e-07;-3.38904073000026e-07;-2.49803995250908e-07;-5.70810030250352e-07;-3.90354461079148e-07];
% b_all_beta_speech = cat(2, C_b_beta, A_b_beta, D_b_beta, B_b_beta);
% 
% C_b_beta = [-1.42078701375940e-07;-1.19231904307760e-07;-2.35238867924441e-07;-4.55334291315551e-07;-3.20561573992192e-08;-2.99765420172657e-07;-1.44468445010949e-07;1.56563032949974e-07;-2.13737826728944e-07;-1.10797695002658e-07;-2.01307379657301e-07;4.05187330138266e-09;-1.61373525665176e-07];
% A_b_beta = [-4.48989001632162e-08;-8.65925464914778e-09;-4.63135675010965e-07;-3.27115122205474e-07;-3.80010295273651e-10;-3.34064471123838e-08;-8.52671910143658e-08;6.34740999173338e-08;-2.12756945301526e-07;-2.94049438875944e-07;-1.69587611737083e-07;1.09755426538892e-08;-2.29333009820592e-07];
% D_b_beta = [-7.97525293068148e-08;-8.65304007441954e-08;-2.75417420629025e-08;-3.03948780842158e-07;-1.79419853195988e-08;-5.77045324294982e-07;-6.03939210323800e-07;2.34135777039904e-07;-2.90000991929262e-08;-4.66210840441886e-07;-4.03806449829259e-07;1.60376328367337e-08;-2.92288848473424e-07];
% B_b_beta = [8.65817377138668e-08;7.30278761592565e-08;-1.47932302245623e-07;-3.67126947874465e-07;-2.97752743023618e-08;-1.95155035464108e-07;-5.68623219159586e-07;2.24347592803028e-07;-5.15426087018543e-09;-1.94076618510319e-07;-6.29394525191639e-07;-6.92138361230267e-08;-4.83799760059985e-08];
% b_all_beta_SSN = cat(2, C_b_beta, A_b_beta, D_b_beta, B_b_beta);
% 
% figure;
% for ii = 1:length(all_good_fit_subject_index)
%     
%     subplot(2,2,ii);
%     results_region = subject_selection_by_breath_table{ii};
%     speech_beta = results_region(:,2);
%     SSN_beta = results_region(:,3);
%     
%     b_results_region = subject_breath_table{ii};
%     
%     
%     b_beta = b_results_region(:,2);
%     
%     
%     axis square;
%     hold on;
%     % 1 d
% %     scatter(abs(b_all_beta_speech(:,ii)), abs(speech_beta), 'r', 'o');
% %     scatter(abs(b_all_beta_SSN(:,ii)), abs(SSN_beta), 'g', '^');
%     plot(speech_beta, SSN_beta, '.', 'MarkerSize',20);
%     
%     xlabel('speech beta');
%     ylabel('SSN beta');
%     xlim([-2E-7, 6E-7]);
%     ylim([-2E-7, 6E-7]);
%     
%     m = 1; % slop
%     b = 0; % intersept
%     h = refline(m,b);
%     h.Color = 'k';
% end


%% save all GLM beta data
GLM_data_name = [pwd, '\GLM_data_exp1.mat'];
save(GLM_data_name,'subject_selection_by_breath_table');




