%% analysis on response

clear;

% load FAR lookup table
load('FAR_table.mat');

order_1= {'TNQ_active_exp1_repeat', ...
    'TNT_active_exp1_repeat', ...
    'TNU_active_exp1_repeat',...
    'TNV_active_exp1_repeat', ...
    'TNX_active_exp1_repeat', ...
    'TNY_active_exp1_repeat', ...
    'TNZ_active_exp1_repeat', ...
    'TOB_active_exp1_repeat', ...
    'TOC_active_exp1_repeat', ...
    'TOF_active_exp1_repeat', ...
    'alexi_active_exp1_repeat', ...
    'sravana_active_exp1_repeat', ...
    'steven_active_exp1_repeat'...
    'ww_active_exp1_repeat'};

order_2= {'TNQ_passive_exp1_repeat', ...
    'TNT_passive_exp1_repeat', ...
    'TNU_passive_exp1_repeat',...
    'TNV_passive_exp1_repeat', ...
    'TNX_passive_exp1_repeat', ...
    'TNY_passive_exp1_repeat', ...
    'TNZ_passive_exp1_repeat', ...
    'TOB_passive_exp1_repeat', ...
    'TOC_passive_exp1_repeat', ...
    'TOF_passive_exp1_repeat', ...
    'alexi_passive_exp1_repeat', ...
    'sravana_passive_exp1_repeat', ...
    'steven_passive_exp1_repeat'...
    'ww_passive_exp1_repeat'};

% 2 conditions
% 1: speech
% 2: noise

out_put = [];
for ii = 1:length(order_1)
    
    name = order_1{ii};
    
    disp(name);
    load(name);
    
    
    all_block_hitrate_SNR_cal = [];
    for jj = 1 : length(all_response_timestamp)
        % loop through each target play time and try to find a match within
        % n s (1.2 s for now) in the response timestamp
        correct_window = 1.2;
        
        target_timestamp = all_target_timestamp{jj};
        response_timestamp = all_response_timestamp{jj};
        block_total_color = length(target_timestamp);
        block_total_response = length(response_timestamp);
        
        block_correct_total = 0;
        
        for k = 1 : size(target_timestamp, 1)
            target_time = target_timestamp(k,2);
            % loop through each response time to find a match
            for l = 1 : size(response_timestamp,1)
                response_time = response_timestamp(l,3);
                % 			disp(target_time);
                % 			disp(response_time);
                % 			disp(response_timestamp);
                % 			disp(l);
                % if within the window
                if (target_time <= response_time) && (response_time <= (target_time + correct_window))
                    block_correct_total = block_correct_total + 1;
                    response_timestamp(l,:) = [];
                    break;
                end
            end
        end
        
        % calculate hitrate
        hitrate = block_correct_total/block_total_color;
        %calcualte d'
        % 	find hits and false alarm
        % 	simulated farate, depends on the number of keywords and button
        % 	pushes for that block and listener
        % 	(use Antje's simulation code to generate a look-up table)
        
        hits = block_correct_total;
        switch block_total_color
            case 3
                FAR_index= 1;
            case 4
                FAR_index = 2;
            case 5
                FAR_index = 3;
        end
        if block_total_response == 0
            farate = 0;% mean(FAR(FAR_index, 1,:));
        else
            if block_total_response > 9
                block_total_response = 9;
            end
            farate = FAR(FAR_index, block_total_response);
        end
        Ntrials = block_total_color;
        
        % d' instead of percentage correct
        % use the look up table to find corresponding farate
        %simulated farate, depends on the number of keywords and button pushes for that block and listener (use Antje's simulation code to generate a look-up table)
        hitrate = hits/Ntrials;% observed percent correct for that block and listener
        % clipping
        hitrate(hitrate > .999) = .999;
        hitrate(hitrate < .001) = .001;
        farate(farate > .999) = .999;
        farate(farate < .001) = .001;
        
        zh = erfinv(1-2.*(1-hitrate)).*sqrt(2);%z-score of the observed hit rate
        zf = erfinv(1-2.*(1-farate)).*sqrt(2);% z-score of the simulated false alarm rate
        dprime = zh - zf;%d prime score for that listener and block
        beta = -(zh+zf)/2;% bias for that listener and block
        
        
        % total color and total hits
        %      correct                                d'      hits                 total color        total response
        temp =[block_correct_total/block_total_color, dprime, block_correct_total, block_total_color, block_total_response];
        all_block_hitrate_SNR_cal = [all_block_hitrate_SNR_cal; temp];
        
    end
    
    % select from all_block_hitrate_SNR
    % correct%  d'  hits    total color   total response
    condition = mean(all_block_hitrate_SNR_cal,1);
    
    % final matrix
    out_put = [out_put; condition];
end
% index_1 = find(condition_sequence_rand == 1);
% index_2 = find(condition_sequence_rand == 2);
% index_3 = find(condition_sequence_rand == 3);
% index_4 = find(condition_sequence_rand == 4);
%
% % select from all_block_hitrate_SNR
% % correct%  d'  hits    total color   total response
% condition_original = mean(all_block_hitrate_SNR(index_1,:),1);
% condition_8_band = mean(all_block_hitrate_SNR(index_2,:),1);
% condition_16_band = mean(all_block_hitrate_SNR(index_3,:),1);
% condition_32_band = mean(all_block_hitrate_SNR(index_4,:),1);
%
% % final matrix
% final_matrix = cat(1, condition_original, condition_32_band, condition_16_band, condition_8_band);
%
