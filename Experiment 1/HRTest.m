clear; close all;
currentFolder = pwd; HRFolder = fileparts(currentFolder);
Homer2Folder = [HRFolder '\FNIRS heartrate\Homer2\'];
addpath(Homer2Folder);

% heart rate extraction
% load one example data
%load('TJT_heart_heat.mat');
load('dc_Hb.mat');
% start_index = first_task_index-50*5;
% end_index = last_task_index+50*45;

%% define basic stuff
% Define window
window_b=[-5,40];
window_t=[-5,40];

% chose one channel and only task portion to study for now
% channels at B
RC = 10;
channels = [12, 9, 11, 7];%, 8];
% channels at D
% RC = 19;
% channels = [20, 23, 21, 24];%, 22];

% dod_ind = dod_raw(start_index:end_index,RC);
% tx = t(start_index:end_index);

dod_ind = dc_Hb;

%Call the Butter IIR LPF fc=4Hz 
Hd = fnirLPFilt;
%filtfilt(SOS,G,x) G is between 1 and L+1, L = 2
dod_filt = filtfilt(Hd.sosMatrix,1,dod_ind);
%data_original2 = dod_filt;
%data_original = resample(dod_filt,tx,100);
data_original = dod_filt;

% Peaks in the HbO signal were found using the FINDPEAKS routine in the
% MATLAB signal processing toolbox, with a minimum spacing equivalent to
% 200 BPM
[pks,locs] = findpeaks(data_original,'MinPeakDistance',0.5*50);%,'MinPeakProminence', 1);

%% Heart Rate Calculation

%Find Sample Interval between peaks
RR_SampleDiff = zeros(1,size(locs,2)-1);
for i = 1:length(locs)-1
    RR_SampleDiff(i) = locs(i+1) - locs(i);
end
%Find Time Interval using Sample Intervals
RR_TimeDiff = RR_SampleDiff./100;

%Dropped Beats
%Intervals Longer than the mean + three standard deviations are assumed to
%be dropped beats and the intervals are divided in half
averageTimeDiff = mean(RR_TimeDiff); threeDev = std(RR_TimeDiff)*3;
maxInterval = averageTimeDiff+threeDev;
for checkInt = 1:size(RR_TimeDiff,2)
    if RR_TimeDiff(checkInt) > maxInterval
        RR_TimeDiff(checkInt) = RR_TimeDiff(checkInt)/2;
    end
end

%Divide Blocks into time blocks.
totalSamples = length(data_original); %totalSecs = totalSamples/100; %100Hz at resampled rate
%Divide into 10 second blocks: 1000 samples per block at 100Hz
%Find which peaks belong in which block. Check if sample is less than 1k, 2k, ... to 589k
samplesPerBlock = 50 * 4; 
numBlocks = ceil(totalSamples/samplesPerBlock);
sampleCheckPrev = 0;
index = cell([1,numBlocks]); Block = cell([1,numBlocks]);
instantHR = zeros(1,length(Block));

for i = 1:numBlocks
    sampleCheck = i*samplesPerBlock;
    index{1,i} = find(locs < sampleCheck & locs > sampleCheckPrev); %gives index of locs_corr

    if isempty(index{1,i}) %if there are no samples that match, replace with NaN
        index{1,i} = NaN;
    else
        if index{1,i}(end) >= size(RR_TimeDiff,2) %Avoid error when index is larger than size
            index{1,i}(end) = size(RR_TimeDiff,2);
        end
        Block{1,i} = RR_TimeDiff(index{1,i}(1):index{1,i}(end)); %Cell filled with Time Intervals
        instantHR(i) = (nanmean(1./Block{1,i})); %Take average of inverse of time intevals
    end
    sampleCheckPrev = sampleCheck;
end

% Interpolation of HRV using interp1(SamplePoints, SampleValues, QueryPoints)
% Interpolating to size of filtered data
x = 1:size(instantHR,2); 
xq = linspace(1,size(instantHR,2),size(instantHR,2)*samplesPerBlock);
HRV = interp1(x, instantHR, xq, 'spline');

%% Plots
%Filtered Data and Peaks
figure(); plot(data_original);
hold on;
y_limit=get(gca,'ylim');
x_limit=get(gca,'xlim');
plot(locs, pks, 'o');
title('Heart Rate Peaks on Filtered Data')

%Plot Calculated HRV
figure(); plot(instantHR);
title('Instantaneous HRV')

%Plot Interpolated HRV against filtered data
figure(); plot(data_original);
hold on;
plot(locs, pks, 'o');
title('Peaks and HRV')
yyaxis right
plot(HRV)