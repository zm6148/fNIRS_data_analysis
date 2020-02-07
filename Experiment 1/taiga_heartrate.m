clear; close all;
currentFolder = pwd; HRFolder = fileparts(currentFolder);
Homer2Folder = [HRFolder '\FNIRS heartrate\Homer2\'];
TimWindowFolder = [HRFolder '\FNIRS heartrate\TimWindow\'];
addpath(genpath(Homer2Folder)); addpath(genpath(TimWindowFolder));

% heart rate extraction
% load one example data
%load('TJT_heart_heat.mat');
load('TJN_heart_heat.mat');
start_index = first_task_index-50*5;
end_index = last_task_index+50*45; 

%% define basic stuff
Fs = 50;
% Define window
window_b=[-5,40];
window_t=[-5,40];

% choose one channel and only task portion to study for now
% channels at B
RC = 10;
channels = [12, 9, 11, 7];%, 8];

% %channels at D
% RC = 19;
% channels = [20, 23, 21, 24];%, 22];

dod_ind = dod_raw(start_index:end_index,RC);
tx = t(start_index:end_index);

%% Find average heart rate
%Filter raw data in broad band centered at 1 Hz (wide filter just to remove sharp peaks and drift)

Hd = LPfilt_1_10;
%filtfilt(SOS,G,x) G is between 1 and L+1, L = 2
dod_filtWide = filtfilt(Hd.sosMatrix,1,dod_ind);

%Find peak of the autocorrelation function for each channel
%[acf,lags,bounds] = autocorr(dod_filtWide)
[r,lags] = xcorr(dod_filtWide);

[pks,locs] = findpeaks(r,'MinPeakDistance',2e3);

figure(1);
plot(r)
hold on
plot(locs,pks,'o')
yyaxis right
plot(lags)

%Average across channels
%Convert lag a peak to Hz, this is the average heart rate


%We are taking the peaks of the positive lag times
indexPks = ceil(length(pks)/2);
indexLocs = locs(indexPks+1:end);

peakLags = lags(indexLocs);
%avgLagPeaks = 1/(mean(peakLags/Fs));

%We are taking the difference between each lag
lag_diff = zeros(1,indexPks-2);
for i = 1:(indexPks)-2
    lag_diff(i) = peakLags(i+1) - peakLags(i);
end

avgLagDiff = 1/(mean(lag_diff/Fs));


%% Filter raw data in narrow band around average HR +/-20% of average
pass_minus = avgLagDiff * 0.8;
pass_plus = avgLagDiff * 1.2;

%Wn = [pass_minus/Fs, pass_plus/Fs];
Wn = [pass_minus, pass_plus];
[b,a] = butter(4,Wn);

% d = designfilt('bandpassiir','FilterOrder',4, ...
%     'HalfPowerFrequency1',pass_minus,'HalfPowerFrequency2',pass_plus, ...
%     'SampleRate',Fs);

dod_filtNarrow = filtfilt(b,a,dod_ind);

%Use findpeaks to determine time of each beat
[pks2,locs2] = findpeaks(dod_filtNarrow);%,'MinPeakDistance',2e3);

figure(2);
plot(dod_filtNarrow)
hold on
plot(locs2,pks2,'o')

%Convert beat times to inter-beat intervals (diff)
%Find Sample Interval between peaks
IBI_Sample = zeros(1,length(locs2)-1);
for i = 1:length(locs2)-1
    IBI_Sample(i) = locs2(i+1) - locs2(i);
end
%Find Time Interval using Sample Intervals
IBI_Time = IBI_Sample./Fs;

%Reject IBIs that are impossibly short or long
%Normal Resting Range: 60 - 100
%Set - Too slow: 30bpm, Too Fast: 200bpm
HR_maxbound = 200/60; %in beats per second
HR_minbound = 30/60;

%Reject by adding or subtracting a beat in between
for checkInt = 1:length(IBI_Time)
    if IBI_Time(checkInt) > HR_maxbound
        IBI_Time(checkInt) = IBI_Time(checkInt)/2;
    elseif IBI_Time(checkInt) < HR_minbound
        IBI_Time(checkInt) = IBI_Time(checkInt)*2;
    end
end

%Reject by Removing beat
% for checkInt = 1:length(IBI_Time)
%     if (IBI_Time(checkInt) > HR_maxbound || IBI_Time(checkInt) < HR_minbound)
%         IBI_Time(checkInt) = NaN;
%         %disp('found')
%     end
% end

figure(3);
plot(IBI_Time)

%% Find channel that has the best scalp-coupling index

%Stats on the IBIs from best channel