%% Vigilance Prediction
% 1.23.2021
% Sarah Goodale - author
% ELIFE paper: fMRI-based detection of alertness predicts electrophysiological and behavioral variability

% if using please cite: 

clear; clc; close all; 
addpath('C:\Users\goodalse\Documents\MATLAB\NIfTI'); %add NIFTI toolbox
addpath('C:\Users\goodalse\Documents\MATLAB\spm12'); %add spm toolbox

% Load Subjects and template
load('E:\neurdy\Goodale\template\752020\full_template.mat'); % template mat file (voxels x 1) created using create_template.m
subs_to_do = importdata('E:\neurdy\Goodale\template\10282020\singlesub_ectSubs_to_do.txt'); % text file with sub names
load('E:\neurdy\Goodale\template\752020\singlesub_exclude_info.mat'); % trials to exclude based on closeness criteria - see paper

% Initialize
TR = 2.1; % scanner TR
ROI_evts_all = []; hit_miss_all = []; vig_pred_subs = []; all_delay = [];

%% Calculate Vigilance Prediction
for n = 1:length(subs_to_do)
sub_id = subs_to_do{n}
%load in subject mat file

%EEG ----------------------------------------------------------------------
cd('E:\neurdy\PROC\EEG_mats_122919'); tmp_eeg = dir(sub_id); 
E = load([tmp_eeg.folder,'\', tmp_eeg.name]); nm = {tmp_eeg.name};
eeg_proc = E.OUT.eeg; % pull task eeg from OUT struct provided on github repo
%data provided is processed in same way as in create_template.m

all_eeg(:,n) = eeg_proc; % concatonate all eeg_proc

%TASK ---------------------------------------------------------------------
if ~isempty(tmp_task) 
delay_response = OUT.RT;
stim_time = OUT.stimtime; % in seconds
st_msec = stim_time.*1000; 
% add time_response

TR_dr = delay_response/TR; % the length of the delay in TRs
TR_st = stim_time/TR; % auditory tone time in TRs

%Get TR's lined up with eeg
if TR_st(1) > 8
    % all stays the same
    st_msec = st_msec(1:end-1);
    time_response = time_response(1:end-1);
    delay_response = delay_response(1:end-1);
    stim_time = stim_time(1:end-1);
    TR_dr = TR_dr(1:end-1); TR_st = TR_st(1:end-1);
    disp('same');
elseif TR_st(1) < 8
    st_msec = st_msec(2:end);
    time_response = time_response(2:end);
    delay_response = delay_response(2:end);
    stim_time = stim_time(2:end);
    TR_dr = TR_dr(2:end); TR_st = TR_st(2:end);
    
    st_msec = st_msec(1:end-1);
    time_response = time_response(1:end-1);
    delay_response = delay_response(1:end-1);
    stim_time = stim_time(1:end-1);
    TR_dr = TR_dr(1:end-1); TR_st = TR_st(1:end-1);
end
end

all_RT{n} = delay_response; % collect all reaction times for figs
all_stimtime{n} = stim_time; % collect all stim times

%FMRI ---------------------------------------------------------------------
cd('E:\neurdy\PROC\Ymats_121419'); % load fMRI
tmp_fmri = dir(sub_id);
load([tmp_fmri.folder,'\', tmp_fmri.name]);
nanY = find(isnan(Y)); 
Y(nanY) = 0; 
zY = zscore(Y)'; % z-score fmri data
all_fmri(:,:,n) = Y;

%TEMPLATE -----------------------------------------------------------------
vig_pred = corr(mean_temp_vecz,zY); % correlate template with fmri data
all_vigpred(n,:) = vig_pred; % concat all vigilance predictions

% FMRI FINE SAMPLING ------------------------------------------------------
fmri_WIN = [-5,0]; % window around event in seconds (0 being the stimulus onset)
dt_fmri = 0.840; % new finer sampling rate

% interpolate timepoints to create finer time axis
 [ROI_evts, ROI_evts_fine] = fmri_finesampling_SG(vig_pred',fmri_WIN,dt_fmri,st_msec);
 ROI_evts_ts = [fmri_WIN(1):dt_fmri:fmri_WIN(2)];

% DISCERN HITS AND MISSES -------------------------------------------------
remove = exclude_info.excludeTrials{1, n}; % remove trials too close in time - see paper
hit_miss_idx = ~isnan(time_response); % did they respond 1:yes, 0:no  where time response is in seconds

% Pull out time series for hits and misses separate based on response and
% threshold
for i = 1:length(hit_miss_idx)
    if hit_miss_idx(i) == 1
        sub_hits(i,:) = ROI_evts(i,:);
        if delay_response(i) < .565 % all trial median threshold
            fast_hits(i,:) = ROI_evts(i,:);
            fast_dr(i) = delay_response(i);
        elseif delay_response(i) > .565 % all trial median threshold
            slow_hits(i,:) = ROI_evts(i,:);
            slow_dr(i) = delay_response(i);
        end
    elseif hit_miss_idx(i) == 0
        sub_miss(i,:) = ROI_evts(i,:);
    end
end

% Find only where they hit
r_idx = find(remove<=size(sub_hits,1));
sub_hits(remove(r_idx),:) = [];
nz_subhits = unique(sub_hits,'rows','stable');
nz_subhits(nz_subhits==0) = NaN;
sub_hit_avg = nanmean(nz_subhits);
allsub_hitavg(n,:) = sub_hit_avg;
total_hit{n} = nz_subhits; % hit time courses across all subjects

% Breakdown slow vs fast
if exist('fast_hits') == 1
    allfastdr{n} = fast_dr; 
r_idx = find(remove<=size(fast_hits,1));
fast_hits(remove(r_idx),:) = [];    
nz_fasthits = unique(fast_hits,'rows','stable');
nz_fasthits(nz_fasthits==0) = NaN;
sub_fasthit_avg = nanmean(nz_fasthits);
allsub_fasthitavg(n,:) = sub_fasthit_avg;
total_fasthit{n} = nz_fasthits;
end

if exist('slow_hits') == 1
    allslowdr{n} = slow_dr;
r_idx = find(remove<=size(slow_hits,1));
slow_hits(remove(r_idx),:) = [];
nz_slowhits = unique(slow_hits,'rows','stable');
nz_slowhits(nz_slowhits==0) = NaN;
sub_slowhit_avg = nanmean(nz_slowhits);
allsub_slowhitavg(n,:) = sub_slowhit_avg;
total_slowhit{n} = nz_slowhits;
end

% Find only where they miss, if they missed
if exist('sub_miss') == 1
r_idx = find(remove<=size(sub_miss,1));
sub_miss(remove(r_idx),:) = [];
nz_submiss = unique(sub_miss,'rows','stable');
nz_submiss(nz_submiss==0) = NaN;
sub_miss_avg = nanmean(nz_submiss);
allsub_missavg(n,:) = sub_miss_avg;
total_miss{n} = nz_submiss; % miss time courses across all subjects
end

% SUBJECT AVGS AND PLOTS --------------------------------------------------
% Subject Overall Averages
sub_avg = mean(ROI_evts);
all_sub_avg(n,:) = sub_avg;

% Cross correlation to see how vigilance prediction matches subject EEG
[cross_corr,lags] = xcov(eeg_proc,vig_pred',20,'coeff');
max_corr = max(cross_corr);
all_maxcorr(n) = double(max_corr);
all_xcorr(:,n) = cross_corr;

% plot vigilance prediction timecourse on eeg
figure; plot(eeg_proc,'b','LineWidth',1); hold on; plot(vig_pred,'r','LineWidth',1);

% cleanup for each iteration
clearvars hit_miss_idx sub_hits sub_miss nz_subhits nz_submiss fast_hits slow_hits high_miss
clearvars nz_fasthits nz_slowhits fast_dr slow_dr hit_xcov miss_xcov high_hit hit_miss_dr

end
