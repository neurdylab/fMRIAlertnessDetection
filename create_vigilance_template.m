%% Build Template
% 1.23.2021
% Sarah Goodale - author
% ELIFE paper: fMRI-based detection of alertness predicts electrophysiological and behavioral variability

% if using please cite: 

clear; clc; close all;
addpath('C:\Users\goodalse\Documents\MATLAB\NIfTI'); %add NIFTI toolbox path
addpath('C:\Users\goodalse\Documents\MATLAB\spm12'); %add spm toolbox path

% Load subjects - we have them as txt file
subs_to_do = importdata('E:\neurdy\Goodale\template\10282020\singlesub_ecrSubs_to_do.txt');

for n = 1:length(subs_to_do)
sub_id = subs_to_do{n}

%Load FMRI data
cd('E:\neurdy\PROC\Ymats_121419'); % our fMRI data loads as variable Y in timecourse x voxel matrices 
tmp_fmri = dir(sub_id);
load([tmp_fmri.folder,'\', tmp_fmri.name]);
nanY = find(isnan(Y)); 
Y(nanY) = 0;

%Load EEG data
cd('E:\neurdy\PROC\EEG_mats_122919');
tmp_eeg = dir(sub_id);
E = load([tmp_eeg.folder,'\', tmp_eeg.name]);

TR = 2.1;
eeg_reg0 = E.OUT.BL_RATS{1}; % calls our alpha/theta ratio EEG data
eeg_reg1 = eeg_reg0(8:end)'; % matches eeg with fMRI after dropping initial files in pre-processing
eeg_reg2 = eeg_reg1 - mean(eeg_reg1); % zero meaning to avoid weirdness at edge when convolution
eeg_conv = conv(eeg_reg2,spm_hrf(TR)); % convolve with spm HRF
eeg_reg = eeg_conv(1:length(eeg_reg1),:); % get rid of extra at the end
eeg_reg_bp = bandpass(eeg_reg,[0.01 0.2],1/TR); % bandpass EEG 0.01-0.2
  
% Create TEMPLATE
template_vec = corr(eeg_reg_bp,Y);
template_vec(find(isnan(template_vec)))=0;
template_vecz = atanh(template_vec'); % normalize data
all_temp_vecz(:,n) = template_vecz;

clearvars template_vec template_vecz templateVol
end

% Create mean template to use to make predictions in another dataset
mean_temp_vecz = double(mean(all_temp_vecz,2)); % average correlations across subjects to get mean template
