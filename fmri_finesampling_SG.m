function [ROI_evts, fax_TR_evt_s] = fmri_finesampling_SG(timecourse,fmri_WIN,dt_fmri, stimTime)
% 1.23.2021
% Catie Chang author with slight edits Sarah Goodale

% ELIFE paper: fMRI-based detection of alertness predicts electrophysiological and behavioral variability

% citation: 

% Inputs
% -- timecourse : timecourse to interpolate ie. vigilance prediction or eeg
% -- fmri_WIN : window to sample from timecourse (default -10:15)
% -- dt_fmri : new finer sampling rate (default 0.840)
% -- stimTime : OUT.stimTime_msec from task data

TR = 2.1;

    nframes = size(timecourse,1)+7; % 'original' nframes  %700
    fax0 = [0:TR:(nframes-1)*TR];
    fax_TR = fax0(8:end); % fmri time axis in TRs
    tmax = fax_TR(end); % end of fmri scan, in sec
    % also starts at 14.7sec
    fax_fine = [fax_TR(1):dt_fmri:fax_TR(end)];
    
for e=1:length(stimTime)
        t_evt =  stimTime(e)/1000; % sec
                    
        % skip first
        if t_evt==0
            continue;
        end
        
        % absolute time of stimulus event
        t_pre = max(t_evt+fmri_WIN(1),0);
        t_post = min(t_evt+fmri_WIN(2),tmax);
        
        % find the inds in the fine fMRI axis that correspond to this stim event
        fax_fine_evt_inds = intersect(find(fax_fine>t_pre),find(fax_fine<t_post));
        fax_fine_evt_s = fax_fine(fax_fine_evt_inds); % corresponding times in sec
                                                  
        % get original (TR-sampled) fMRI indices within the same event intervals
        fax_TR_evt_inds = intersect(find(fax_TR>t_pre),find(fax_TR<t_post));
        fax_TR_evt_s = fax_TR(fax_TR_evt_inds); % times in sec
      
        % get this section of the ROI data
        RZ_evt = timecourse(fax_TR_evt_inds,:);
              
        % resample onto finer grid
        RZ_evt_i = interp1(fax_TR_evt_s,RZ_evt,fax_fine_evt_s,'spline');
    
        %figure; plot(fax_TR_evt_s,RZ_evt,'-b*'); hold on; plot(fax_fine_evt_s, RZ_evt_i,'-r*');

        % the one-point hack (adjust length based on interval size...
        % (if -10:15 make 30, if -5:0 make 6, if 0:10 make 12)
        len_this = size(RZ_evt_i,2); 
        if len_this < 6
            diff = 6-len_this;
            RZ_evt_i = [RZ_evt_i repmat(RZ_evt_i(end),1,diff)];
        end
        % stack events on 3rd dim
         ROI_evts(e,:) = RZ_evt_i;
         
        
    end