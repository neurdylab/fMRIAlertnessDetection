%% Remake Main Paper Figures
% 1.23.2021
% Sarah Goodale - author
% ELIFE paper: fMRI-based detection of alertness predicts electrophysiological and behavioral variability

% citation: 


%% Figure 3B -------------------------------------------------------------- 
% Removes the columns that are 0 bc subjects did not have any misses.
dims = size(allsub_missavg);
idx = isnan(allsub_missavg(:,1));
allsub_missavg_new = allsub_missavg;
allsub_missavg_new(idx,:) = [];

figure;
semshading_ribbonplot(allsub_hitavg,.3,'',ROI_evts_ts); hold on;
semshading_ribbonplot(allsub_missavg_new,.3,'b',ROI_evts_ts);
title('Mean and Std Error of All subject ROI Hit and Miss');

%% Figure 3C and 4A and supp 4-4 ------------------------------------------
% Plots individual hit and misses in dot plot switch allsub_missavg and
% allsub_hitavg to slow and fast hits respectively.
allsub_missavg(allsub_missavg==0) = NaN; % replace 0s or misses to NaN
avg_miss = nanmean(allsub_missavg,2); % avg across
avg_hit = mean(allsub_hitavg,2);

rts = [avg_hit avg_miss];
x = [.25 .75];
figure;
for i = 1:length(avg_hit)
  
   plot(x,rts(i,:),'-o','LineWidth',2,'Color',[0,0,0],'MarkerSize',14,'MarkerFaceColor',[0.30,0.75,0.93],'MarkerEdgeColor',[0.00,0.45,0.74]); hold on;  
  
   if ~isnan(allsub_missavg(i,:))
   [h,p(i)] = ttest2(allsub_hitavg(i,:),allsub_missavg(i,:));
   end
   
end

% set x-lims for display
xmin = 0;
xmax = 1;
xlim([xmin xmax]);

% choose y-ticks
set(gca,'ytick',[-0.4:0.1:0.2]); % 3C
ylim([-0.4,0.2]);  
% set(gca,'ytick',[-0.25:0.05:0.25]); % 4A
% ylim([-0.25,0.25]);
set(gca,'TickDir','out');
set(gca,'TickLength',[0.015,0.015]);

% numbers slightly bigger and gray
set(gca,'fontsize',16);
set(gca,'YColor',0.25*[1 1 1]);

% x-axis invisible
set(gca,'Xcolor',[1 1 1]);
box off; % right-side y-axis off

%% Figure 4B/C and supp 4-1, 4-3 change mean and std err vectors accordingly

%example
allsub_missavg(allsub_missavg==0) = NaN;
allsub_fasthitavg(allsub_fasthitavg==0) = NaN;
all_avg_miss = nanmean(allsub_missavg,'all'); 
std_miss = std(nanmean(allsub_missavg,2),'omitnan');
all_avg_fast = nanmean(allsub_fasthitavg,'all'); std_fast = std(nanmean(allsub_fasthitavg,2),'omitnan');
all_avg_slow = mean(allsub_slowhitavg,'all'); std_slow = std(nanmean(allsub_fasthitavg,2),'omitnan');

meanvec = [all_avg_fast all_avg_slow all_avg_miss];
errvec = [std_fast/sqrt(12) std_slow/sqrt(12) std_miss/sqrt(12)];


% meanvec = [0.028101492	-0.050842697	-0.130593206]; 
% errvec = [0.016723986	0.023996924	0.03062143];

figure;
set(gcf,'color','w')

hh = bar([1,1.25,1.6],diag(meanvec),'stacked');
hold on;
ee = errorbar([1,1.25,1.6],meanvec,errvec,'.', ...
              'color',0.3*[1 1 1],'linewidth',1);

% change bar colors
hh(1).FaceColor = [0.85,0.51,0.10];
hh(2).FaceColor = [0.97,0.76,0.27];
hh(3).FaceColor = [0.41,0.60,0.79];
% remove bar outline and shrink width
for b=1:3
    hh(b).EdgeColor = [1 1 1];
    hh(b).BarWidth=0.65;
end
% set x-lims for display
xmin = 0.8;
xmax = 1.9;
xlim([xmin xmax]);
axis square;

% zero-line gray and thicker
hL = line([xmin xmax],[0 0]);
hL.Color = 0.9*[1 1 1];
hL.LineWidth = 1.5;

% fMRI choose y-ticks
set(gca,'ytick',[-0.18:0.04:0.14]);
ylim([-0.18,0.14]);  
set(gca,'TickDir','out');
set(gca,'TickLength',[0.015,0.015]);

% % EEG choose y-ticks
% set(gca,'ytick',[0:0.3:1.8]);
% ylim([0,1.8]);  
% set(gca,'TickDir','out');
% set(gca,'TickLength',[0.015,0.015]);

% numbers slightly bigger and gray
set(gca,'fontsize',14);
set(gca,'YColor',0.5*[1 1 1]);

% x-axis invisible
set(gca,'Xcolor',[1 1 1]);
box off; % right-side y-axis off

% size
set(gcf,'pos',[ 431   482   375   306])

%% Cross-Correlation Plots - figure 2C and supp figures 2-1,2-3,and 4-4
figure;
lags_TR = TR.*lags; % x-axis is set to seconds by multiplying by TR
semshading_ribbonplot(all_xcorr',.3,'',lags_TR); hold on; % all_xcorr of cross correlation method you're testing
semshading_ribbonplot(normal_all_xcorr',.3,'b',lags_TR); % all_xcorr of cross correlation method of original eeg and vig_pred
ylim([-0.3 .6]); % adjust ylim based on plot
set(gca,'ytick',[-.3:0.3:0.6]); 
set(gca,'fontsize',15);
set(gca,'YColor',0.35*[1 1 1]);
box off;

set(gca,'Xcolor',0.35*[1 1 1]);
set(gca,'xtick',[-50:10:50]);
ylabel('cross-correlation'); xlabel('time lag (s)');