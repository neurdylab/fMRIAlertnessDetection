function semshading_ribbonplot(Y,a,c,lags)

if exist('a')==0 || isempty(c) % a sets the color
    c='r'; 
end
if exist('lags')==0 || isempty(lags) % lags is the x-axis
    lags=1:size(Y,2);
end

Ymean = nanmean(Y,1); % mean over first dimension
Ysem = nanstd(Y,[],1)/sqrt(size(Y,1)); % SEM shading

if exist('a')==0 || isempty(a) % plotting the shading SEM
    fill([lags fliplr(lags)],[Ymean+Ysem fliplr(Ymean-Ysem)],c,'linestyle','none');
    c='k';
else
    fill([lags fliplr(lags)],[Ymean+Ysem fliplr(Ymean-Ysem)],c, 'FaceAlpha', a,'linestyle','none');
end
hold on; plot(lags,Ymean,c,'linewidth',1.5); % plotting the mean line

end