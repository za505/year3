%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all
%re-analysis of 06/06/2021 data
%% first, I want to plot how intensity changes overtime
%the change when switching from PBS + 5% NLS to *PBS [+ Mg2]
colony1=['/Users/zarina/Documents/MATLAB/MatlabReady/06062021_analysis/06062021_Exp1_colony1'];
colony2=['/Users/zarina/Documents/MATLAB/MatlabReady/06062021_analysis/06062021_Exp1_colony2'];

cd(colony1)
load('06062021_Exp1_dm.mat', 'icell_green', 'icell_mcherry', 'nGreen', 'nCherry')
colony1_green=icell_green;
colony1_mcherry=icell_mcherry;

cd(colony2)
load('06062021_Exp1_dm.mat', 'icell_green', 'icell_mcherry', 'nGreen', 'nCherry', 'time', 'tpt', 'tidx')
colony2_green=icell_green;
colony2_mcherry=icell_mcherry;

combo_green=[colony1_green; colony2_green];
combo_mcherry=[colony1_mcherry; colony2_mcherry];

n_green=height(combo_green);
n_mcherry=height(combo_mcherry);

expGreen=cell(n_green,1);
expCherry=cell(n_mcherry,1);

time=time(tidx);
cutoff=find(time==tpt);
time2=time(cutoff:end)-tpt;

norm_green=[];
norm_mcherry=[];

cd('..')

for i=1:n_green
    i
    norm_green=[norm_green;combo_green(i,cutoff:end)/combo_green(i,cutoff)];
    [xData, yData]=prepareCurveData(time2, norm_green(i,:));
    expGreen{i,1}=fit(xData, yData, 'exp1');
    plot(expGreen{i,1}, xData, yData)
    title(['Fitting mNeonGreen Fluorescence Intensity, n= ' num2str(i)]) 
    xlabel('Time (s)')
    ylabel('Normalized Intensity (A.U.)')

    saveas(gcf, ['06062021_Exp1_mNeonGreenFit_' num2str(i) '.png'])
    saveas(gcf, ['06062021_Exp1_mNeonGreenFit_' num2str(i) '.fig'])
    
    close
end

for i=1:n_mcherry
    i
    norm_mcherry=[norm_mcherry; combo_mcherry(i,cutoff:end)./combo_mcherry(i,cutoff)];
    [xData, yData]=prepareCurveData(time2, norm_mcherry(i,:));
    expCherry{i,1}=fit(xData, yData, 'exp1');
    plot(expCherry{i,1}, xData, yData)
    title(['Fitting mCherry Fluorescence Intensity, n= ' num2str(i)]) 
    xlabel('Time (s)')
    ylabel('Normalized Intensity (A.U.)')
    
    saveas(gcf, ['06062021_Exp1_mCherryFit_' num2str(i) '.png'])
    saveas(gcf, ['06062021_Exp1_mCherryFit_' num2str(i) '.fig'])
    
    close
end

poorfit_green=[2,7,10,16,18,19];
poorfit_mcherry=[3,7,9,10,11,12];

tconst_green=[];
tconst_mcherry=[];

figure,hold on
for i=1:n_green
    if ismember(i, poorfit_green)==0
        plot(time2./60, norm_green(i,:))
        %tconst_green=[tconst_green -1/expGreen{i}.b];
    end
end
title('mNeonGreen Fluorescence Intensity Traces') 
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')  
saveas(gcf, ['06062021_Exp1_mNeonGreenTraces.png'])
saveas(gcf, ['06062021_Exp1_mNeonGreenTraces.fig'])
pause, close

figure,hold on
for i=1:n_mcherry
    if ismember(i, poorfit_mcherry)==0
        plot(time2./60, norm_mcherry(i,:)) 
        %tconst_mcherry=[tconst_mcherry -1/expCherry{i}.b];
    end
end
title('mCherry Fluorescence Intensity Traces') 
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')  
saveas(gcf, ['06062021_Exp1_mCherryTraces.png'])
saveas(gcf, ['06062021_Exp1_mCherryTraces.fig'])
pause, close

histogram(tconst_mcherry./60, length(tconst_mcherry))
title('time constant (tau) for mCherry')
xlabel('tau (min-)')
saveas(gcf, ['06062021_Exp1_mCherryHistogram.png'])
saveas(gcf, ['06062021_Exp1_mCherryHistogram.fig'])
pause, close

histogram(tconst_green./60, length(tconst_green))
title('time constant (tau) for mNeonGreen')
xlabel('tau (min-)')
saveas(gcf, ['06062021_Exp1_greenHistogram.png'])
saveas(gcf, ['06062021_Exp1_greenHistogram.fig'])
pause, close
