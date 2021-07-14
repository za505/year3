%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all
%re-analysis of 06/16/2021 data
%% first, I want to compare the change in cell wall length in Exp1 and Exp2
%the change when switching from PBS + 5% NLS to *PBS [+ Mg2]
exp1=['/Users/zarina/Documents/MATLAB/MatlabReady/06162021_analysis/06162021_Exp1'];
exp2=['/Users/zarina/Documents/MATLAB/MatlabReady/06162021_analysis/06162021_Exp2'];

cd(exp1)
load(['06162021_Exp1_BTphase.mat'], 'lcell');
exp1_lcell=lcell;

cd(exp2)
load(['06162021_Exp2_BTphase.mat'], 'lcell', 'time');
exp2_lcell=lcell;

ti=find(time==420); %time initial 
tf=find(time==540); %time final

d1=((exp1_lcell(:, tf)-exp1_lcell(:, ti))./exp1_lcell(:, ti))*100;
d2=((exp2_lcell(:, tf)-exp2_lcell(:, ti))./exp2_lcell(:, ti))*100;

h1=histogram(d1); hold on
h2=histogram(d2);

%h1.Normalization = 'probability';
h1.BinWidth = 0.5;
%h1.NumBins = 10;
%h2.Normalization = 'probability';
h2.BinWidth = 0.5;
%h2.NumBins = 10;
title('Percent Change in Cell Wall Length')
xlabel('Percent Change (%)')
legend({'No Mg^{2+}', 'Mg^{2+}'})

cd('..')
saveas(gcf, 'cwlHist.png')
saveas(gcf, 'cwlHist.fig')

close
% %% Function
% function [d1]=lengthChange=(l1,l2)
% 
%     d1=l1-l2;
%     d1=d1/l2;
% 
% end

%% given that the pixels are saturated in these movies until a certain frame, I want to edit those plots
cd(exp1)
load(['06162021_Exp1_BTfluo.mat'], 'unsatFrame', 'icell', 'ratio', 'xtime', 'labels', 'time', 'ncells');
exp1_unsatFrame=unsatFrame;
exp1_icell=icell;
exp1_ratio=ratio;
exp1_xtime=xtime;
exp1_labels=labels;
exp1_time=time;
exp1_ncells=ncells;

cd(exp2)
load(['06162021_Exp2_BTfluo.mat'], 'unsatFrame', 'icell', 'ratio', 'xtime', 'labels', 'time', 'ncells');
exp2_unsatFrame=unsatFrame;
exp2_icell=icell;
exp2_ratio=ratio;
exp2_xtime=xtime;
exp2_labels=labels;
exp2_time=time;
exp2_ncells=ncells;

%Plot data
cd('/Users/zarina/Downloads/NYU/Year3_2021_Summer/Updates/GroupMeeting')

figure, hold on, 
for i=1:exp1_ncells
    plot(exp1_time(exp1_unsatFrame:end), exp1_icell{1}(i,exp1_unsatFrame:end))
end
xlabel('Time (h)')
ylabel('Intensity (A.U.)')
%xlim([-0.2 Inf])
fig2pretty
for i=1:length(exp1_labels)
    if exp1_xtime(i)>=exp1_time(exp1_unsatFrame)
        xline(exp1_xtime(i), '--k', exp1_labels(i)) 
    end
end
saveas(gcf, ['06162021_Exp1_intensity.fig'])
saveas(gcf, ['06162021_Exp1_intensity.png'])

figure, hold on, 
for i=1:exp1_ncells
    plot(exp1_time(exp1_unsatFrame:end), exp1_ratio{1}(i,exp1_unsatFrame:end))
end
xlabel('Time (h)')
ylabel('Cellular Intensity/Background Intensity (A.U.)')
%xlim([-0.2 Inf])
fig2pretty
for i=1:length(exp1_labels)
    if exp1_xtime(i)>=exp1_time(exp1_unsatFrame)
        xline(exp1_xtime(i), '--k', exp1_labels(i)) 
    end
end
saveas(gcf, ['06162021_Exp1_ratio.fig'])
saveas(gcf, ['06162021_Exp1_ratio.png'])

figure, hold on, 
for i=1:exp2_ncells
    plot(exp2_time(exp2_unsatFrame:end), exp2_icell{1}(i,exp2_unsatFrame:end))
end
xlabel('Time (h)')
ylabel('Intensity (A.U.)')
%xlim([-0.2 Inf])
fig2pretty
for i=1:length(exp2_labels)
    if exp2_xtime(i)>=exp2_time(exp2_unsatFrame)
        xline(exp2_xtime(i), '--k', exp2_labels(i)) 
    end
end
saveas(gcf, ['06162021_Exp2_intensity.fig'])
saveas(gcf, ['06162021_Exp2_intensity.png'])

figure, hold on, 
for i=1:exp2_ncells
    plot(exp2_time(exp2_unsatFrame:end), exp2_ratio{1}(i,exp2_unsatFrame:end))
end
xlabel('Time (h)')
ylabel('Cellular Intensity/Background Intensity (A.U.)')
%xlim([-0.2 Inf])
fig2pretty
for i=1:length(exp2_labels)
    if exp2_xtime(i)>=exp2_time(exp2_unsatFrame)
        xline(exp2_xtime(i), '--k', exp2_labels(i)) 
    end
end
saveas(gcf, ['06162021_Exp2_ratio.fig'])
saveas(gcf, ['06162021_Exp2_ratio.png'])