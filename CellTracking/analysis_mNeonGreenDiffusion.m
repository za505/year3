%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

%inputs
lengthTraces=0;
intensityTraces=0;
adjTraces=0;
normTraces=0;
compTraces=0;
positionCheck=0;
modelFit=1;
modelCheck=1;
modelInterp=1;
%% load the data
dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/10312021_analysis';
cd([dirsave '/MatFiles'])

datadir1=dir(['10232021_Exp1' '*dm.mat']); %LB, frame rate=1 min, replicate 1
datadir2=dir(['10232021_Exp2' '*dm.mat']); %PBS, frame rate=1 min, replicate 1
datadir3=dir(['10262021_Exp1' '*dm.mat']); %LB, frame rate=1 min, replicate 2
datadir4=dir(['10262021_Exp2' '*dm.mat']); %PBS, frame rate=1 min, replicate 2
datadir5a=dir(['10282021_Exp1' '*dm.mat']); %LB, frame rate=5 min, adjusted frame rate
datadir5b=dir(['10282021_Exp1' '*dm2.mat']); %LB, frame rate=5 min, unadjusted
datadir6=dir(['10302021_Exp1' '*dm.mat']); %LB, frame rate=30 s
datadir7=dir(['10302021_Exp2' '*dm.mat']); %LB, frame rate=15 s

[intensity1, adjintensity1, normintensity1, positions1, lCell1, time1]=dataInput(datadir1);
[intensity2, adjintensity2, normintensity2, positions2, lCell2, time2]=dataInput(datadir2);
[intensity3, adjintensity3, normintensity3, positions3, lCell3, time3]=dataInput(datadir3);
[intensity4, adjintensity4, normintensity4, positions4, lCell4, time4]=dataInput(datadir4);
[intensity5a, adjintensity5a, normintensity5a, positions5a, lCell5a, time5a]=dataInput(datadir5a);
[intensity5b, adjintensity5b, normintensity5b, positions5b, lCell5b, time5b]=dataInput(datadir5b);
[intensity6, adjintensity6, normintensity6, positions6, lCell6, time6]=dataInput(datadir6);
[intensity7, adjintensity7, normintensity7, positions7, lCell7, time7]=dataInput(datadir7);

%% sanity check: length traces
if lengthTraces==1
    
    cd([dirsave '/LengthTraces'])
    [l_avg1, l_std1]=lengthPlot(time1, lCell1,  'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1'); 
    [l_avg2, l_std2]=lengthPlot(time2, lCell2, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2'); 
    [l_avg3, l_std3]=lengthPlot(time3, lCell3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1'); 
    [l_avg4, l_std4]=lengthPlot(time4, lCell4,  'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2'); 
    [l_avg5a, l_std5a]=lengthPlot(time5a, lCell5a, 'LB perfusion, frame rate = 5 min, adjusted', 1, '10282021_Exp1_adj'); 
    [l_avg5b, l_std5b]=lengthPlot(time5b, lCell5b,  'LB perfusion, frame rate = 5 min, unadjusted', 1, '10282021_Exp1_unadj'); 
    [l_avg6, l_std6]=lengthPlot(time6, lCell6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');  
    [l_avg7, l_std7]=lengthPlot(time7, lCell7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2'); 
    
end

%% sanity check: fluorescence traces
if intensityTraces==1
    cd([dirsave '/IntensityTraces'])
    [intensity_avg1, intensity_std1]=intensityPlot(time1, intensity1, 'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1');
    [intensity_avg2, intensity_std2]=intensityPlot(time2, intensity2, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2');
    [intensity_avg3, intensity_std3]=intensityPlot(time3, intensity3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1');
    [intensity_avg4, intensity_std4]=intensityPlot(time4, intensity4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2');
    [intensity_avg5a, intensity_std5a]=intensityPlot(time5a, intensity5a, 'LB perfusion, frame rate = 5 min, adjusted', 1, '10282021_Exp1_adj');
    [intensity_avg5b, intensity_std5b]=intensityPlot(time5b, intensity5b, 'LB perfusion, frame rate = 5 min, unadjusted', 1, '10282021_Exp1_unadj');
    [intensity_avg6, intensity_std6]=intensityPlot(time6, intensity6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');
    [intensity_avg7, intensity_std7]=intensityPlot(time7, intensity7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2');
end

%% sanity check: adjusted fluorescence traces
if adjTraces==1
    cd([dirsave '/AdjTraces'])
    [adjintensity_avg1, adjintensity_std1]=adjintensityPlot(time1, adjintensity1, 'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1');
    [adjintensity_avg2, adjintensity_std2]=adjintensityPlot(time2, adjintensity2,  'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2');
    [adjintensity_avg3, adjintensity_std3]=adjintensityPlot(time3, adjintensity3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1');
    [adjintensity_avg4, adjintensity_std4]=adjintensityPlot(time4, adjintensity4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2');
    [adjintensity_avg5a, adjintensity_std5a]=adjintensityPlot(time5a, adjintensity5a, 'LB perfusion, frame rate = 5 min, adjusted', 1, '10282021_Exp1_adj');
    [adjintensity_avg5b, adjintensity_std5b]=adjintensityPlot(time5b, adjintensity5b, 'LB perfusion, frame rate = 5 min, unadjusted', 1, '10282021_Exp1_unadj');
    [adjintensity_avg6, adjintensity_std6]=adjintensityPlot(time6, adjintensity6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');
    [adjintensity_avg7, adjintensity_std7]=adjintensityPlot(time7, adjintensity7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2');
end

%% normalized fluorescence traces
if normTraces==1
    cd([dirsave '/NormTraces'])
    [normintensity_avg1, normintensity_std1]=normintensityPlot(time1, normintensity1, 'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1');
    [normintensity_avg2, normintensity_std2]=normintensityPlot(time2, normintensity2, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2');
    [normintensity_avg3, normintensity_std3]=normintensityPlot(time3, normintensity3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1');
    [normintensity_avg4, normintensity_std4]=normintensityPlot(time4, normintensity4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2');
    [normintensity_avg5a, normintensity_std5a]=normintensityPlot(time5a, normintensity5a, 'LB perfusion, frame rate = 5 min, adjusted', 1, '10282021_Exp1_adj');
    [normintensity_avg5b, normintensity_std5b]=normintensityPlot(time5b, normintensity5b, 'LB perfusion, frame rate = 5 min, unadjusted', 1, '10282021_Exp1_unadj');
    [normintensity_avg6, normintensity_std6]=normintensityPlot(time6, normintensity6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');
    [normintensity_avg7, normintensity_std7]=normintensityPlot(time7, normintensity7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2');
end

%% compare fluorescence traces
if compTraces==1
    %set the color for the traces
    colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
    
    figure(1), hold on
    for i=1:height(normintensity1)
        plot(time1, normintensity1(i,:), 'Color', colorcode{1})
    end
    for i=1:height(normintensity3)
        plot(time3, normintensity3(i,:), 'Color', colorcode{2})
    end
    for i=1:height(normintensity5a(1:24))
        plot(time5a(1:24), normintensity5a(1:24), 'Color', colorcode{3})
    end
    for i=1:height(normintensity6)
        plot(time6, normintensity6(i,:), 'Color', colorcode{4})
    end
    for i=1:height(normintensity7)
        plot(time7, normintensity7(i,:), 'Color', colorcode{5})
    end
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title('Diffusion in LB')
    text(75, .8, 'frame rate = 1 min, replicate 1');
    text(68, .82, '_____', 'Color', colorcode{1});
    text(75, .75, 'frame rate = 1 min, replicate 2');
    text(68, .77, '_____', 'Color', colorcode{2});
    text(75, .7, 'frame rate = 5 min, adjusted');
    text(68, .72, '_____', 'Color', colorcode{3});
    text(75, .65, 'frame rate = 30 s');
    text(68, .67, '_____', 'Color', colorcode{4});
    text(75, .6, 'frame rate = 15 s');
    text(68, .62, '_____', 'Color', colorcode{5});
    saveas(gcf,'LB_normintensity.png')
    saveas(gcf, 'LB_normintensity.fig')

    figure(2), hold on
    for i=1:height(normintensity2)
        plot(time2, normintensity2(i,:), 'Color', colorcode{1})
    end
    for i=1:height(normintensity4)
        plot(time4, normintensity4(i,:), 'Color', colorcode{2})
    end
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title('Diffusion After Incubation in PBS')
    text(35, .8, 'replicate 1');
    text(32, .82, '_____', 'Color', 'Color', colorcode{1});
    text(35, .75, 'replicate 2');
    text(32, .77, '_____', 'Color', 'Color', colorcode{2});
    saveas(gcf,'PBS_normintensity.png')
    saveas(gcf, 'PBS_normintensity.fig')

end
%% determine whether position affects fluorescence traces
if positionCheck==1
    %set the color for the traces
    colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
    %set the position of the text in the legend
    ty=[.8, .75, .7, .65, .6, .55];

    figure, hold on
    for i=1:height(normintensity1)
        for h=2:length(positions1)
            if i>positions1(h-1)&i<=positions1(h)
                plot(time1, normintensity1(i,:), 'Color', colorcode{h-1})
            end
        end
    end
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title('Diffusion in LB, replicate 1')
    for i=1:length(positions1)-1
        text(time1(end)-10, ty(i), ['position ' num2str(i)]);
        text(time1(end)-13, ty(i)+.02, '_____', 'Color', colorcode{i}, 'FontWeight', 'bold');
    end
    saveas(gcf,'10232021_Exp1_colonyComp.png')
    saveas(gcf, '10232021_Exp1_colonyComp.fig')
    
    figure, hold on
    for i=1:height(normintensity2)
        for h=2:length(positions2)
            if i>positions2(h-1)&i<=positions2(h)
                plot(time2, normintensity2(i,:), 'Color', colorcode{h-1})
            end
        end
    end
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title('Diffusion After PBS Incubation, replicate 1')
    for i=1:length(positions2)-1
        text(time2(end)-10, ty(i), ['position ' num2str(i)]);
        text(time2(end)-13, ty(i)+.02, '_____', 'Color', colorcode{i}, 'FontWeight', 'bold');
    end
    saveas(gcf,'10232021_Exp2_colonyComp.png')
    saveas(gcf, '10232021_Exp2_colonyComp.fig')
    
    figure, hold on
    for i=1:height(normintensity3)
        for h=2:length(positions3)
            if i>positions3(h-1)&i<=positions3(h)
                plot(time3, normintensity3(i,:), 'Color', colorcode{h-1})
            end
        end
    end
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title('Diffusion in LB, replicate 2')
    for i=1:length(positions3)-1
        text(time3(end)-10, ty(i), ['position ' num2str(i)]);
        text(time3(end)-13, ty(i)+.02, '_____', 'Color', colorcode{i}, 'FontWeight', 'bold');
    end
    saveas(gcf,'10262021_Exp1_colonyComp.png')
    saveas(gcf, '10262021_Exp1_colonyComp.fig')
    
    figure, hold on
    for i=1:height(normintensity4)
        for h=2:length(positions4)
            if i>positions4(h-1)&i<=positions4(h)
                plot(time4, normintensity4(i,:), 'Color', colorcode{h-1})
            end
        end
    end
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title('Diffusion After PBS Incubation, replicate 2')
    for i=1:length(positions4)-1
        text(time4(end)-10, ty(i), ['position ' num2str(i)]);
        text(time4(end)-13, ty(i)+.02, '_____', 'Color', colorcode{i}, 'FontWeight', 'bold');
    end
    saveas(gcf,'10262021_Exp2_colonyComp.png')
    saveas(gcf, '10262021_Exp2_colonyComp.fig')
end
%% fit the data to an exponential model
if modelFit==1

    cd([dirsave '/ModelFit'])
    [fit1, tau1, yhat1]=expModel(time1, normintensity1, 'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1');
    [fit2, tau2, yhat2]=expModel(time2, normintensity2,  'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2');
    [fit3, tau3, yhat3]=expModel(time3, normintensity3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1');
    [fit4, tau4, yhat4]=expModel(time4, normintensity4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2');
    [fit5a, tau5a, yhat5a]=expModel(time5a(1:24), normintensity5a(:, 1:24), 'LB perfusion, frame rate = 5 min, adjusted', 1, '10282021_Exp1_adj');
    [fit5b, tau5b, yhat5b]=expModel(time5b, normintensity5b, 'LB perfusion, frame rate = 5 min, unadjusted', 1, '10282021_Exp1_unadj');
    [fit6, tau6, yhat6]=expModel(time6, normintensity6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');
    [fit7, tau7, yhat7]=expModel(time7, normintensity7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2');
    
end
%% determine the quality of the fit
if modelCheck==1
    cd([dirsave '/ModelCheck'])

    [residuals1, est1]=residualPlot(normintensity1, yhat1, time1, ['LB perfusion, frame rate = 1 min, replicate 1'], 1, ['10232021_Exp1']);
   [residuals2, est2]=residualPlot(normintensity2, yhat2, time2, ['PBS 1-hour Incubation, frame rate = 1 min, replicate 1'], 1, ['10232021_Exp2']);
   [residuals3, est3]=residualPlot(normintensity3, yhat3, time3, ['LB perfusion, frame rate = 1 min, replicate 2'], 1, ['10262021_Exp1']);
   [residuals4, est4]=residualPlot(normintensity4, yhat4, time4, ['PBS 1-hour Incubation, frame rate = 1 min, replicate 2'], 1, ['10262021_Exp2']);

   [residuals5a, est5a]=residualPlot(normintensity5a(:, 1:24), yhat5a, time5a(1:24), ['LB perfusion, frame rate = 5 minutes, adjusted'], 1, ['10282021_Exp1_adj']);

   [residuals5b, est5b]=residualPlot(normintensity5b, yhat5b, time5b, ['LB perfusion, frame rate = 5 minutes, unadjusted'], 1, ['10282021_Exp1_unadj']);

   [residuals6, est6]=residualPlot(normintensity6, yhat6, time6, ['LB perfusion, frame rate = 30 s'], 1, ['10302021_Exp1']);

   [residuals7, est7]=residualPlot(normintensity7, yhat7, time7, ['LB perfusion, frame rate = 15 s'], 1, ['10302021_Exp2']);

%     % compare the time constants
%     figure, hold on
%     histogram(tau1, 'BinWidth', 1)
%     histogram(tau3, 'BinWidth', 1)
%     histogram(tau2, 'BinWidth', 1)
%     histogram(tau4, 'BinWidth', 1)
%     legend({'LB, rep1', 'LB, rep2', 'PBS, rep1', 'PBS, rep2'})
%     title('Time Constants')
%     xlabel('\tau (min^{-1})')
%     ylabel('Count')
%     % %saveas(gcf,'timescale_hist.png')
%     % %saveas(gcf, 'timescale_hist.fig')
% 
%     % find the mean and standard deviation of the time constants
%     tau_means = [mean(tau1), mean(tau3), mean(tau2), mean(tau4)];
%     tau_std = [std(tau1), std(tau3), std(tau2), std(tau4)];
% 
%     figure
%     scatter(categorical({'LB, rep1', 'LB, rep2', 'PBS, rep1', 'PBS, rep2'}), tau_means)
end

%% interpolate t_{1/2}
if modelInterp==1
    cd([dirsave '/ModelInterp'])
    [xq1, vq1, thalf1]=interpPlot(normintensity1, time1, 'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1');
    [xq2, vq2, thalf2]=interpPlot(normintensity2, time2, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2');
    [xq3, vq3, thalf3]=interpPlot(normintensity3, time3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1');
    [xq4, vq4, thalf4]=interpPlot(normintensity4, time4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2');
    [xq5a, vq5a, thalf5a]=interpPlot(normintensity5a(:, 1:24), time5a(1:24), 'LB perfusion, frame rate = 5 minutes, adjusted', 1, '10282021_Exp1_adj');
    [xq5b, vq5b, thalf5b]=interpPlot(normintensity5b, time5b, 'LB perfusion, frame rate = 5 minutes, unadjusted', 1, '10282021_Exp1_unadj');
    [xq6, vq6, thalf6]=interpPlot(normintensity6, time6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');
    [xq7, vq7, thalf7]=interpPlot(normintensity7, time7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2');
end

%% plots
cd(dirsave)

groups = categorical({'PBS, rep1', 'PBS, rep2', 'LB, rep1', 'LB, rep2', 'LB, 5 min unadj', 'LB, 5 min adj', 'LB, 30 s', 'LB, 15 s'});
tau = {tau2; tau4; tau1; tau3; tau5a; tau5b; tau6; tau7};
order = [1, 2, 3, 4, 5, 6, 7, 8];

figure, hold on
for i=1:height(tau)
    a=scatter(repelem(order(i), length(tau{i})), tau{i}, 'MarkerEdgeColor', '#0072BD');
    b=errorbar(order(i), mean(tau{i}), -std(tau{i}), std(tau{i}), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
end
xlim([0 9])
ylim([0 Inf])
xticklabels([' ', groups, ' '])
xlabel('Experiment')
ylabel('Tau (min^{-1})')
title('Distribution of \tau Values')
% saveas(gcf, 'tau.png')
% saveas(gcf, 'tau.fig')

thalf = {thalf2; thalf4; thalf1; thalf3; thalf5a; thalf5b; thalf6; thalf7};

figure, hold on
for i=1:height(tau)
    a=scatter(repelem(order(i), length(thalf{i})), thalf{i}, 'MarkerEdgeColor', '#0072BD');
    b=errorbar(order(i), mean(thalf{i}), -std(thalf{i}), std(thalf{i}), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
end
xlim([0 9])
ylim([0 Inf])
xticklabels([' ', groups, ' '])
xlabel('Experiment')
ylabel('t_{1/2} (min^{-1})')
title('Distribution of t_{1/2} Values')
% saveas(gcf, 'thalf.png')
% saveas(gcf, 'thalf.fig')

%% Functions
function [intensity, adjintensity, normintensity, positions, lCell, time]=dataInput(datadir)
    intensity=[];
    adjintensity=[];
    normintensity=[];
    lCell=[];
    positions=[0];
    for i=1:length(datadir)
        load(datadir(i).name, 'icell_intensity', 'adj_intensity', 'norm_intensity', 'lcell', 'time')
        intensity=[intensity; icell_intensity];
        adjintensity=[adjintensity; adj_intensity];
        normintensity=[normintensity; norm_intensity];
        lCell=[lCell; lcell];
        positions(i+1)=height(norm_intensity)+positions(i);
    end
end

function [l_avg, l_std]=lengthPlot(time, lCell, text, save, saveAs)

    l_avg = mean(lCell, 1, 'omitnan');
    l_std = std(lCell, 0, 1, 'omitnan');
    
    figure(1), hold on
    for i=1:height(lCell)
        plot(time, lCell(i,:))    
    end
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    title(text)

    if save==1
        saveas(gcf, [saveAs '_lTraces.png'])
        saveas(gcf, [saveAs '_lTraces.fig'])
    end

    figure(2), hold on
    ciplot((l_avg-l_std),(l_avg+l_std),time,[0.75 0.75 1])
    plot(time,l_avg,'-r')
    ylim([0 Inf])
    xlabel('Time (minutes)')
    ylabel('Length Average (\mum)')
    title(text)
    
    if save==1
        saveas(gcf, [saveAs '_lavgTraces.png'])
        saveas(gcf, [saveAs '_lavgTraces.fig'])
    end

    pause, close all
end

function [intensity_avg, intensity_std]=intensityPlot(time, intensity, text, save, saveAs)
    
    intensity_avg=mean(intensity, 1, 'omitnan');
    intensity_std=std(intensity, 0, 1, 'omitnan');
    
        figure, hold on
        for i=1:height(intensity)
            plot(time, intensity(i, :), 'Color', '#77AC30')
        end
        xlabel('Time (minutes)')
        ylabel('Fluorescence (A.U.)')
        title(text)
        if save==1
            saveas(gcf, [saveAs '_intensityTraces.png'])
            saveas(gcf, [saveAs '_intensityTraces.fig'])
        end
        
        figure(2), hold on
        ciplot((intensity_avg-intensity_std),(intensity_avg+intensity_std),time,[0.75 0.75 1])
        plot(time,intensity_avg,'-r')
        xlabel('Time (minutes)')
        ylabel('Average Fluorescence (A.U.)')
        title(text)
        
        if save==1
            saveas(gcf, [saveAs '_intensityAvg.png'])
            saveas(gcf, [saveAs '_intensityAvg.fig'])
        end

        pause, close all
end

function [adjintensity_avg, adjintensity_std]=adjintensityPlot(time, adjintensity, text, save, saveAs)
    
    adjintensity_avg=mean(adjintensity, 1, 'omitnan');
    adjintensity_std=std(adjintensity, 0, 1, 'omitnan');
    
    figure, hold on
    for i=1:height(adjintensity)
        plot(time, adjintensity(i, :), 'Color', '#77AC30')
    end
    xlabel('Time (minutes)')
    ylabel('Adjusted Fluorescence (A.U.)')
    title(text)
    if save==1
        saveas(gcf, [saveAs '_adjintensityTraces.png'])
        saveas(gcf, [saveAs '_adjintensityTraces.fig'])
    end

    figure(2), hold on
    ciplot((adjintensity_avg-adjintensity_std),(adjintensity_avg+adjintensity_std),time,[0.75 0.75 1])
    plot(time,adjintensity_avg,'-r')
    xlabel('Time (minutes)')
    ylabel('Adjusted Fluorescence (A.U.)')
    title(text)
    
    if save==1
        saveas(gcf, [saveAs '_adjintensityAvg.png'])
        saveas(gcf, [saveAs '_adjintensityAvg.fig'])
    end

    pause, close all
end

function [normintensity_avg, normintensity_std]=normintensityPlot(time, normintensity, text, save, saveAs)
    
    normintensity_avg=mean(normintensity, 1, 'omitnan');
    normintensity_std=std(normintensity, 0, 1, 'omitnan');
    
    figure, hold on
    for i=1:height(normintensity)
        plot(time, normintensity(i, :), 'Color', '#77AC30')
    end
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title(text)
    if save==1
        saveas(gcf, [saveAs '_normintensityTraces.png'])
        saveas(gcf, [saveAs '_normintensityTraces.fig'])
    end

    figure(2), hold on
    ciplot((normintensity_avg-normintensity_std),(normintensity_avg+normintensity_std),time,[0.75 0.75 1])
    plot(time,normintensity_avg,'-r')
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title(text)

    if save==1
        saveas(gcf, [saveAs '_normintensityAvg.png'])
        saveas(gcf, [saveAs '_normintensityAvg.fig'])
    end
        
    pause, close all
end

function [fit, tau, yhat]=expModel(time, normintensity, text, save, saveAs)
  
    modelfun=@(tau,x)exp(-x./tau);
    tau0=0.1;

    fit=[];
    tau=[];
    yhat=normintensity;
    for i=1:height(normintensity)
        if ~isnan(normintensity(i,:))
            fit=[fit, i];
            tau(i)=nlinfit(time, normintensity(i,:), modelfun, tau0);
            yhat(i,:)=modelfun(tau(i), time);   
        end
    end

    figure, hold on
    for i=1:height(normintensity)
        plot(time, normintensity(i,:),'Color', '#77AC30')
        plot(time, yhat(i,:), '--k')
    end
    ylim([0 Inf])
    legend({'data', 'model'})
    title(text)
    xlabel('Time (min)')
    ylabel('Normalized Intensity (A.U.)')
    if save==1
        saveas(gcf, [saveAs '_expFit.png'])
        saveas(gcf, [saveAs '_expFit.fig'])
    end
        
    pause, close all
end

function [residuals, est]=residualPlot(normintensity, yhat, time, text, save, saveAs)

    residuals=abs(normintensity-yhat);
    est=residuals./normintensity;
    est(est==Inf)=0;


    figure, hold on
    for i=1:height(residuals)
        plot(time, residuals(i,:), 'Color', '#7E2F8E')
    end
    xlabel('Time (minutes)')
    ylabel('$\frac{|y-\hat{y}|}{y}$','Interpreter','latex', 'FontSize', 20)
    title(text)

    
    if save==1
        saveas(gcf, [saveAs '_residuals.png'])
        saveas(gcf, [saveAs '_residuals.fig'])
    end

    pause, close all
end

function [xq, vq, thalf]=interpPlot(normintensity,time, text, save, saveAs)

    xq=[0:0.0167:time(end)];
    vq=nan(height(normintensity), length(xq));
    
    for i=1:height(normintensity)
        vq(i,:)=interp1(time,normintensity(i,:),xq);
    end
    
    figure, hold on
    for i=1:height(normintensity)
        plot(xq,vq(i,:),':.', 'Color', '#D95319');
        plot(time,normintensity(i,:),'o', 'Color','#0072BD');
    end
    ylim([0 Inf])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluoresence (A.U.)')
    legend({'interpolation', 'data'})
    title([text ', Linear Interpolation']);
    [~, tidx]=min(abs(vq-0.5), [], 2);
    thalf=xq(tidx);
    
    if save==1
        saveas(gcf, [saveAs '_interp.png'])
        saveas(gcf, [saveAs '_interp.fig'])
    end
    
    close all
end

function ciplot(lower,upper,x,colour,Trans);
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Raymond Reynolds 24/11/06

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<4
    colour='b';
end

if nargin<3
    x=1:length(lower);
end

if nargin<5
    Trans=1;
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

fill([x fliplr(x)],[upper fliplr(lower)],colour,'LineStyle','none','FaceAlpha',Trans)
end