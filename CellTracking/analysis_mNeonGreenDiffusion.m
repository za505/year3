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
modelFit=1;
modelCheck=0;
modelInterp=0;
tauPlot=0;
bleachCorrect=1;

colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
%% load the data
dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/11032021_analysis';
cd([dirsave '/MatFiles'])

datadir1=dir(['10232021_Exp1' '*dm.mat']); %LB, frame rate=1 min, replicate 1
datadir2=dir(['10232021_Exp2' '*dm.mat']); %PBS, frame rate=1 min, replicate 1
datadir3=dir(['10262021_Exp1' '*dm.mat']); %LB, frame rate=1 min, replicate 2
datadir4=dir(['10262021_Exp2' '*dm.mat']); %PBS, frame rate=1 min, replicate 2
datadir5a=dir(['10282021_Exp1' '*dm.mat']); %LB, frame rate=5 min, adjusted frame rate
datadir5b=dir(['10282021_Exp1' '*dm2.mat']); %LB, frame rate=5 min, unadjusted
datadir6=dir(['10302021_Exp1' '*dm.mat']); %LB, frame rate=30 s
datadir7=dir(['10302021_Exp2' '*dm.mat']); %LB, frame rate=15 s

[intensity1, adjintensity1, normintensity1, positions1, lCell1, time1, t1, lidx1]=dataInput(datadir1);
[intensity2, adjintensity2, normintensity2, positions2, lCell2, time2, t2, lidx2]=dataInput(datadir2);
[intensity3, adjintensity3, normintensity3, positions3, lCell3, time3, t3, lidx3]=dataInput(datadir3);
[intensity4, adjintensity4, normintensity4, positions4, lCell4, time4, t4, lidx4]=dataInput(datadir4);
[intensity5a, adjintensity5a, normintensity5a, positions5a, lCell5a, time5a, t5a, lidx5a]=dataInput(datadir5a);
[intensity5b, adjintensity5b, normintensity5b, positions5b, lCell5b, time5b, t5b, lidx5b]=dataInput(datadir5b);
[intensity6, adjintensity6, normintensity6, positions6, lCell6, time6, t6, lidx6]=dataInput(datadir6);
[intensity7, adjintensity7, normintensity7, positions7, lCell7, time7, t7, lidx7]=dataInput(datadir7);

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


%         figure, hold on
%         for i=1:height(intensity1)
%             plot(time1, intensity1(i, :), 'Color', '#77AC30')
%         end
%         for i=1:height(intensity3)
%             plot(time3, intensity3(i, :), 'Color', '#77AC30')
%         end
%         ylim([0 Inf])
%         xlabel('Time (minutes)')
%         ylabel('Fluorescence (A.U.)')
%         
%                 figure, hold on
%         for i=1:height(intensity2)
%             plot(time2, intensity2(i, :), 'Color', colorcode{2})
%         end
%         for i=1:height(intensity3)
%             plot(time4, intensity4(i, :), 'Color', colorcode{2})
%         end
%         ylim([0 Inf])
%         xlabel('Time (minutes)')
%         ylabel('Fluorescence (A.U.)')

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
    [normintensity_avg1, normintensity_std1]=normintensityPlot(t1, normintensity1, 'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1');
    [normintensity_avg2, normintensity_std2]=normintensityPlot(t2, normintensity2, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2');
    [normintensity_avg3, normintensity_std3]=normintensityPlot(t3, normintensity3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1');
    [normintensity_avg4, normintensity_std4]=normintensityPlot(t4, normintensity4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2');
    [normintensity_avg5a, normintensity_std5a]=normintensityPlot(t5a, normintensity5a, 'LB perfusion, frame rate = 5 min, adjusted', 1, '10282021_Exp1_adj');
    [normintensity_avg5b, normintensity_std5b]=normintensityPlot(t5b, normintensity5b, 'LB perfusion, frame rate = 5 min, unadjusted', 1, '10282021_Exp1_unadj');
    [normintensity_avg6, normintensity_std6]=normintensityPlot(t6, normintensity6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');
    [normintensity_avg7, normintensity_std7]=normintensityPlot(t7, normintensity7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2');
end

%% fit the data to an exponential model
if modelFit==1

    cd([dirsave '/ModelFit'])
    [fit1, tau1, yhat1]=expModel(t1, normintensity1, 'LB perfusion, frame rate = 1 min, replicate 1', 0, '10232021_Exp1');
    [fit2, tau2, yhat2]=expModel(t2, normintensity2,  'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 0, '10232021_Exp2');
    [fit3, tau3, yhat3]=expModel(t3, normintensity3, 'LB perfusion, frame rate = 1 min, replicate 2', 0, '10262021_Exp1');
    [fit4, tau4, yhat4]=expModel(t4, normintensity4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 0, '10262021_Exp2');
%     [fit5a, tau5a, yhat5a]=expModel(t5a(1:17), normintensity5a(:, 1:17), 'LB perfusion, frame rate = 5 min, adjusted', 0, '10282021_Exp1_adj');
%     [fit5b, tau5b, yhat5b]=expModel(t5b, normintensity5b, 'LB perfusion, frame rate = 5 min, unadjusted', 0, '10282021_Exp1_unadj');
%     [fit6, tau6, yhat6]=expModel(t6, normintensity6, 'LB perfusion, frame rate = 30 s', 0, '10302021_Exp1');
%     [fit7, tau7, yhat7]=expModel(t7, normintensity7, 'LB perfusion, frame rate = 15 s', 0, '10302021_Exp2');
    
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
%% determine the quality of the fit
if modelCheck==1
    
   cd([dirsave '/ModelCheck'])
   [residuals1, est1]=residualPlot(normintensity1, yhat1, t1, fit1, 'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1');
   [residuals2, est2]=residualPlot(normintensity2, yhat2, t2, fit2, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2');
   [residuals3, est3]=residualPlot(normintensity3, yhat3, t3, fit3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1');
   [residuals4, est4]=residualPlot(normintensity4, yhat4, t4, fit4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2');
   [residuals5a, est5a]=residualPlot(normintensity5a(:, 1:17), yhat5a, t5a(1:17), fit5a, 'LB perfusion, frame rate = 5 minutes, adjusted', 1, '10282021_Exp1_adj');
   [residuals5b, est5b]=residualPlot(normintensity5b, yhat5b, t5b, fit5b, 'LB perfusion, frame rate = 5 minutes, unadjusted', 1, '10282021_Exp1_unadj');
   [residuals6, est6]=residualPlot(normintensity6, yhat6, t6, fit6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');
   [residuals7, est7]=residualPlot(normintensity7, yhat7, t7, fit7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2');

end

%% interpolate t_{1/2}
if modelInterp==1
    cd([dirsave '/ModelInterp'])
    [xq1, vq1, thalf1]=interpPlot(normintensity1, t1, 'LB perfusion, frame rate = 1 min, replicate 1', 1, '10232021_Exp1');
    [xq2, vq2, thalf2]=interpPlot(normintensity2, t2, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 1', 1, '10232021_Exp2');
    [xq3, vq3, thalf3]=interpPlot(normintensity3, t3, 'LB perfusion, frame rate = 1 min, replicate 2', 1, '10262021_Exp1');
    [xq4, vq4, thalf4]=interpPlot(normintensity4, t4, 'PBS 1-hour Incubation, frame rate = 1 min, replicate 2', 1, '10262021_Exp2');
    [xq5a, vq5a, thalf5a]=interpPlot(normintensity5a(:, 1:17), t5a(1:17), 'LB perfusion, frame rate = 5 minutes, adjusted', 1, '10282021_Exp1_adj');
    [xq5b, vq5b, thalf5b]=interpPlot(normintensity5b, t5b, 'LB perfusion, frame rate = 5 minutes, unadjusted', 1, '10282021_Exp1_unadj');
    [xq6, vq6, thalf6]=interpPlot(normintensity6, t6, 'LB perfusion, frame rate = 30 s', 1, '10302021_Exp1');
    [xq7, vq7, thalf7]=interpPlot(normintensity7, t7, 'LB perfusion, frame rate = 15 s', 1, '10302021_Exp2');
end

%%  plot tau and thalf as a function of frame rate
if tauPlot==1
    cd(dirsave)
    
    %plot tau as a function of initial adjusted fluor. (fluor - background)
    figure, hold on
    scatter([adjintensity1(:,1); adjintensity3(:,1)], [tau1, tau3], 'MarkerFaceColor', colorcode{1})
    scatter(adjintensity6(:,1), tau6, 'MarkerFaceColor', colorcode{2})
    scatter(adjintensity7(:,1), tau7, 'MarkerFaceColor', colorcode{3})
    scatter([adjintensity2(:,1); adjintensity4(:,1)], [tau2, tau4], 'MarkerFaceColor', colorcode{4})
    xlabel('Initial Adjusted Intensity (A.U.)')
    ylabel('Tau (min^{-1})')
    legend({'LB, frame rate = 1 min', 'LB, frame rate = 30 s', 'LB, frame rate = 15 s', 'PBS, frame rate = 1 min'})
    saveas(gcf, 'tau_initialF.png')
    saveas(gcf, 'tau_initialF.fig')    


%     groups = categorical({'PBS, rep1', 'PBS, rep2', 'LB, rep1', 'LB, rep2', 'LB, 5 min unadj', 'LB, 5 min adj', 'LB, 30 s', 'LB, 15 s'});
%     tau = {tau2; tau4; tau1; tau3; tau5a; tau5b; tau6; tau7};
%     order = [1, 2, 3, 4, 5, 6, 7, 8];
% 
%     figure, hold on
%     for i=1:height(tau)
%         a=scatter(repelem(order(i), length(tau{i})), tau{i}, 'MarkerEdgeColor', '#0072BD');
%         b=errorbar(order(i), mean(tau{i}), -std(tau{i}), std(tau{i}), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
%     end
%     xlim([0 9])
%     ylim([0 Inf])
%     xticklabels([' ', groups, ' '])
%     xlabel('Experiment')
%     ylabel('Tau (min^{-1})')
%     title('Distribution of \tau Values')
%     % saveas(gcf, 'tau.png')
%     % saveas(gcf, 'tau.fig')
% 
%     thalf = {thalf2; thalf4; thalf1; thalf3; thalf5a; thalf5b; thalf6; thalf7};
% 
%     figure, hold on
%     for i=1:height(tau)
%         a=scatter(repelem(order(i), length(thalf{i})), thalf{i}, 'MarkerEdgeColor', '#0072BD');
%         b=errorbar(order(i), mean(thalf{i}), -std(thalf{i}), std(thalf{i}), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
%     end
%     xlim([0 9])
%     ylim([0 Inf])
%     xticklabels([' ', groups, ' '])
%     xlabel('Experiment')
%     ylabel('t_{1/2} (min)')
%     title('Distribution of t_{1/2} Values')
%     % saveas(gcf, 'thalf.png')
%     % saveas(gcf, 'thalf.fig')

    %there are three frames rates tested for LB perfusion:1 min, 30 s, and 15s
    frameRate = [repelem(0.25, length(tau7)), repelem(0.5, length(tau6)), repelem(1, length([tau1, tau3]))];
    frameRate_tau = [tau7, tau6, tau1, tau3];
    frameRate_thalf = [thalf7, thalf6, thalf1, thalf3];

    linearCoef_tau = polyfit(frameRate,frameRate_tau,1);
    linearFit_tau = polyval(linearCoef_tau,[0 frameRate]);

    linearCoef_thalf = polyfit(frameRate,frameRate_thalf,1);
    linearFit_thalf = polyval(linearCoef_thalf,[0 frameRate]);
% 
% %     figure, hold on
% %     for i=1:length(frameRate_tau)
% %         a=scatter(frameRate, frameRate_tau, 'MarkerEdgeColor', '#0072BD');
% %         b=errorbar(0.25, mean(tau7), -std(tau7), std(tau7), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
% %         c=errorbar(0.5, mean(tau6), -std(tau6), std(tau6), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
% %         d=errorbar(1, mean([tau1, tau3]), -std([tau1, tau3]), std([tau1, tau3]), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
% %     end
% %     plot([0 frameRate], linearFit_tau, '-r')
% %     ylim([0 35])
% %     xlim([0, 1.1])
% %     xlabel('Frame Rate (min)')
% %     ylabel('\tau (min^{-1})')
% %     title('\tau as a function of frame rate')
% %     % saveas(gcf, 'frameRate_tau.png')
% %     % saveas(gcf, 'frameRate_tau.fig')
% % 
% %     figure, hold on
% %     for i=1:length(frameRate_thalf)
% %         a=scatter(frameRate, frameRate_thalf, 'MarkerEdgeColor', '#0072BD');
% %         b=errorbar(0.25, mean(thalf7), -std(thalf7), std(thalf7), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
% %         c=errorbar(0.5, mean(thalf6), -std(thalf6), std(thalf6), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
% %         d=errorbar(1, mean([thalf1, thalf3]), -std([thalf1, thalf3]), std([thalf1, thalf3]), 'Color', 'black', 'Marker', 'x', 'LineWidth', 1, 'MarkerSize', 8);
% %     end
% %     plot([0 frameRate], linearFit_thalf, '-r')
% %     ylim([0 30])
% %     xlim([0, 1.1])
% %     xlabel('Frame Rate (min)')
% %     ylabel('t_{1/2} (min)')
% %     title('t_{1/2} as a function of frame rate')
% %     % saveas(gcf, 'frameRate_thalf.png')
%     % saveas(gcf, 'frameRate_thalf.fig')
%     
%     figure
%     scatter([intensity7(:,1)', intensity6(:,1)', intensity1(:,1)', intensity3(:,1)'], framerate_tau)
%     
%     
end
%% Photobleach Correction Model
if bleachCorrect==1
    %alpha = linearCoef_tau(1);
%     [Cnew1, dCP1, unb_frac1, Cbl_exp1, gamma1] = photoCorrect(normintensity1, t1, tau1);
%     [Cnew2, dCP2, unb_frac2, Cbl_exp2, gamma2] = photoCorrect(normintensity2, t2, tau2);
%     [Cnew3, dCP3, unb_frac3, Cbl_exp3, gamma3] = photoCorrect(normintensity3, t3, tau3);
%     [Cnew4, dCP4, unb_frac4, Cbl_exp4, gamma4] = photoCorrect(normintensity4, t4, tau4);
%     
% %     %perform normalization
% %     Cnew2 = Cnew2 - Cnew2(:, end);
% %     Cnew2 = Cnew2./Cnew2(:,1);
% %     
% %     Cnew4 = Cnew4 - Cnew4(:, end);
% %     Cnew4 = Cnew4./Cnew4(:,1);
% %     
% %     Combine datasets 
%     LB_combo = [Cnew1(:, 1:93); Cnew3];
%     PBS_combo = [Cnew2; Cnew4(:, 1:48)];
%     
%     LB_avg = mean(LB_combo, 1, 'omitnan');
%     LB_std = std(LB_combo, 0, 1, 'omitnan');
% 
%     PBS_avg = mean(PBS_combo, 1, 'omitnan');
%     PBS_std = std(PBS_combo, 0, 1, 'omitnan');
% 
%     figure('DefaultAxesFontSize',18), hold on
%     ciplot((LB_avg-LB_std),(LB_avg+LB_std),t3,[0.75 0.75 1])
%     plot(t3(1:51),LB_avg(1,1:51),'Color', [0.00 0.45 0.74], 'LineWidth', 2)
%     ciplot((PBS_avg-PBS_std),(PBS_avg+PBS_std),t2,[0.99 0.76 0.93])
%     plot(t2,PBS_avg,'Color', [0.64 0.08 0.18], 'LineWidth', 2)
%     ylim([0 1.3])
%     xlabel('Time (minutes)')
%     ylabel('Normalized Fluorescence (A.U.)')
%     legend({ 'LB Standard Deviation', 'LB Average','PBS Standard Deviation', 'PBS Average'})
%legend({ 'LB Average', 'PBS Average'})

%     figure, hold on
%     for n=1:height(normintensity1)
%         plot(t1, normintensity1(n,:), 'Color', '#FF80CC')
%         plot(t1, Cnew1(n, :), 'Color', '#A2142F')
%     end
%     for n=1:height(normintensity3)
%         plot(t3, normintensity3(n,:), 'Color', '#FF80CC')
%         plot(t3, Cnew3(n, :), 'Color', '#A2142F')
%     end
    
%     figure, hold on
%     for n=1:height(normintensity2)
%         plot(t2, normintensity2(n,:), 'Color', colorcode{6})
%         plot(t2, Cnew2(n, :), 'Color', colorcode{1})
%     end
%     for n=1:height(normintensity4)
%         plot(t4, normintensity4(n,:), 'Color', colorcode{6})
%         plot(t4, Cnew4(n, :), 'Color', colorcode{1})
%     end
    
    
%     figure, hold on
%     ciplot((mean(Cnew1, 1, 'omitnan')-std(Cnew1, 0, 1, 'omitnan')),(mean(Cnew1, 1, 'omitnan')+std(Cnew1, 0, 1, 'omitnan')),t1,[0.75 0.75 1])
%     plot(t1, mean(Cnew1, 1, 'omitnan'), 'Color', '-r')
%     ylim([0 Inf])
%     xlabel('Time (minutes)')
%     ylabel('Normalized Fluorescence (A.U.)')
    
%     figure, hold on
%     ciplot((mean(Cnew2, 1, 'omitnan')-std(Cnew2, 0, 1, 'omitnan')),(mean(Cnew2, 1, 'omitnan')+std(Cnew2, 0, 1, 'omitnan')),t2,[0.75 0.75 1])
%     plot(t2, mean(Cnew2, 1, 'omitnan'), '-r')
%     ylim([0 Inf])
%     xlabel('Time (minutes)')
%     ylabel('Normalized Fluorescence (A.U.)')

%    Combine datasets 
      LB_combo = [normintensity1(:, 1:48); normintensity3(:,1:48)];
      [Cnew1, dCP1, unb_frac1, Cbl_exp1, gamma1] = photoCorrect(LB_combo, t3(1:48), [tau1, tau3], 1);
            
      PBS_combo = [normintensity2; normintensity4(:, 1:48)];
      [CnewA, dCPA, unb_fracA, Cbl_expA, gammaA] = photoCorrect(PBS_combo, t2, [tau2, tau4], 1);
      [CnewB, dCPB, unb_fracB, Cbl_expB, gammaB] = photoCorrect(PBS_combo, t2, [tau2, tau4], 2);
      [CnewC, dCPC, unb_fracC, Cbl_expC, gammaC] = photoCorrect(PBS_combo, t2, [tau2, tau4], 3);
      
    figure, hold on
    %ciplot((mean(CnewA, 1, 'omitnan')-std(CnewA, 0, 1, 'omitnan')),(mean(CnewA, 1, 'omitnan')+std(CnewA, 0, 1, 'omitnan')),t2,[0.99 0.76 0.93])
    %plot(t2, mean(CnewA, 1, 'omitnan'), 'Color', [0.64 0.08 0.18], 'LineWidth', 2)
    %ciplot((mean(CnewB, 1, 'omitnan')-std(CnewB, 0, 1, 'omitnan')),(mean(CnewB, 1, 'omitnan')+std(CnewB, 0, 1, 'omitnan')),t2,[0.75 0.75 1])
    %plot(t2, mean(CnewB, 1, 'omitnan'), 'Color', [0.00 0.45 0.74], 'LineWidth', 2)
    %ciplot((mean(CnewC, 1, 'omitnan')-std(CnewC, 0, 1, 'omitnan')),(mean(CnewC, 1, 'omitnan')+std(CnewC, 0, 1, 'omitnan')),t2,[0.6471 0.9020 0.3137])
    plot(t3(1:48), mean(Cnew1, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 2)
    plot(t2, mean(CnewC, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 2)
    ylim([0 1.3])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    legend({'Average Control', 'Average PBS-incubated'})
    %legend({ 'Original Correction Standard Deviation', 'Original Correction','Average Correction Standard Deviation', 'Average Correction','Forward Correction Standard Deviation', 'Forward Correction'})
    %legend({'Original Correction','Average Correction','Forward Correction'})
end

%% Functions
function [intensity, adjintensity, normintensity, positions, lCell, time, t, lidx]=dataInput(datadir)
    
    %to determine the initial lysis frame, go to the 647 image stack of
    %colony 1, which is usually to the right-most part of the chip
    cd(datadir(1).folder)
    load(datadir(1).name, 'channels')

    for t=5:9

        %to determine the lysis frame, go to the 647 image stack
        new_channel = strrep(channels,'mNeonGreen','647');
        cd(new_channel{1}); 
        fluo_directory=dir('*.tif');

        imagename=fluo_directory(t).name;
        im=imread(imagename);
        nim=imbinarize(im);

        if mean(mean(nim))<1 %frames without 647 will have a mean of 1
            lidx=t;
            break
        end
    end

    %pre-allocate variables
    intensity=[];
    lCell=[];
    positions=[0];
    
    for i=1:length(datadir)
        cd(datadir(i).folder)
        load(datadir(i).name, 'icell_intensity', 'lcell', 'time')
        count=0;
        
        for n=1:height(icell_intensity)
            if ~isnan(icell_intensity(n, lidx-1))&~isnan(icell_intensity(n, end))
                intensity=[intensity; icell_intensity(n, :)];
                lCell=[lCell; lcell(n, :)];
                count=count+1;
            end
        end
        
        positions(i+1)=count+positions(i); %this is to distinguish which rows are from which colony
     
    end
    
        adjintensity = intensity-intensity(:, end);

        %I don't know why there is a spike in fluor. upon detergent
        %perfusion, but maybe it's messing with the downstream analysis, so
        %I'll just set that point to NaN
        adj_intensity=adjintensity;
        %adj_intensity(:, lidx)=1;
        
        %normalize to pre-lysis frame instead of initial frame. The
        %earliest lysis frame is lidx
        normintensity=adj_intensity(:, lidx-1:end)./adj_intensity(:,lidx-1);
        normintensity(normintensity<0)=0;
        
        %new time vector (0 = just prior to detergent perfusion)
        t=time(lidx-1:end)-time(lidx-1);
        
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
    
%         figure, hold on
%         for i=1:height(intensity)
%             plot(time, intensity(i, :), 'Color', '#77AC30')
%         end
%         ylim([0 Inf])
%         xlabel('Time (minutes)')
%         ylabel('Fluorescence (A.U.)')
%         title(text)
%         if save==1
%             saveas(gcf, [saveAs '_intensityTraces.png'])
%             saveas(gcf, [saveAs '_intensityTraces.fig'])
%         end
        
%         figure, hold on
%         ciplot((intensity_avg-intensity_std),(intensity_avg+intensity_std),time,[0.75 0.75 1])
%         plot(time,intensity_avg,'-r')
%         ylim([0 Inf])
%         xlabel('Time (minutes)')
%         ylabel('Average Fluorescence (A.U.)')
%         title(text)
        
%         if save==1
%             saveas(gcf, [saveAs '_intensityAvg.png'])
%             saveas(gcf, [saveAs '_intensityAvg.fig'])
%         end

        %pause, close all
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
    ylim([0 Inf])
    title(text)
    if save==1
        saveas(gcf, [saveAs '_adjintensityTraces.png'])
        saveas(gcf, [saveAs '_adjintensityTraces.fig'])
    end

    figure(2), hold on
    ciplot((adjintensity_avg-adjintensity_std),(adjintensity_avg+adjintensity_std),time,[0.75 0.75 1])
    plot(time,adjintensity_avg,'-r')
    ylim([0 Inf])
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
    
%     figure, hold on
%    for i=1:height(normintensity)
%         plot(time, normintensity(i, :), 'Color', '#77AC30')
%     end
%     ylim([0 Inf])
%     xlabel('Time (minutes)')
%     ylabel('Normalized Fluorescence (A.U.)')
%     title(text)
%     if save==1
%         saveas(gcf, [saveAs '_normintensityTraces.png'])
%         saveas(gcf, [saveAs '_normintensityTraces.fig'])
%     end
% 
%     figure(2), hold on
%     ciplot((normintensity_avg-normintensity_std),(normintensity_avg+normintensity_std),time,[0.75 0.75 1])
%     plot(time,normintensity_avg,'-r')
%     ylim([0 Inf])
%     xlabel('Time (minutes)')
%     ylabel('Normalized Fluorescence (A.U.)')
%     title(text)
% 
%     if save==1
%         saveas(gcf, [saveAs '_normintensityAvg.png'])
%         saveas(gcf, [saveAs '_normintensityAvg.fig'])
%     end
%         
%     pause, close all
end

function [fit, tau, yhat]=expModel(time, normintensity, text, save, saveAs)
  
    modelfun=@(tau,x)exp(-x./tau);
    tau0=1;

    fit=[];
    tau=[];
    yhat=[];
    for i=1:height(normintensity)
            idx=find(normintensity(i,:)<=1);
            tau_temp=nlinfit(time(idx), normintensity(i,idx), modelfun, tau0);
            y_hat=modelfun(tau_temp, time);   
            
            residuals=abs(normintensity(i,:)-y_hat);
            est=residuals./normintensity(i,:);
            est(est==Inf)=0;
            
            %if max(est)<5
                fit=[fit, i];
                tau=[tau, tau_temp];
                yhat=[yhat; y_hat]; 
            %end
    end

%     figure, hold on
%     for i=1:length(fit)
%         plot(time, normintensity(fit(i),:),'Color', '#77AC30')
%         plot(time, yhat(i,:), '--k')
%     end
%     ylim([0 Inf])
%     legend({'data', 'model'})
%     title(text)
%     xlabel('Time (min)')
%     ylabel('Normalized Intensity (A.U.)')
%     if save==1
%         saveas(gcf, [saveAs '_expFit.png'])
%         saveas(gcf, [saveAs '_expFit.fig'])
%     end
%         
%     pause, close all
end

function [residuals, est]=residualPlot(normintensity, yhat, time, fit, text, save, saveAs)

    residuals=abs(normintensity(fit,:)-yhat);
    est=residuals./normintensity(fit, :);
    est(est==Inf)=0;


%     figure, hold on
%     for i=1:height(est)
%         plot(time, est(i,:), 'Color', '#7E2F8E')
%     end
%     xlabel('Time (minutes)')
%     ylabel('$\frac{|y-\hat{y}|}{y}$','Interpreter','latex', 'FontSize', 20)
%     title(text)
% 
%     
%     if save==1
%         saveas(gcf, [saveAs '_est.png'])
%         saveas(gcf, [saveAs '_est.fig'])
%     end
% 
%     pause, close all
end

function [xq, vq, thalf]=interpPlot(normintensity,time, text, save, saveAs)

    xq=[0:0.0167:time(end)];
    vq=nan(height(normintensity), length(xq));
    
    for i=1:height(normintensity)
        vq(i,:)=interp1(time,normintensity(i,:),xq);
    end
    
%     figure, hold on
%     for i=1:height(normintensity)
%         plot(xq,vq(i,:),':.', 'Color', '#D95319');
%         plot(time,normintensity(i,:),'o', 'Color','#0072BD');
%     end
%     ylim([0 Inf])
%     xlabel('Time (minutes)')
%     ylabel('Normalized Fluoresence (A.U.)')
%     legend({'interpolation', 'data'})
%     title([text ', Linear Interpolation']);

    [~, tidx]=min(abs(vq-0.5), [], 2);
    thalf=xq(tidx);
%     
%     if save==1
%         saveas(gcf, [saveAs '_interp.png'])
%         saveas(gcf, [saveAs '_interp.fig'])
%     end
%     
%     close all
end

function [Cnew, dCP, unb_frac, Cbl_exp, gamma] = photoCorrect(normintensity, t, tau, parameter)
    
    %calculate dt (there's variation of dt during initial values, so it's
    %easier to use end values)
    dt=abs(t(end-1)-t(end));
    beta = 1/((16.0806)*dt + 0.4234); %1/tau
    alpha=1/beta/dt; %tau*dt
    
    %assume that the initial 'measured' fluorescence values and corrected
    %fluor. values will be equal. I prefer to pre-allocate with nan in case 
    %some values are missing in the raw data
    Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
    Cnew(:, 1)=normintensity(:, 1);

    %this is the dCP, or loss attributable to permeability
    dCP=nan(height(normintensity), length(t)-1);

    unb_frac=nan(size(normintensity)); %fraction of unvbleached fluor. 
    unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

    Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
    Cbl_exp(:, 1)=0;
    
    for n=1:height(normintensity)
        
        if parameter==1
            for i=1:length(t)-1

                dCB = normintensity(n,i)/(alpha); %this is the amount of photobleaching that occured in our measured value
    %             if dCB<0 %no negative dCB for normintensity1
    %                 disp(['negative dCB for cell ' num2str(n) ', t = ' num2str(i)])
    %             end

                dCT = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
    %             if dCT>0 %many positive dCT for normintensity1
    %                 disp(['positive dCT for cell ' num2str(n) ', t = ' num2str(i)])
    %                 %dCT=-dCT;
    %                 %continue 
    %             end

                dCP(n, i) = dCT + dCB; %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores

            end
        elseif parameter==2
            
              for i=1:length(t)-1

                dCB = ((normintensity(n,i)+normintensity(n,i+1))/2)/(alpha); %this is the amount of photobleaching that occured in our measured value
    %             if dCB<0 %no negative dCB for normintensity1
    %                 disp(['negative dCB for cell ' num2str(n) ', t = ' num2str(i)])
    %             end

                dCT = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
    %             if dCT>0 %many positive dCT for normintensity1
    %                 disp(['positive dCT for cell ' num2str(n) ', t = ' num2str(i)])
    %                 %dCT=-dCT;
    %                 %continue 
    %             end

                dCP(n, i) = dCT + dCB; %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores

              end
              
        elseif parameter==3
            
              for i=1:length(t)-1

                dCB = normintensity(n,i+1)/(alpha); %this is the amount of photobleaching that occured in our measured value
    %             if dCB<0 %no negative dCB for normintensity1
    %                 disp(['negative dCB for cell ' num2str(n) ', t = ' num2str(i)])
    %             end

                dCT = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
    %             if dCT>0 %many positive dCT for normintensity1
    %                 disp(['positive dCT for cell ' num2str(n) ', t = ' num2str(i)])
    %                 %dCT=-dCT;
    %                 %continue 
    %             end

                dCP(n, i) = dCT + dCB; %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores

              end    
        end
        
    end
    
    gamma = nan(height(normintensity), 1);
    for n=1:height(normintensity)
        p = dCP(n,end)/dt;
        gamma(n,1) = p/Cnew(n, end-1);
    end
    
end

function ciplot(lower,upper,x,colour,Trans)
     
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