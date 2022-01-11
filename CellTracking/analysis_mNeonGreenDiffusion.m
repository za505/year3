%Author: Zarina Akbary
%Date: 12/17/2021
%Purpose: To combine and analyze decayMeasure (dm.mat) data from mNeonGreen diffusion
%experiments. 

clear, close all
%% Generating Figure 1a
%to demonstrate that the cell wall is rate-limiting

% colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#E872D0'};
% colorcode2={'#A1D5F7', '#FFBB9E', '#FFE3A3', '#EAA6F7', '#CBF09C', '#C7EEFF', '#FFB3C1', '#FACAF0'};

colorcode={[0 0.45 0.74], [0.85 0.33 0.1], [0.49 0.18 0.56], [0.93 0.69 0.13], [0.47 0.67 0.19], [0.3 0.75 0.93], [0.64 0.08 0.18], [0.91 0.45 0.82]};
colorcode2={[0.63 0.84 0.97], [1 0.73 0.62], [0.92 0.65 0.97], [1 0.72 0.74], [0.8 0.94 0.61], [0.78 0.93 1], [1 0.7 0.76], [0.98 0.79 0.94]};

%first, direct the code to the location of the dm.mat files 
dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/12112021_analysis';
cd([dirsave '/MatFiles'])

LBa=dir(['10232021_Exp1' '*dm.mat']); %LB, rep 1
LBb=dir(['10262021_Exp1' '*dm.mat']); %LB, rep 2
PBS60a=dir(['10232021_Exp2' '*dm.mat']); %PBS 60 min, rep 1
PBS60b=dir(['10262021_Exp2' '*dm.mat']); %PBS 60 min, rep 2
PBS20a=dir(['11192021_Exp1' '*dm.mat']); %PBS 20 min, rep 1
PBS20b=dir(['11302021_Exp1' '*dm.mat']); %PBS 20 min, rep 2
PBS2a=dir(['11192021_Exp2' '*dm.mat']); %PBS 2 min, rep 1
PBS2b=dir(['12082021_Exp3' '*dm.mat']); %PBS 2 min, rep 2

%% Perform analysis with traces normalized to pre-lysis frame
%calculate normalized fluorescence traces
[normintensity_LBa, intensity_LBa, time_LBa, tme_LBa]=dataNormalize(LBa, 5);
[normintensity_LBb, intensity_LBb, time_LBb, tme_LBb]=dataNormalize(LBb, 4);
[normintensity_PBS60a, intensity_PBS60a, time_PBS60a, tme_PBS60a]=dataNormalize(PBS60a, 2);
[normintensity_PBS60b, intensity_PBS60b, time_PBS60b, tme_PBS60b]=dataNormalize(PBS60b, 4);
[normintensity_PBS20a, intensity_PBS20a, time_PBS20a, tme_PBS20a]=dataNormalize(PBS20a, 3);
[normintensity_PBS20b, intensity_PBS20b, time_PBS20b, tme_PBS20b]=dataNormalize(PBS20b, 3);
[normintensity_PBS2a, intensity_PBS2a, time_PBS2a, tme_PBS2a]=dataNormalize(PBS2a, 9);
[normintensity_PBS2b, intensity_PBS2b, time_PBS2b, tme_PBS2b]=dataNormalize(PBS2b, 9);

%combine datasets
time_LB=time_LBb;
time_PBS60=time_PBS60a;
time_PBS20=time_PBS20a;
time_PBS2=time_PBS2a;

tme_LB=tme_LBb;
tme_PBS60=tme_PBS60a;
tme_PBS20=tme_PBS20a;
tme_PBS2=tme_PBS2a;

intensity_LB=[intensity_LBa(:, 1:length(time_LB)); intensity_LBb];
intensity_PBS60=[intensity_PBS60a; intensity_PBS60b(:, 1:length(time_PBS60))];
intensity_PBS20=[intensity_PBS20a; intensity_PBS20b(:, 1:length(time_PBS20))];
intensity_PBS2=[intensity_PBS2a; intensity_PBS2b(:, 1:length(time_PBS2))];

normintensity_LB=[normintensity_LBa(:, 1:length(tme_LB)); normintensity_LBb];
normintensity_PBS60=[normintensity_PBS60a; normintensity_PBS60b(:, 1:length(tme_PBS60))];
normintensity_PBS20=[normintensity_PBS20a; normintensity_PBS20b(:, 1:length(tme_PBS20))];
normintensity_PBS2=[normintensity_PBS2a; normintensity_PBS2b(:, 1:length(tme_PBS2))];

%% correct for photobleaching
[Cnew_LB, beta_LB, dCB_LB, dCT_LB, dCP_LB, CblExp_LB, unbFrac_LB, midx_LB]=photoCorrect(tme_LB, normintensity_LB, 28.9210);
[Cnew_PBS60, beta_PBS60, dCB_PBS60, dCT_PBS60, dCP_PBS60, CblExp_PBS60, unbFrac_PBS60, midx_PBS60]=photoCorrect(tme_PBS60, normintensity_PBS60, 28.9210);
[Cnew_PBS20, beta_PBS20, dCB_PBS20, dCT_PBS20, dCP_PBS20, CblExp_PBS20, unbFrac_PBS20, midx_PBS20]=photoCorrect(tme_PBS20, normintensity_PBS20, 28.9210);
[Cnew_PBS2, beta_PBS2, dCB_PBS2, dCT_PBS2, dCP_PBS2, CblExp_PBS2, unbFrac_PBS2, midx_PBS2]=photoCorrect(tme_PBS2, normintensity_PBS2, 28.9210);

%% check to see that the photobleach correction has no anomalies 
%do dCT, dCB, and dCP all approach zero? Are the values for dCT and dCP
%negative? Are the values for dCB positive?
% correctionCheck(tme_LB, dCB_LB, dCT_LB, dCP_LB, CblExp_LB, unbFrac_LB)
% correctionCheck(tme_PBS60, dCB_PBS60, dCT_PBS60, dCP_PBS60, CblExp_PBS60, unbFrac_PBS60)  %the Cbl peaks when unbFrac=0.5
% correctionCheck(tme_PBS20, dCB_PBS20, dCT_PBS20, dCP_PBS20, CblExp_PBS20, unbFrac_PBS20)  %the Cbl peaks when unbFrac=0.5
% correctionCheck(tme_PBS2, dCB_PBS2, dCT_PBS2, dCP_PBS2, CblExp_PBS2, unbFrac_PBS2)

%% fit corrected plot to exponential decay function to calculate tau
tau0_LB=30;
tau0_PBS60=1;
tau0_PBS20=10;
tau0_PBS2=25;
[tau_LB, yhat_LB, fidx_LB]=expFit(tme_LB, Cnew_LB, beta_LB, tau0_LB);
[tau_PBS60, yhat_PBS60, fidx_PBS60]=expFit(tme_PBS60, Cnew_PBS60, beta_PBS60, tau0_PBS60);
[tau_PBS20, yhat_PBS20, fidx_PBS20]=expFit(tme_PBS20, Cnew_PBS20, beta_PBS20, tau0_PBS20);
[tau_PBS2, yhat_PBS2, fidx_PBS2]=expFit(tme_PBS2, Cnew_PBS2, beta_PBS2, tau0_PBS2);

%figure out how well tau fits
%comparePlot(Cnew_LB, yhat_LB, tme_LB, fidx_LB), pause, close 
%comparePlot(Cnew_PBS60, yhat_PBS60, tme_PBS60, fidx_PBS60), pause, close
%comparePlot(Cnew_PBS20, yhat_PBS20, tme_PBS20, fidx_PBS20), pause, close
%comparePlot(Cnew_PBS2, yhat_PBS2, tme_PBS2, fidx_PBS2), pause, close
%% Plot the combined values
cd(dirsave)
% %plot intensity values
% figure, hold on
% %ciplot(mean(intensity_LB, 1, 'omitnan')-std(intensity_LB, 0, 1, 'omitnan'), mean(intensity_LB, 1, 'omitnan')+std(intensity_LB, 0, 1, 'omitnan'), time_LB, colorcode2{1})
% plot(time_LB, mean(intensity_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% %ciplot(mean(intensity_PBS60, 1, 'omitnan')-std(intensity_PBS60, 0, 1, 'omitnan'), mean(intensity_PBS60, 1, 'omitnan')+std(intensity_PBS60, 0, 1, 'omitnan'), time_PBS60, colorcode2{3})
% plot(time_PBS60, mean(intensity_PBS60, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% 
% %ciplot(mean(intensity_PBS20, 1, 'omitnan')-std(intensity_PBS20, 0, 1, 'omitnan'), mean(intensity_PBS20, 1, 'omitnan')+std(intensity_PBS20, 0, 1, 'omitnan'), time_PBS20, colorcode2{5})
% plot(time_PBS20, mean(intensity_PBS20, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% 
% %ciplot(mean(intensity_PBS2, 1, 'omitnan')-std(intensity_PBS2, 0, 1, 'omitnan'), mean(intensity_PBS2, 1, 'omitnan')+std(intensity_PBS2, 0, 1, 'omitnan'), time_PBS2, colorcode2{7})
% plot(time_PBS2, mean(intensity_PBS2, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)
% %legend('LB, rep 1, std', 'LB, rep 1, mean', 'LB, rep 2, std', 'LB, rep 2, mean', 'PBS 1 hour, rep 1, std', 'PBS 1 hour, rep 1, mean', 'PBS 1 hour, rep 2, std', 'PBS 1 hour, rep 2, mean', 'PBS 20 min, rep 1, std', 'PBS 20 min, rep 1, mean', 'PBS 20 min, rep 2, std', 'PBS 20 min, rep 2, mean', 'PBS 2 min, rep 1, std', 'PBS 2 min, rep 1, mean', 'PBS 2 min, rep 2, std', 'PBS 2 min, rep 2, mean')
% %legend('LB, std', 'LB, mean', 'PBS 1 hour, std', 'PBS 1 hour, mean',  'PBS 20 min, std', 'PBS 20 min, mean', 'PBS 2 min, std', 'PBS 2 min, mean')
% legend('LB, mean', 'PBS 1 hour, mean',  'PBS 20 min, mean', 'PBS 2 min, mean')
% xlabel('Time (minutes)')
% ylabel('Fluorescence (A.U.)')
% 
% %plot normalized values
% figure, hold on
% %ciplot(mean(normintensity_LB, 1, 'omitnan')-std(normintensity_LB, 0, 1, 'omitnan'), mean(normintensity_LB, 1, 'omitnan')+std(normintensity_LB, 0, 1, 'omitnan'), tme_LB, colorcode2{1})
% plot(tme_LB, mean(normintensity_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% %ciplot(mean(normintensity_PBS60, 1, 'omitnan')-std(normintensity_PBS60, 0, 1, 'omitnan'), mean(normintensity_PBS60, 1, 'omitnan')+std(normintensity_PBS60, 0, 1, 'omitnan'), tme_PBS60, colorcode2{3})
% plot(tme_PBS60, mean(normintensity_PBS60, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% 
% %ciplot(mean(normintensity_PBS20, 1, 'omitnan')-std(normintensity_PBS20, 0, 1, 'omitnan'), mean(normintensity_PBS20, 1, 'omitnan')+std(normintensity_PBS20, 0, 1, 'omitnan'), tme_PBS20, colorcode2{5})
% plot(tme_PBS20, mean(normintensity_PBS20, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% 
% %ciplot(mean(normintensity_PBS2, 1, 'omitnan')-std(normintensity_PBS2, 0, 1, 'omitnan'), mean(normintensity_PBS2, 1, 'omitnan')+std(normintensity_PBS2, 0, 1, 'omitnan'), tme_PBS2, colorcode2{7})
% plot(tme_PBS2, mean(normintensity_PBS2, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)
% 
% %legend('LB, rep 1, std', 'LB, rep 1, mean', 'LB, rep 2, std', 'LB, rep 2, mean', 'PBS 1 hour, rep 1, std', 'PBS 1 hour, rep 1, mean', 'PBS 1 hour, rep 2, std', 'PBS 1 hour, rep 2, mean', 'PBS 20 min, rep 1, std', 'PBS 20 min, rep 1, mean', 'PBS 20 min, rep 2, std', 'PBS 20 min, rep 2, mean', 'PBS 2 min, rep 1, std', 'PBS 2 min, rep 1, mean', 'PBS 2 min, rep 2, std', 'PBS 2 min, rep 2, mean')
% %legend('LB, std', 'LB, mean', 'PBS 1 hour, std', 'PBS 1 hour, mean',  'PBS 20 min, std', 'PBS 20 min, mean', 'PBS 2 min, std', 'PBS 2 min, mean')
% legend('LB, mean', 'PBS 1 hour, mean',  'PBS 20 min, mean', 'PBS 2 min, mean')
% xlabel('Time (minutes)')
% ylabel('Normalized Fluorescence (A.U.)')
% ylim([0 1.05])
% 
%plot corrected values
alpha=0.5;
figure, hold on
ciplot(mean(Cnew_LB(:, 1:50), 1, 'omitnan')-std(Cnew_LB(:, 1:50), 0, 1, 'omitnan'), mean(Cnew_LB(:, 1:50), 1, 'omitnan')+std(Cnew_LB(:, 1:50), 0, 1, 'omitnan'), tme_LB(1:50), colorcode2{1}, alpha)
plot(tme_LB(1:50), mean(Cnew_LB(:, 1:50), 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(Cnew_PBS2, 1, 'omitnan')-std(Cnew_PBS2, 0, 1, 'omitnan'), mean(Cnew_PBS2, 1, 'omitnan')+std(Cnew_PBS2, 0, 1, 'omitnan'), tme_PBS2, colorcode2{7}, alpha)
plot(tme_PBS2, mean(Cnew_PBS2, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

ciplot(mean(Cnew_PBS20, 1, 'omitnan')-std(Cnew_PBS20, 0, 1, 'omitnan'), mean(Cnew_PBS20, 1, 'omitnan')+std(Cnew_PBS20, 0, 1, 'omitnan'), tme_PBS20, colorcode2{5}, alpha)
plot(tme_PBS20, mean(Cnew_PBS20, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(Cnew_PBS60, 1, 'omitnan')-std(Cnew_PBS60, 0, 1, 'omitnan'), mean(Cnew_PBS60, 1, 'omitnan')+std(Cnew_PBS60, 0, 1, 'omitnan'), tme_PBS60, colorcode2{3}, alpha)
plot(tme_PBS60, mean(Cnew_PBS60, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

%legend('LB, rep 1, std', 'LB, rep 1, mean', 'LB, rep 2, std', 'LB, rep 2, mean', 'PBS 1 hour, rep 1, std', 'PBS 1 hour, rep 1, mean', 'PBS 1 hour, rep 2, std', 'PBS 1 hour, rep 2, mean', 'PBS 20 min, rep 1, std', 'PBS 20 min, rep 1, mean', 'PBS 20 min, rep 2, std', 'PBS 20 min, rep 2, mean', 'PBS 2 min, rep 1, std', 'PBS 2 min, rep 1, mean', 'PBS 2 min, rep 2, std', 'PBS 2 min, rep 2, mean')
legend({'LB, std', 'LB, mean', 'PBS 2 min, std', 'PBS 2 min, mean', 'PBS 20 min, std', 'PBS 20 min, mean', 'PBS 1 hour, std', 'PBS 1 hour, mean'}, 'Location', 'southeast')
%legend('LB, mean', 'PBS 1 hour, mean',  'PBS 20 min, mean', 'PBS 2 min, mean')
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
ylim([0 1.05])
% %saveas(gcf, 'correctedIntensity.png')
% %saveas(gcf, 'correctedIntensity.fig')

%plot beta values vs PBS incubation
linearCoef1 = polyfit([zeros(1, length(beta_LB)), repelem(2, length(beta_PBS2)), repelem(20, length(beta_PBS20)), repelem(60, length(beta_PBS60))],[beta_LB', beta_PBS2', beta_PBS20', beta_PBS60'],1);
linearFit1= polyval(linearCoef1,[0 2 20 60]);

linearCoef1 = polyfit([0 2 20 60], [mean(beta_LB, 'omitnan'), mean(beta_PBS2, 'omitnan'), mean(beta_PBS20, 'omitnan'), mean(beta_PBS60, 'omitnan')], 1);
linearFit1 = polyval(linearCoef1, [0 2 20 60]);

%plot beta vs PBS incubation
figure, hold on
scatter(zeros(1, length(beta_LB)), beta_LB', 'MarkerFaceColor', colorcode2{1}, 'MarkerEdgeColor', colorcode2{1})
scatter(repelem(60, length(beta_PBS60)), beta_PBS60', 'MarkerFaceColor', colorcode2{1}, 'MarkerEdgeColor', colorcode2{1})
scatter(repelem(20, length(beta_PBS20)), beta_PBS20', 'MarkerFaceColor', colorcode2{1}, 'MarkerEdgeColor', colorcode2{1})
scatter(repelem(2, length(beta_PBS2)), beta_PBS2', 'MarkerFaceColor', colorcode2{1}, 'MarkerEdgeColor', colorcode2{1})
scatter([0 2 20 60], [mean(beta_LB, 'omitnan'), mean(beta_PBS2, 'omitnan'), mean(beta_PBS20, 'omitnan'), mean(beta_PBS60, 'omitnan')], 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')

errorbar(0, mean(beta_LB, 'omitnan'), std(beta_LB, 'omitnan'), 'Color', 'black')
errorbar(2, mean(beta_PBS2, 'omitnan'), std(beta_PBS2, 'omitnan'), 'Color', 'black')
errorbar(20, mean(beta_PBS20, 'omitnan'), std(beta_PBS20, 'omitnan'), 'Color', 'black')
errorbar(60, mean(beta_PBS60, 'omitnan'), std(beta_PBS60, 'omitnan'), 'Color', 'black')

plot([0 2 20 60], linearFit1, '--b')

xlim([-2 62])
ylim([0 1.1])
xlabel('Time in PBS (minutes)')
ylabel('Beta (A.U.)')

% %saveas(gcf, 'beta_vs_PBS.png')
% %saveas(gcf, 'beta_vs_PBS.fig')

% %plot tau as a function of PBS incubation
linearCoef2 = polyfit([0 2 20 60], [mean(tau_LB, 'omitnan'), mean(tau_PBS2, 'omitnan'), mean(tau_PBS20, 'omitnan'), mean(tau_PBS60, 'omitnan')], 1);
linearFit2 = polyval(linearCoef2, [0 2 20 60]);

figure, hold on
scatter(zeros(1, length(tau_LB)), tau_LB, 'MarkerFaceColor', colorcode2{1}, 'MarkerEdgeColor', colorcode2{1})
scatter(repelem(60, length(tau_PBS60)), tau_PBS60, 'MarkerFaceColor', colorcode2{1}, 'MarkerEdgeColor', colorcode2{1})
scatter(repelem(20, length(tau_PBS20)), tau_PBS20, 'MarkerFaceColor', colorcode2{1}, 'MarkerEdgeColor', colorcode2{1})
scatter(repelem(2, length(tau_PBS2)), tau_PBS2, 'MarkerFaceColor', colorcode2{1}, 'MarkerEdgeColor', colorcode2{1})
scatter([0 2 20 60], [mean(tau_LB, 'omitnan'), mean(tau_PBS2, 'omitnan'), mean(tau_PBS20, 'omitnan'), mean(tau_PBS60, 'omitnan')], 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')

errorbar(0, mean(tau_LB, 'omitnan'), std(tau_LB, 'omitnan'), 'Color', 'black')
errorbar(2, mean(tau_PBS2, 'omitnan'), std(tau_PBS2, 'omitnan'), 'Color', 'black')
errorbar(20, mean(tau_PBS20, 'omitnan'), std(tau_PBS20, 'omitnan'), 'Color', 'black')
errorbar(60, mean(tau_PBS60, 'omitnan'), std(tau_PBS60, 'omitnan'), 'Color', 'black')

plot([0 2 20 60], linearFit2, '--b')

xlim([-2 62])
%ylim([0 1.1])
xlabel('Time in PBS (minutes)')
ylabel('Tau (minutes)')

% %saveas(gcf, 'tau_vs_PBS.png')
% %saveas(gcf, 'tau_vs_PBS.fig')
%% plot correction variables
% figure, plot(tme_LBa(1:end-1), dCB_LBa)
% figure, plot(tme_LBa(1:end-1), dCT_LBa)
% figure, plot(tme_LBa(1:end-1), dCP_LBa)
% figure, plot(tme_LBa, Cbl_exp_LBa)
% figure, plot(tme_LBa, unb_frac_LBa)
% 
% figure, plot(tme_PBS20a(1:end-1), dCB_PBS20a)
% figure, plot(tme_PBS20a(1:end-1), dCT_PBS20a)
% figure, plot(tme_PBS20a(1:end-1), dCP_PBS20a)
% figure, plot(tme_PBS20a, Cbl_exp_PBS20a)
% figure, plot(tme_PBS20a, unb_frac_PBS20a)

% %plot normalized traces vs corrected traces
% cd([dirsave '/postLysis'])
% comparePlot(normintensity_LBa, Cnew_LBa, tme_LBa)
% %saveas(gcf, '10232021_Exp1_rawCorrectedTraces.fig'), %saveas(gcf, '10232021_Exp1_rawCorrectedTraces.png'), close
% comparePlot(normintensity_LBb, Cnew_LBb, tme_LBb)
% %saveas(gcf, '10262021_Exp1_rawCorrectedTraces.fig'), %saveas(gcf, '10262021_Exp1_rawCorrectedTraces.png'), close
% comparePlot(normintensity_PBS60a, Cnew_PBS60a, tme_PBS60a)
% %saveas(gcf, '10232021_Exp2_rawCorrectedTraces.fig'), %saveas(gcf, '10232021_Exp2_rawCorrectedTraces.png'), close
% comparePlot(normintensity_PBS60b, Cnew_PBS60b, tme_PBS60b)
% %saveas(gcf, '10262021_Exp2_rawCorrectedTraces.fig'), %saveas(gcf, '10262021_Exp2_rawCorrectedTraces.png'), close
% comparePlot(normintensity_PBS20a, Cnew_PBS20a, tme_PBS20a)
% %saveas(gcf, '11192021_Exp1_rawCorrectedTraces.fig'), %saveas(gcf, '11192021_Exp1_rawCorrectedTraces.png'), close
% comparePlot(normintensity_PBS20b, Cnew_PBS20b, tme_PBS20b)
% %saveas(gcf, '11302021_Exp1_rawCorrectedTraces.fig'), %saveas(gcf, '11302021_Exp1_rawCorrectedTraces.png'), close
% comparePlot(normintensity_PBS2a, Cnew_PBS2a, tme_PBS2a)
% %saveas(gcf, '11192021_Exp2_rawCorrectedTraces.fig'), %saveas(gcf, '11192021_Exp2_rawCorrectedTraces.png'), close
% comparePlot(normintensity_PBS2b, Cnew_PBS2b, tme_PBS2b)
% %saveas(gcf, '12082021_Exp3_rawCorrectedTraces.fig'), %saveas(gcf, '12082021_Exp3_rawCorrectedTraces.png'), close

%% Functions
function [normintensity, intensity, time, tme]=dataNormalize(datadir, imstart)
        
        %pre-allocate variables
        intensity=[];
      
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'time')

            for n=1:height(icell_intensity)

                %make sure there are fluor readings during the initial frame and the final frame, otherwise the adjust. and norm. will look off
                if ~isnan(icell_intensity(n, imstart))&~isnan(icell_intensity(n, end)) 
                    intensity=[intensity; icell_intensity(n, :)];
                end
            end
            
            if i==1
                tme=time(imstart:end)-time(imstart); %new time vector
            end

        end
        
        %interpolate the fluor values during detergent perfusion
        omit=[2:4]+imstart;
        idx=setdiff(1:length(time), omit);
        v=intensity(:, idx);
        vq=nan(height(v), length(omit));
        for n=1:height(v)
            vq(n, :)=interp1(idx, v(n, :), omit);
        end
        intensity(:, omit)=vq;

        %adjust the background
        adjintensity = intensity-intensity(:, end);

        %normalize to the initial pre-lysis frame
        normintensity=adjintensity(:, imstart:end)./adjintensity(:,imstart);
        normintensity(normintensity<0)=0;
        
end

function [Cnew, beta, dCB, dCT, dCP, Cbl_exp, unb_frac, midx]=photoCorrect(tme, normintensity, parameter)
        
        %Correct for photobleaching
        %calculate dt (it's easier to use end values)
        dt=tme(end)-tme(end-1);
        
        %this value is the slope and intercept calculated from the 1.2, 2,
        %and 3 second controls
        alpha=parameter + 0.2482; %multipy dt?
        
        %assume that the initial 'measured' fluorescence values and corrected
        %fluor. values will be equal. I prefer to pre-allocate with nan in case 
        %some values are missing in the raw data
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:, 1);
        
        %pre-allocate variables
        dCB=nan(height(normintensity), length(tme)-1);
        dCT=nan(height(normintensity), length(tme)-1);
        dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability

        unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
        unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

        Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
        Cbl_exp(:, 1)=0;
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                dCB(n,i) = normintensity(n,i)/alpha; %this is the amount of photobleaching that occured in our measured value

                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)*unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
        end
        
        %calculate beta from corrected values
        beta=Cnew(:, end);
        
        %remove outlier values (there cannot be more fluor than what the
        %cell started with)
        midx=find(beta<=0.9);
        beta=beta(midx);
        Cnew=Cnew(midx, :);
        beta_mean=mean(beta, 'omitnan');
        
end

function correctionCheck(tme, dCB, dCT, dCP, Cbl_exp, unb_frac)

figure('Name', 'dCT vs time')
plot(tme(1:end-1), dCT)
xlabel('Time (minutes)')
ylabel('dCT')
title('Does dCT go to zero? Is it always negative?')
pause, close

figure('Name', 'dCB vs time')
plot(tme(1:end-1), dCB)
xlabel('Time (minutes)')
ylabel('dCB')
title('Does dCB go to zero? Is it always positive?')
pause, close

figure('Name', 'dCP vs time')
plot(tme(1:end-1), dCP)
xlabel('Time (minutes)')
ylabel('dCP')
title('Does dCP go to zero? Is it always negative?')
pause, close

figure('Name', 'unb vs Clb vs time'), hold on
plot(tme, unb_frac, '-b')
plot(tme, Cbl_exp, '-k')
xlabel('Time (minutes)')
title('Does the unbleached fraction reach zero at the same time the bleached concentration peaks?')
pause, close
end

function [tau, yhat, fidx]=expFit(tme, Cnew, beta0, tau0)
        %fit normalized traces to exponential decay function
        %y=(1-beta)*exp(-t./tau)+beta
        %modelfun=@(tau,t)(1-beta)*exp(-t./tau)+beta;
        
        %pre-allocate variables
        tau=[];
        yhat=[];
        fidx=[];
        
        for i=1:height(Cnew)
     
            if ~isnan(beta0(i,1))
                fidx=[fidx i];
                %modelfun=@(tau,t)(1-parameter)*exp(-t./tau)+parameter;
                parameter=beta0(i,1);
                modelfun=@(tau,t)(1-parameter)*exp(-t./tau)+parameter;
                tau_temp=nlinfit(tme, Cnew(i,:), modelfun, tau0);
                y_hat=modelfun(tau_temp, tme);   

                tau=[tau, tau_temp];
                yhat=[yhat; y_hat]; 
            end
        end       
end

function comparePlot(Cnew, yhat, tme, fidx)

    figure, hold on
    for i=1:length(fidx)
        n=fidx(i);
        plot(tme, Cnew(n, :), '-g')
        plot(tme, yhat(i, :), '-b')
    end
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')

end

%% Graveyard

% %Author: Zarina Akbary
% %Date: 05/04/2021
% %Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
% %experiments. 
% 
% clear, close all
% %% matching fluorescent traces to images
% %there are two questions I am trying to answer 1) do the cells with
% %increasing fluorescent traces look different from cells with decreasing
% %fluorescence? 2) what is the alpha that yields a straight line for the 1.2
% %s frame rate LB control?
% 
% dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/12032021_analysis';
% cd([dirsave '/MatFiles'])
% 
% datadir1=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
% datadir2=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
% datadir3=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s
% 
% %% calculate corrected fluorescence traces
% %testing whether the correction is better when traces are normalized to
% % %initial post-lysis frame
% [normintensity1, adjintensity1, Cnew1, tme1, tau1, yhat1]=photoCorrect(datadir1, 9, 28.9210);
% [normintensity2, adjintensity2, Cnew2, tme2, tau2, yhat2]=photoCorrect(datadir2, 9, 28.9210);
% [normintensity3, adjintensity3, Cnew3, tme3, tau3, yhat3]=photoCorrect(datadir3, 9, 28.9210);
% 
% % %% plot normalized traces vs corrected traces
% comparePlot(normintensity1, Cnew1, tme1)
% comparePlot(normintensity2, Cnew2, tme2)
% comparePlot(normintensity3, Cnew3, tme3)
% 
% comparePlot(normintensity1, yhat1, tme1)
% comparePlot(normintensity2, yhat2, tme2)
% comparePlot(normintensity3, yhat3, tme3)
% 
% % comparePlot(adjintensity1(1:9, 8:358), adjintensity2(:, 8:358), tme2)
% %% calculate alpha based on tau values
% % find the mean and standard deviation of the time constants
% tau_means = [mean(tau1, 'omitnan'), mean(tau2, 'omitnan'), mean(tau3, 'omitnan')];
% tau_std = [std(tau1, 0, 'omitnan'), std(tau2, 0, 'omitnan'), std(tau3, 0, 'omitnan')];
% 
% linearCoef1 = polyfit([repelem(1.2/60, length(tau1)), repelem(2/60, length(tau2)), repelem(3/60, length(tau3))],[tau1, tau2, tau3],1);
% linearFit1= polyval(linearCoef1,[0, repelem(1.2/60, length(tau1)), repelem(2/60, length(tau2)), repelem(3/60, length(tau3))]);
% 
% cd(dirsave)
% % figure, hold on
% % scatter(repelem(1.2/60, length(tau1)), tau1, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
% % scatter(repelem(2/60, length(tau2)), tau2, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
% % scatter(repelem(3/60, length(tau3)), tau3, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
% % 
% % scatter([1.2/60,2/60, 3/60], tau_means, 'MarkerFaceColor', 'black')
% % errorbar(1.2/60, tau_means(1), tau_std(1), 'Color', 'black')
% % errorbar(2/60, tau_means(2), tau_std(2), 'Color', 'black')
% % errorbar(3/60, tau_means(3), tau_std(3), 'Color', 'black')
% % 
% % plot([0, repelem(1.2/60, length(tau1)), repelem(2/60, length(tau2)), repelem(3/60, length(tau3))], linearFit1, '--b')
% % xlim([0, 0.06])
% % ylabel('\tau (min^{-1})')
% % xticks([0 1.2/60 2/60 3/60])
% % xticklabels({'0', '1.2 s', '2 s', '3 s'})
% % title('Tau vs Frame Rate')
% % %saveas(gcf, 'alpha.fig')
% % %saveas(gcf, 'alpha.png')
% 
% beta1=1./tau1;
% beta2=1./tau2;
% beta3=1./tau3;
% 
% beta_means = [mean(beta1, 'omitnan'), mean(beta2, 'omitnan'), mean(beta3, 'omitnan')];
% linearCoef2 = polyfit([repelem(1.2/60, length(beta1)), repelem(2/60, length(beta2)), repelem(3/60, length(beta3))],[beta1, beta2, beta3],1);
% linearFit2= polyval(linearCoef2,[0, repelem(1.2/60, length(beta1)), repelem(2/60, length(beta2)), repelem(3/60, length(beta3))]);
% 
% figure, hold on
% scatter([1.2/60,2/60, 3/60], beta_means, 'MarkerFaceColor', 'black')
% plot([0, repelem(1.2/60, length(tau1)), repelem(2/60, length(tau2)), repelem(3/60, length(tau3))], linearFit2, '--b')
% %% Functions
% function [normintensity, adjintensity, Cnew, tme, tau, yhat]=photoCorrect(datadir, imstart, parameter)
%         
%         %pre-allocate variables
%         intensity=[];
%       
%         %go through the data for each position
%         for i=1:length(datadir)
%             
%             %load decayMeasure .mat file
%             cd(datadir(i).folder)
%             load(datadir(i).name, 'icell_intensity', 'time')
% 
%             for n=1:height(icell_intensity)
% 
%                 %make sure there are fluor readings during the initial frame and the final frame, otherwise the adjust. and norm. will look off
%                 if ~isnan(icell_intensity(n, imstart))&~isnan(icell_intensity(n, end)) 
%                     intensity=[intensity; icell_intensity(n, :)];
%                 end
%             end
% 
%             if i==1
%                 tme=time(imstart:end)-time(imstart); %new time vector
%             end
% 
%         end
% 
%         %adjust the background
%         adjintensity = intensity-intensity(:, end);
% 
%         %normalize to the initial post-lysis frame
%         normintensity=adjintensity(:, imstart:end)./adjintensity(:,imstart);
%         normintensity(normintensity<0)=0;
%         
%         %fit normalized traces to exponential decay function
%         modelfun=@(tau,x)exp(-x./tau);
%         tau0=1;
%         
%         %pre-allocate variables
%         fit=[];
%         tau=[];
%         yhat=[];
% 
%         for i=1:height(normintensity)
% %             %this will skip the detergent perfusion reading
% %             idx=find(normintensity(i,:)<=1); 
%             tau_temp=nlinfit(tme, normintensity(i,:), modelfun, tau0);
%             y_hat=modelfun(tau_temp, tme);   
% 
% %             fit=[fit, i];
%             tau=[tau, tau_temp];
%             yhat=[yhat; y_hat]; 
% 
%         end
%     
%         %calculate dt (there's variation of dt during initial values, so it's
%         %easier to use end values)
%         dt=tme(end)-tme(end-1);
%         %beta = 1/((16.9441)*dt + 0.1823); %1/tau
%         beta = 1/((parameter)*dt + 0.1823); %1/tau
%         alpha=1/beta/dt; %tau*dt
%         
%         %assume that the initial 'measured' fluorescence values and corrected
%         %fluor. values will be equal. I prefer to pre-allocate with nan in case 
%         %some values are missing in the raw data
%         Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
%         Cnew(:, 1)=normintensity(:, 1);
% 
%         %pre-allocate variables
%         dCB=nan(height(normintensity), length(tme)-1);
%         dCT=nan(height(normintensity), length(tme)-1);
%         dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability
% 
%         unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
%         unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point
% 
%         Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
%         Cbl_exp(:, 1)=0;
% 
%         mx=[];
%         for n=1:height(normintensity)
%            
%             for i=1:length(tme)-1
%                 
%                 if normintensity(n,i)<0.05
%                     continue
%                 end
%                 
%                 dCB(n,i) = normintensity(n,i)/(alpha); %this is the amount of photobleaching that occured in our measured value
% 
%                 dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
% 
%                 dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability
% 
%                 dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached
% 
%                 Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);
% 
%                 Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
%                 
%                 unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
%                 
%             end  
%             
%             %when does the loop break (Cnew no longer calculated)
%             mx=[mx,i];
%         end
%         
%         %mean time stamp for Cnew
%         idx=round(mean(mx));
%         
%         gamma = nan(height(normintensity), 1);
%     
%         for n=1:height(normintensity)
%             p = dCP(n,end)/dt;
%             gamma(n,1) = p/Cnew(n, end-1);
%         end
%         
% end
% 
% function comparePlot(normintensity, Cnew, tme)
% 
%     figure, hold on
%     for n=1:height(normintensity)
%         plot(tme, normintensity(n, :), '-b')
%         plot(tme, Cnew(n, :), '-g')
%     end
%     xlabel('Time (minutes)')
%     ylabel('Normalized Fluorescence (A.U.)')
% 
% end
% 
% function compareImage(datadir, imstart, idx, B, positions, Cnew)
%       
%     mdx=length(Cnew)/2;
%     
%     for i=1:length(datadir)
%         
%         if i==1
%             s1=1;
%             s2=positions(i);
%         else
%             s1=positions(i-1);
%             s2=positions(i);
%         end
%         
%         %first, load the decay measure data
%         cd(datadir(i).folder)
%         load(datadir(i).name, 'fluo_directory', 'channels')
%         
%         cd(channels{1});
%         imagename=fluo_directory{1}(imstart + idx).name;
%         im=imread(imagename);
%         
%         figure, imshow(im, []), hold on
%         for n=s1:s2
%             if Cnew(n, mdx)>1
%                 plot(B{n,mdx}(:,1),B{n,mdx}(:,2),'-r')
%             else
%                 plot(B{n,mdx}(:,1),B{n,mdx}(:,2),'-b')
%             end
%         end
%         
%         pause, close all
%     end
% end
