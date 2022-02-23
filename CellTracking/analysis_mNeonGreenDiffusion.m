%Author: Zarina Akbary
%Date: 02/12/2022
%Purpose: To combine and analyze decayMeasure (dm.mat) data from mNeonGreen diffusion
%experiments. 

clear, close all

%Inputs
%dirsave = directory where all the dm.mat files are located and where plots
%are saved

%Outputs
%intensity = cell x time matrix of raw intensity values; concatenated
%icell_intensity variable from dm.mat files

%time = 1 x time matrix of time (minutes) values that correspond to
%intensity traces

%bgintensity = m x time matrix of background intensity values; concatenated
%bg_intensity variable from dm.mat files

%tme = 1 x n matrix of time (minutes) values that correspond to
%the pre- or post-lysis time point (imstart) and imend

%adjintensity = cell x tme maxtrix of intensity values minus bg intensity
%values from imstart to imend

%normintensity = cell x tme matrix of adjintensity values normalized to the
%initial value

%% User Input 

%color codes must be in RGB [0-1] format to be used in ciplot
colorcode={[204 0 0], [204 102 0], [204 204 0], [102 204 0], [0 204 204], [0 0 204], [102 0 204], [204 0 204], [204 0 102], [255 102 102], [255 178 102], [102 255 102], [102 255 255], [102 178 255], [178 102 255], [255 102 255],[255 102 178]};
colorcode2={[255 51 51], [255 153 51], [255 255 51], [153 255 51], [51 255 255], [51 51 255], [153 51 255], [255 51 255], [255 51 153], [255 204 204], [255 229 204], [204 255 204], [204 255 255], [204 229 255], [229 204 255], [255 204 255],[255 204 229]};

colorcode=cellfun(@(x)(x./255), colorcode, 'UniformOutput', false);
colorcode2=cellfun(@(x)(x./255), colorcode2, 'UniformOutput', false);

transparency = 0.3; %this is the alpha argument for the ciplot function

%location of the dm.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/02102022_reanalysis';
cd([dirsave '/MatFiles'])

LB1a=dir(['10232021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 1
LB1b=dir(['10262021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 2
LB1c=dir(['02212022_Exp2' '*dm.mat']); %LB, fr = 1 min, intensity = 20%
LB5a=dir(['02122022_Exp1' '*dm.mat']); %LB, fr = 5 min
LB10a=dir(['02122022_Exp2' '*dm.mat']); %LB, fr = 10 min
LB20a=dir(['02092022_Exp1' '*dm.mat']); %LB, fr = 20 min
LB20b=dir(['02212022_Exp1' '*dm.mat']); %LB, fr = 20 min, intensity = 20%

PBS2a=dir(['11192021_Exp2' '*dm.mat']); %PBS 2 min, rep 1
PBS2b=dir(['12082021_Exp3' '*dm.mat']); %PBS 2 min, rep 2
PBS20a=dir(['11192021_Exp1' '*dm.mat']); %PBS 20 min, rep 1
PBS20b=dir(['11302021_Exp1' '*dm.mat']); %PBS 20 min, rep 2
PBS60a=dir(['10232021_Exp2' '*dm.mat']); %PBS 60 min, rep 1
PBS60b=dir(['10262021_Exp2' '*dm.mat']); %PBS 60 min, rep 2
PBS120a=dir(['01142022_Exp1' '*dm.mat']); %PBS 120 min, based on WT not BT

LBMga=dir(['01172022_Exp1' '*dm.mat']); %LB, 20 mM Mg2+
LBEa=dir(['01172022_Exp2' '*dm.mat']); %LB, 10 mM EDTA
LBtuna=dir(['01242022_Exp1' '*dm.mat']); %LB, 0.5 ug/mL tunicamycin
LBvana=dir(['01262022_Exp1' '*dm.mat']); %LB, 1 ug/mL vancomycin
LBspnta=dir(['02122022_Exp3' '*dm.mat']); %LB, spent media (10 min frame rate)
LBspntb=dir(['02192022_Exp1' '*dm.mat']); %LB, spent media (10 min frame rate)

LB1s=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
LB2s=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
LB3s=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s

LB1t=dir(['02192022_Exp2' '*dm.mat']); %LB, frame rate = 1.76 s
LB2t=dir(['02192022_Exp3' '*dm.mat']); %LB, frame rate = 2.3 s
LB3t=dir(['02192022_Exp4' '*dm.mat']); %LB, frame rate = 3 s

% i100_1s=dir(['02172022_Exp1_i100_fr001_dm.mat']); %LB, frame rate = 1.2 s
% i100_2s=dir(['02172022_Exp1_i100_fr002_dm.mat']); %LB, frame rate = 2 s
% i100_3s=dir(['02172022_Exp1_i100_fr003_dm.mat']); %LB, frame rate = 3 s
% 
% i075_1s=dir(['02172022_Exp1_i075_fr001_dm.mat']); %LB, frame rate = 1.2 s
% i075_2s=dir(['02172022_Exp1_i075_fr002_dm.mat']); %LB, frame rate = 2 s
% i075_3s=dir(['02172022_Exp1_i075_fr003_dm.mat']); %LB, frame rate = 3 s
% 
% i050_1s=dir(['02172022_Exp1_i050_fr001_dm.mat']); %LB, frame rate = 1.2 s
% i050_2s=dir(['02172022_Exp1_i050_fr002_dm.mat']); %LB, frame rate = 2 s
% i050_3s=dir(['02172022_Exp1_i050_fr003_dm.mat']); %LB, frame rate = 3 s

%% calculate normalized fluorescence traces
[intensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lcell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize(LB1a, 6); %02/15/2022, correct imstart
[intensity_LB1b, adjintensity_LB1b, normintensity_LB1b, lcell_LB1b, time_LB1b, tme_LB1b, imstart_LB1b]=dataNormalize(LB1b, 5); %02/15/2022, correct imstart
[intensity_LB5a, adjintensity_LB5a, normintensity_LB5a, lcell_LB5a, time_LB5a, tme_LB5a, imstart_LB5a]=dataNormalize(LB5a, 6); %02/15/2022, correct imstart
[intensity_LB10a, adjintensity_LB10a, normintensity_LB10a, lcell_LB10a, time_LB10a, tme_LB10a, imstart_LB10a]=dataNormalize(LB10a, 6); %02/15/2022, correct imstart
[intensity_LB20a, adjintensity_LB20a, normintensity_LB20a, lcell_LB20a, time_LB20a, tme_LB20a, imstart_LB20a]=dataNormalize(LB20a, 6); %02/15/2022, correct imstart

[intensity_PBS2a, adjintensity_PBS2a, normintensity_PBS2a, lcell_PBS2a, time_PBS2a, tme_PBS2a, imstart_PBS2a]=dataNormalize(PBS2a, 10); %02/15/2022, correct imstart
[intensity_PBS2b, adjintensity_PBS2b, normintensity_PBS2b, lcell_PBS2b, time_PBS2b, tme_PBS2b, imstart_PBS2b]=dataNormalize(PBS2b, 10); %02/15/2022, correct imstart
[intensity_PBS20a, adjintensity_PBS20a, normintensity_PBS20a, lcell_PBS20a, time_PBS20a, tme_PBS20a, imstart_PBS20a]=dataNormalize(PBS20a, 4); %02/15/2022, correct imstart
[intensity_PBS20b, adjintensity_PBS20b, normintensity_PBS20b, lcell_PBS20b, time_PBS20b, tme_PBS20b, imstart_PBS20b]=dataNormalize(PBS20b, 4); %02/15/2022, correct imstart
[intensity_PBS60a, adjintensity_PBS60a, normintensity_PBS60a, lcell_PBS60a, time_PBS60a, tme_PBS60a, imstart_PBS60a]=dataNormalize(PBS60a, 3); %02/15/2022, correct imstart
[intensity_PBS60b, adjintensity_PBS60b, normintensity_PBS60b, lcell_PBS60b, time_PBS60b, tme_PBS60b, imstart_PBS60b]=dataNormalize(PBS60b, 5); %02/15/2022, correct imstart
[intensity_PBS120a, adjintensity_PBS120a, normintensity_PBS120a, lcell_PBS120a, time_PBS120a, tme_PBS120a, imstart_PBS120a]=dataNormalize(PBS120a, 4); %02/15/2022, correct imstart

[intensity_LBMga, adjintensity_LBMga, normintensity_LBMga, lcell_LBMga, time_LBMga, tme_LBMga, imstart_LBMga]=dataNormalize(LBMga, 6);  %02/15/2022, correct imstart
[intensity_LBEa, adjintensity_LBEa, normintensity_LBEa, lcell_LBEa, time_LBEa, tme_LBEa, imstart_LBEa]=dataNormalize(LBEa, 6); 
[intensity_LBtuna, adjintensity_LBtuna, normintensity_LBtuna, lcell_LBtuna, time_LBtuna, tme_LBtuna, imstart_LBtuna]=dataNormalize(LBtuna, 3); 
[intensity_LBvana, adjintensity_LBvana, normintensity_LBvana, lcell_LBvana, time_LBvana, tme_LBvana, imstart_LBvana]=dataNormalize(LBvana, 5); %the true pre-lysis frame is t=4, but there are no fluor readings 
[intensity_LBspnta, adjintensity_LBspnta, normintensity_LBspnta, lcell_LBspnta, time_LBspnta, tme_LBspnta, imstart_LBspnta]=dataNormalize(LBspnta, 6); %02/15/2022, correct imstart
[intensity_LBspntb, adjintensity_LBspntb, normintensity_LBspntb, lcell_LBspntb, time_LBspntb, tme_LBspntb, imstart_LBspntb]=dataNormalize(LBspntb, 6);

[intensity_LB1s, bgintensity_LB1s, adjintensity_LB1s, normintensity_LB1s, lcell_LB1s, time_LB1s, tme_LB1s, imstart_LB1s]=controlNormalize(LB1s); 
[intensity_LB2s, bgintensity_LB2s, adjintensity_LB2s, normintensity_LB2s, lcell_LB2s, time_LB2s, tme_LB2s, imstart_LB2s]=controlNormalize(LB2s); 
[intensity_LB3s, bgintensity_LB3s, adjintensity_LB3s, normintensity_LB3s, lcell_LB3s, time_LB3s, tme_LB3s, imstart_LB3s]=controlNormalize(LB3s); 

[intensity_LB1t, bgintensity_LB1t, adjintensity_LB1t, normintensity_LB1t, lcell_LB1t, time_LB1t, tme_LB1t, imstart_LB1t]=controlNormalize(LB1t); 
[intensity_LB2t, bgintensity_LB2t, adjintensity_LB2t, normintensity_LB2t, lcell_LB2t, time_LB2t, tme_LB2t, imstart_LB2t]=controlNormalize(LB2t); 
[intensity_LB3t, bgintensity_LB3t, adjintensity_LB3t, normintensity_LB3t, lcell_LB3t, time_LB3t, tme_LB3t, imstart_LB3t]=controlNormalize(LB3t);

% [intensity_i1001s, bgintensity_i1001s, adjintensity_i1001s, normintensity_i1001s, lcell_i1001s, time_i1001s, tme_i1001s, imstart_i1001s]=controlNormalize(i100_1s); 
% [intensity_i1002s, bgintensity_i1002s, adjintensity_i1002s, normintensity_i1002s, lcell_i1002s, time_i1002s, tme_i1002s, imstart_i1002s]=controlNormalize(i100_2s); 
% [intensity_i1003s, bgintensity_i1003s, adjintensity_i1003s, normintensity_i1003s, lcell_i1003s, time_i1003s, tme_i1003s, imstart_i1003s]=controlNormalize(i100_3s);
% 
% [intensity_i0751s, bgintensity_i0751s, adjintensity_i0751s, normintensity_i0751s, lcell_i0751s, time_i0751s, tme_i0751s, imstart_i0751s]=controlNormalize(i075_1s); 
% [intensity_i0752s, bgintensity_i0752s, adjintensity_i0752s, normintensity_i0752s, lcell_i0752s, time_i0752s, tme_i0752s, imstart_i0752s]=controlNormalize(i075_2s); 
% [intensity_i0753s, bgintensity_i0753s, adjintensity_i0753s, normintensity_i0753s, lcell_i0753s, time_i0753s, tme_i0753s, imstart_i0753s]=controlNormalize(i075_3s); 
% 
% [intensity_i0501s, bgintensity_i0501s, adjintensity_i0501s, normintensity_i0501s, lcell_i0501s, time_i0501s, tme_i0501s, imstart_i0501s]=controlNormalize(i050_1s); 
% [intensity_i0502s, bgintensity_i0502s, adjintensity_i0502s, normintensity_i0502s, lcell_i0502s, time_i0502s, tme_i0502s, imstart_i0502s]=controlNormalize(i050_2s); 
% [intensity_i0503s, bgintensity_i0503s, adjintensity_i0503s, normintensity_i0503s, lcell_i0503s, time_i0503s, tme_i0503s, imstart_i0503s]=controlNormalize(i050_3s); 

[intensity_LB1c, adjintensity_LB1c, normintensity_LB1c, lcell_LB1c, time_LB1c, tme_LB1c, imstart_LB1c]=dataNormalize(LB1c, 6); 
[intensity_LB20b, adjintensity_LB20b, normintensity_LB20b, lcell_LB20b, time_LB20b, tme_LB20b, imstart_LB20b]=dataNormalize(LB20b, 6);
%% do the strain and growth rate vary between the untreated and spent cells?
[strain_LB1a]=strainCalc(lcell_LB1a, imstart_LB1a);
[strain_LB10a]=strainCalc(lcell_LB10a, imstart_LB10a);
[strain_LBspnta]=strainCalc(lcell_LBspnta, imstart_LBspnta);
[strain_LBspntb]=strainCalc(lcell_LBspntb, imstart_LBspntb);

[growthRate_LB1a]=gRateCalc(time_LB1a, lcell_LB1a, imstart_LB1a);
[growthRate_LB10a]=gRateCalc(time_LB10a, lcell_LB10a, imstart_LB10a);
[growthRate_LBspnta]=gRateCalc(time_LBspnta, lcell_LBspnta, imstart_LBspnta);
[growthRate_LBspntb]=gRateCalc(time_LBspntb, lcell_LBspntb, imstart_LBspntb);

%% plot to compare
figure, hold on
plot(time_LB1a(1:imstart_LB1a+3), mean(strain_LB1a, 1, 'omitnan'), 'Color', colorcode{1})
plot(time_LB10a(1:imstart_LB10a+3), mean(strain_LB10a, 1, 'omitnan'), 'Color', colorcode{4})
plot(time_LBspnta(1:imstart_LBspnta+3), mean(strain_LBspnta, 1, 'omitnan'), 'Color', colorcode{5})
plot(time_LBspntb(1:imstart_LBspntb+3), mean(strain_LBspntb, 1, 'omitnan'), 'Color', colorcode{8})

figure, hold on
plot(time_LB1a(1:imstart_LB1a+2), mean(growthRate_LB1a, 1, 'omitnan'), 'Color', colorcode{1})
plot(time_LB10a(1:imstart_LB10a+2), mean(growthRate_LB10a, 1, 'omitnan'), 'Color', colorcode{4})
plot(time_LBspnta(1:imstart_LBspnta+2), mean(growthRate_LBspnta, 1, 'omitnan'), 'Color', colorcode{5})
plot(time_LBspntb(1:imstart_LBspntb+2), mean(growthRate_LBspntb, 1, 'omitnan'), 'Color', colorcode{8})

%% Controls
% plot the raw fluorescence intensity 
cd([dirsave '/Controls'])
figure, hold on
plot(time_LB1s, intensity_LB1s, '-r');
plot(time_LB1s, bgintensity_LB1s, '--r')
plot(time_LB2s, intensity_LB2s, '-b');
plot(time_LB2s, bgintensity_LB2s, '--b')
plot(time_LB3s, intensity_LB3s, '-g'); 
plot(time_LB3s, bgintensity_LB3s, '--g')
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
%saveas(gcf, 'control_intensity.png');

% plot the adjusted intensity
figure, hold on
plot(time_LB1s, adjintensity_LB1s, '-r');
plot(time_LB2s, adjintensity_LB2s, '-b');
plot(time_LB3s, adjintensity_LB3s, '-g'); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
%saveas(gcf, 'control_adjintensity.png');

% plot the normalized intensity
figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r');
plot(tme_LB2s, normintensity_LB2s, '-b');
plot(tme_LB3s, normintensity_LB3s, '-g'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'control_normintensity.png');

%% plot the normalized intensity for the agarose pad controls
figure, hold on
plot(tme_i1001s, normintensity_i1001s, '-r');
plot(tme_i1002s, normintensity_i1002s, '-b');
plot(tme_i1003s, normintensity_i1003s, '-g'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_i0751s, normintensity_i0751s, '-r');
plot(tme_i0752s, normintensity_i0752s, '-b');
plot(tme_i0753s, normintensity_i0753s, '-g'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_i0501s, normintensity_i0501s, '-r');
plot(tme_i0502s, normintensity_i0502s, '-b');
plot(tme_i0503s, normintensity_i0503s, '-g'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

%% plot the adjusted intensity for the agarose pad controls
figure, hold on
plot(time_i1001s, adjintensity_i1001s, '-r');
plot(time_i1002s, adjintensity_i1002s, '-b');
plot(time_i1003s, adjintensity_i1003s, '-g'); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence')

figure, hold on
plot(time_i0751s, adjintensity_i0751s, '-r');
plot(time_i0752s, adjintensity_i0752s, '-b');
plot(time_i0753s, adjintensity_i0753s, '-g'); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence')

figure, hold on
plot(time_i0501s, adjintensity_i0501s, '-r');
plot(time_i0502s, adjintensity_i0502s, '-b');
plot(time_i0503s, adjintensity_i0503s, '-g'); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence')

%% calculate alpha
[tau1, yhat_LB1s]=tauCalc(tme_LB1s, normintensity_LB1s, 1);
[tau2, yhat_LB2s]=tauCalc(tme_LB2s, normintensity_LB2s, 1);
[tau3, yhat_LB3s]=tauCalc(tme_LB3s, normintensity_LB3s, 1);

% [tau_i1001s, yhat_i1001s]=tauCalc(tme_i1001s, normintensity_i1001s, 0.8);
% [tau_i1002s, yhat_i1002s]=tauCalc(tme_i1002s, normintensity_i1002s, 2);
% [tau_i1003s, yhat_i1003s]=tauCalc(tme_i1003s, normintensity_i1003s, 2); 
% 
% [tau_i0751s, yhat_i0751s]=tauCalc(tme_i0751s, normintensity_i0751s, 2);
% [tau_i0752s, yhat_i0752s]=tauCalc(tme_i0752s, normintensity_i0752s, 2);
% [tau_i0753s, yhat_i0753s]=tauCalc(tme_i0753s, normintensity_i0753s, 2); 
% 
% [tau_i0501s, yhat_i0501s]=tauCalc(tme_i0501s, normintensity_i0501s, 0.8);
% [tau_i0502s, yhat_i0502s]=tauCalc(tme_i0502s, normintensity_i0502s, 0.01);
% [tau_i0503s, yhat_i0503s]=tauCalc(tme_i0503s, normintensity_i0503s, 0.01); 

%% plots for tau
figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r')
plot(tme_LB1s, yhat_LB1s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'tauFit_LB1s.png')

figure, hold on
plot(tme_LB2s, normintensity_LB2s, '-b')
plot(tme_LB2s, yhat_LB2s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'tauFit_LB2s.png')

figure, hold on
plot(tme_LB3s, normintensity_LB3s, '-g')
plot(tme_LB3s, yhat_LB3s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'tauFit_LB3s.png')

tau_means = [mean(tau1, 'omitnan'), mean(tau2, 'omitnan'), mean(tau3, 'omitnan')];
tau_std = [std(tau1, 0, 'omitnan'), std(tau2, 0, 'omitnan'), std(tau3, 0, 'omitnan')];

dt1=tme_LB1s(end)-tme_LB1s(end-1);
dt2=tme_LB2s(end)-tme_LB2s(end-1);
dt3=tme_LB3s(end)-tme_LB3s(end-1);

linearCoef1 = polyfit([repelem(dt1, length(tau1)), repelem(dt2, length(tau2)), repelem(dt3, length(tau3))],[tau1', tau2', tau3'],1);
linearFit1= polyval(linearCoef1,[0 dt1 dt2 dt3]);

figure, hold on
scatter(repelem(dt1, length(tau1)), tau1, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(dt2, length(tau2)), tau2, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(dt3, length(tau3)), tau3, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')

scatter([dt1, dt2, dt3], tau_means, 'MarkerFaceColor', 'black')
errorbar(dt1, tau_means(1), tau_std(1), 'Color', 'black')
errorbar(dt2, tau_means(2), tau_std(2), 'Color', 'black')
errorbar(dt3, tau_means(3), tau_std(3), 'Color', 'black')

plot([0 dt1 dt2 dt3], linearFit1, '--b')
xlim([0, 0.06])
ylabel('\tau (min^{-1})')
xticks([0 dt1 dt2 dt3])
xticklabels({'0', '1.2 s', '2 s', '3 s'})
title('Tau vs Frame Rate')
%saveas(gcf, 'tau_vs_frameRate.png')

%% plot for agarose pad controls (i100)
figure, hold on
plot(tme_i1001s, normintensity_i1001s, '-r')
plot(tme_i1001s, yhat_i1001s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_i1002s, normintensity_i1002s, '-b')
plot(tme_i1002s, yhat_i1002s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_i1003s, normintensity_i1003s, '-g')
plot(tme_i1003s, yhat_i1003s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

%% plot for agarose pad controls (i075)
figure, hold on
plot(tme_i0751s, normintensity_i0751s, '-r')
plot(tme_i0751s, yhat_i0751s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_i0752s, normintensity_i0752s, '-b')
plot(tme_i0752s(1:23), yhat_i0752s(1:23), '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_i0753s, intensity_i0753s, '-g')
plot(tme_i0753s(1:14), yhat_i0753s(1:14), '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

%% plot for agarose pad controls (i050)
figure, hold on
plot(tme_i0501s, normintensity_i0501s, '-r')
plot(tme_i0501s, yhat_i0501s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_i0502s, normintensity_i0502s, '-b')
plot(tme_i0502s(1:23), yhat_i0502s(1:23), '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_i0503s, intensity_i0503s, '-g')
plot(tme_i0503s(1:14), yhat_i0503s(1:14), '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

%% correct for photobleaching
alpha=32.2114;
intercept=0.1614;

[Cnew_LB1s, dCB_LB1s, dCT_LB1s, dCP_LB1s, CblExp_LB1s, unbFrac_LB1s]=photoCorrect(tme_LB1s, normintensity_LB1s, alpha, intercept);
[Cnew_LB2s, dCB_LB2s, dCT_LB2s, dCP_LB2s, CblExp_LB2s, unbFrac_LB2s]=photoCorrect(tme_LB2s, normintensity_LB2s, alpha, intercept);
[Cnew_LB3s, dCB_LB3s, dCT_LB3s, dCP_LB3s, CblExp_LB3s, unbFrac_LB3s]=photoCorrect(tme_LB3s, normintensity_LB3s, alpha, intercept);

[Cnew_LB1a, dCB_LB1a, dCT_LB1a, dCP_LB1a, CblExp_LB1a, unbFrac_LB1a]=photoCorrect(tme_LB1a, normintensity_LB1a, alpha, intercept);
[Cnew_LB1b, dCB_LB1b, dCT_LB1b, dCP_LB1b, CblExp_LB1b, unbFrac_LB1b]=photoCorrect(tme_LB1b, normintensity_LB1b, alpha, intercept);
[Cnew_LB5a, dCB_LB5a, dCT_LB5a, dCP_LB5a, CblExp_LB5a, unbFrac_LB5a]=photoCorrect(tme_LB5a, normintensity_LB5a, alpha, intercept);
[Cnew_LB10a, dCB_LB10a, dCT_LB10a, dCP_LB10a, CblExp_LB10a, unbFrac_LB10a]=photoCorrect(tme_LB10a, normintensity_LB10a, alpha, intercept);
[Cnew_LB20a, dCB_LB20a, dCT_LB20a, dCP_LB20a, CblExp_LB20a, unbFrac_LB20a]=photoCorrect(tme_LB20a, normintensity_LB20a, alpha, intercept);

[Cnew_PBS2a, dCB_PBS2a, dCT_PBS2a, dCP_PBS2a, CblExp_PBS2a, unbFrac_PBS2a]=photoCorrect(tme_PBS2a, normintensity_PBS2a, alpha, intercept);
[Cnew_PBS2b, dCB_PBS2b, dCT_PBS2b, dCP_PBS2b, CblExp_PBS2b, unbFrac_PBS2b]=photoCorrect(tme_PBS2b, normintensity_PBS2b, alpha, intercept);
[Cnew_PBS20a, dCB_PBS20a, dCT_PBS20a, dCP_PBS20a, CblExp_PBS20a, unbFrac_PBS20a]=photoCorrect(tme_PBS20a, normintensity_PBS20a, alpha, intercept);
[Cnew_PBS20b, dCB_PBS20b, dCT_PBS20b, dCP_PBS20b, CblExp_PBS20b, unbFrac_PBS20b]=photoCorrect(tme_PBS20b, normintensity_PBS20b, alpha, intercept);
[Cnew_PBS60a, dCB_PBS60a, dCT_PBS60a, dCP_PBS60a, CblExp_PBS60a, unbFrac_PBS60a]=photoCorrect(tme_PBS60a, normintensity_PBS60a, alpha, intercept);
[Cnew_PBS60b, dCB_PBS60b, dCT_PBS60b, dCP_PBS60b, CblExp_PBS60b, unbFrac_PBS60b]=photoCorrect(tme_PBS60b, normintensity_PBS60b, alpha, intercept);
[Cnew_PBS120a, dCB_PBS120a, dCT_PBS120a, dCP_PBS120a, CblExp_PBS120a, unbFrac_PBS120a]=photoCorrect(tme_PBS120a, normintensity_PBS120a, alpha, intercept);

[Cnew_LBMga, dCB_LBMga, dCT_LBMga, dCP_LBMga, CblExp_LBMga, unbFrac_LBMga]=photoCorrect(tme_LBMga, normintensity_LBMga, alpha, intercept);
[Cnew_LBEa, dCB_LBEa, dCT_LBEa, dCP_LBEa, CblExp_LBEa, unbFrac_LBEa]=photoCorrect(tme_LBEa, normintensity_LBEa, alpha, intercept);
[Cnew_LBtuna, dCB_LBtuna, dCT_LBtuna, dCP_LBtuna, CblExp_LBtuna, unbFrac_LBtuna]=photoCorrect(tme_LBtuna, normintensity_LBtuna, alpha, intercept);
[Cnew_LBvana, dCB_LBvana, dCT_LBvana, dCP_LBvana, CblExp_LBvana, unbFrac_LBvana]=photoCorrect(tme_LBvana, normintensity_LBvana, alpha, intercept);
[Cnew_LBspnta, dCB_LBspnta, dCT_LBspnta, dCP_LBspnta, CblExp_LBspnta, unbFrac_LBspnta]=photoCorrect(tme_LBspnta, normintensity_LBspnta, alpha, intercept);
[Cnew_LBspntb, dCB_LBspntb, dCT_LBspntb, dCP_LBspntb, CblExp_LBspntb, unbFrac_LBspntb]=photoCorrect(tme_LBspntb, normintensity_LBspntb, alpha, intercept);

%% new photobleaching control
alpha=152.6455;
intercept=-1.0976;

[Cnew_LB1t, dCB_LB1t, dCT_LB1t, dCP_LB1t, CblExp_LB1t, unbFrac_LB1t]=photoCorrect(tme_LB1t, normintensity_LB1t, alpha, intercept);
[Cnew_LB2t, dCB_LB2t, dCT_LB2t, dCP_LB2t, CblExp_LB2t, unbFrac_LB2t]=photoCorrect(tme_LB2t, normintensity_LB2t, alpha, intercept);
[Cnew_LB3t, dCB_LB3t, dCT_LB3t, dCP_LB3t, CblExp_LB3t, unbFrac_LB3t]=photoCorrect(tme_LB3t, normintensity_LB3t, alpha, intercept);

[Cnew_LB1c, dCB_LB1c, dCT_LB1c, dCP_LB1c, CblExp_LB1c, unbFrac_LB1c]=photoCorrect(tme_LB1c, normintensity_LB1c, alpha, intercept);
[Cnew_LB20b, dCB_LB20b, dCT_LB20b, dCP_LB20b, CblExp_LB20b, unbFrac_LB20b]=photoCorrect(tme_LB20b, normintensity_LB20b, alpha, intercept);
%% plot the corrected traces for the control
figure, hold on
plot(tme_LB1s, Cnew_LB1s, '-r');
plot(tme_LB2s, Cnew_LB2s, '-b');
plot(tme_LB3s, Cnew_LB3s, '-g'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'control_corrected.png')

figure, hold on
plot(tme_LB1t, Cnew_LB1t, '-r');
plot(tme_LB2t, Cnew_LB2t, '-b');
plot(tme_LB3t, Cnew_LB3t, '-g'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
ciplot(mean(Cnew_LB1s, 1, 'omitnan')-std(Cnew_LB1s, 0, 1, 'omitnan'), mean(Cnew_LB1s, 1, 'omitnan')+std(Cnew_LB1s, 0, 1, 'omitnan'), tme_LB1s, colorcode2{1}, transparency)
plot(tme_LB1s, mean(Cnew_LB1s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(Cnew_LB2s, 1, 'omitnan')-std(Cnew_LB2s, 0, 1, 'omitnan'), mean(Cnew_LB2s, 1, 'omitnan')+std(Cnew_LB2s, 0, 1, 'omitnan'), tme_LB2s, colorcode2{2}, transparency)
plot(tme_LB2s, mean(Cnew_LB2s, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)

ciplot(mean(Cnew_LB3s, 1, 'omitnan')-std(Cnew_LB3s, 0, 1, 'omitnan'), mean(Cnew_LB3s, 1, 'omitnan')+std(Cnew_LB3s, 0, 1, 'omitnan'), tme_LB3s, colorcode2{3}, transparency)
plot(tme_LB3s, mean(Cnew_LB3s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

ciplot(mean(Cnew_LB1t, 1, 'omitnan')-std(Cnew_LB1t, 0, 1, 'omitnan'), mean(Cnew_LB1t, 1, 'omitnan')+std(Cnew_LB1t, 0, 1, 'omitnan'), tme_LB1t, colorcode2{4}, transparency)
plot(tme_LB1t, mean(Cnew_LB1t, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

ciplot(mean(Cnew_LB2t, 1, 'omitnan')-std(Cnew_LB2t, 0, 1, 'omitnan'), mean(Cnew_LB2t, 1, 'omitnan')+std(Cnew_LB2t, 0, 1, 'omitnan'), tme_LB2t, colorcode2{5}, transparency)
plot(tme_LB2t, mean(Cnew_LB2t, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(Cnew_LB3t, 1, 'omitnan')-std(Cnew_LB3t, 0, 1, 'omitnan'), mean(Cnew_LB3t, 1, 'omitnan')+std(Cnew_LB3t, 0, 1, 'omitnan'), tme_LB3t, colorcode2{6}, transparency)
plot(tme_LB3t, mean(Cnew_LB3t, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
%% plot the LB frame rate experiments
cd([dirsave '/LB'])

figure, hold on
plot(time_LB1a, intensity_LB1a, 'Color', colorcode{1});
plot(time_LB1b, intensity_LB1b, 'Color', colorcode{3});
plot(time_LB5a, intensity_LB5a, 'Color', colorcode{4}); 
plot(time_LB10a, intensity_LB10a, 'Color', colorcode{6}); 
plot(time_LB20a, intensity_LB20a, 'Color', colorcode{5}); 
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
%saveas(gcf, 'frameRate_intensity.png');

figure, hold on
plot(time_LB1a, adjintensity_LB1a, 'Color', colorcode{1});
plot(time_LB1b, adjintensity_LB1b, 'Color', colorcode{3});
plot(time_LB5a, adjintensity_LB5a, 'Color', colorcode{4});
plot(time_LB10a, adjintensity_LB10a, 'Color', colorcode{6}); 
plot(time_LB20a, adjintensity_LB20a,  'Color', colorcode{5}); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
%saveas(gcf, 'frameRate_adjintensity.png');

figure, hold on
plot(tme_LB1a, normintensity_LB1a, 'Color', colorcode{1});
plot(tme_LB1b, normintensity_LB1b, 'Color', colorcode{3});
plot(tme_LB5a, normintensity_LB5a, 'Color', colorcode{4});
plot(tme_LB10a, normintensity_LB10a, 'Color', colorcode{6}); 
plot(tme_LB20a, normintensity_LB20a, 'Color', colorcode{5}); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'frameRate_normintensity.png');

figure, hold on
plot(tme_LB1a, Cnew_LB1a, 'Color', colorcode{1});
plot(tme_LB1b, Cnew_LB1b, 'Color', colorcode{3});
plot(tme_LB5a, Cnew_LB5a, 'Color', colorcode{4});
plot(tme_LB10a, Cnew_LB10a, 'Color', colorcode{6}); 
plot(tme_LB20a, Cnew_LB20a, 'Color', colorcode{5}); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'frameRate_corrected.png');

figure, hold on
%ciplot(mean(Cnew_LB1a, 1, 'omitnan')-std(Cnew_LB1a, 0, 1, 'omitnan'), mean(Cnew_LB1a, 1, 'omitnan')+std(Cnew_LB1a, 0, 1, 'omitnan'), tme_LB1a, colorcode2{1}, transparency)
plot(tme_LB1a, mean(Cnew_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(Cnew_LB1b, 1, 'omitnan')-std(Cnew_LB1b, 0, 1, 'omitnan'), mean(Cnew_LB1b, 1, 'omitnan')+std(Cnew_LB1b, 0, 1, 'omitnan'), tme_LB1b, colorcode2{3}, transparency)
plot(tme_LB1b, mean(Cnew_LB1b, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

%ciplot(mean(Cnew_LB5a, 1, 'omitnan')-std(Cnew_LB5a, 0, 1, 'omitnan'), mean(Cnew_LB5a, 1, 'omitnan')+std(Cnew_LB5a, 0, 1, 'omitnan'), tme_LB5a, colorcode2{4}, transparency)
plot(tme_LB5a, mean(Cnew_LB5a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

%ciplot(mean(Cnew_LB10a, 1, 'omitnan')-std(Cnew_LB10a, 0, 1, 'omitnan'), mean(Cnew_LB10a, 1, 'omitnan')+std(Cnew_LB10a, 0, 1, 'omitnan'), tme_LB10a, colorcode2{5}, transparency)
plot(tme_LB10a, mean(Cnew_LB10a, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(Cnew_LB20a, 1, 'omitnan')-std(Cnew_LB20a, 0, 1, 'omitnan'), mean(Cnew_LB20a, 1, 'omitnan')+std(Cnew_LB20a, 0, 1, 'omitnan'), tme_LB20a, colorcode2{6}, transparency)
plot(tme_LB20a, mean(Cnew_LB20a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
legend({'LB 1 min, rep 1', 'LB 1 min, rep 2', 'LB 5 min', 'LB 10 min', 'LB 20 min'})
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'frameRate_correctedAvg.png');

%% find where the mean traces plateau
[avg_LB1a, imend_LB1a]=traceAvg(Cnew_LB1a);
[avg_LB1b, imend_LB1b]=traceAvg(Cnew_LB1b);
[avg_LB5a, imend_LB5a]=traceAvg(Cnew_LB5a);
[avg_LB10a, imend_LB10a]=traceAvg(Cnew_LB10a);
[avg_LB20a, imend_LB20a]=traceAvg(Cnew_LB20a);

%combine to see the entire time trace
tme_LB = union(tme_LB1a(1:imend_LB1a), tme_LB1b(1:imend_LB1b));
tme_LB = union(tme_LB, tme_LB5a(1:imend_LB5a));
tme_LB = union(tme_LB, tme_LB10a(1:imend_LB10a));
tme_LB = union(tme_LB, tme_LB20a(1:imend_LB20a));

%create a new matrix to combine the mean traces
avg_LB = nan(5, length(tme_LB));

[~, ia, ib] = intersect(tme_LB1a(1:imend_LB1a), tme_LB);
avg_LB(1, ib)=avg_LB1a(1, ia);

[~, ia, ib] = intersect(tme_LB1b(1:imend_LB1b), tme_LB);
avg_LB(2, ib)=avg_LB1b(1, ia);

[~, ia, ib] = intersect(tme_LB5a(1:imend_LB5a), tme_LB);
avg_LB(3, ib)=avg_LB5a(1, ia);

[~, ia, ib] = intersect(tme_LB10a(1:imend_LB10a), tme_LB);
avg_LB(4, ib)=avg_LB10a(1, ia);

[~, ia, ib] = intersect(tme_LB20a(1:imend_LB20a), tme_LB);
avg_LB(5, ib)=avg_LB20a(1, ia);

%now calculate the standard error of the mean
smple_sz =  sum(~isnan(avg_LB), 1);
ste_LB = std(avg_LB, 0, 1, 'omitnan')./sqrt(smple_sz);

%% plot the mean and ste of the mean for the corrected untreated LB trace
figure, hold on
ciplot(mean(avg_LB, 1, 'omitnan')-ste_LB, mean(avg_LB, 1, 'omitnan')+ste_LB, tme_LB, colorcode2{1}, transparency)
plot(tme_LB, mean(avg_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'frameRate_concatAvg.png');

%% fit the corrected untreated average
[tau_LB, yhat_LB]=tauCalc(tme_LB, mean(avg_LB, 1, 'omitnan'), 1);

%find when fluor goes to zero
xpred = 0:420;
modelfun = @(tau, x)exp(-x./tau);
ypred = modelfun(tau_LB, xpred);

%% plot to see the fit
figure, hold on
plot(tme_LB, mean(avg_LB, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
%plot(tme_LB, yhat_LB, '--k', 'LineWidth', 1)
plot(xpred, ypred, '--m', 'LineWidth', 1)
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'tauFit_LBavg.png')

%% plot the LB treatment experiments
figure, hold on
plot(time_LBMga, intensity_LBMga, 'Color', colorcode{1});
plot(time_LBEa, intensity_LBEa, 'Color', colorcode{3});
plot(time_LBtuna, intensity_LBtuna,'Color', colorcode{4}); 
plot(time_LBvana, intensity_LBvana, 'Color', colorcode{5}); 
plot(time_LBspnta, intensity_LBspnta, 'Color', colorcode{6}); 
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
%saveas(gcf, 'treated_intensity.png');

figure, hold on
plot(time_LBMga, adjintensity_LBMga, 'Color', colorcode{1});
plot(time_LBEa, adjintensity_LBEa, 'Color', colorcode{3});
plot(time_LBtuna, adjintensity_LBtuna, 'Color', colorcode{4}); 
plot(time_LBvana, adjintensity_LBvana, 'Color', colorcode{5}); 
plot(time_LBspnta, adjintensity_LBspnta, 'Color', colorcode{6}); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
%saveas(gcf, 'treated_adjintensity.png');

figure, hold on
plot(tme_LBMga, normintensity_LBMga,'Color', colorcode{1});
plot(tme_LBEa, normintensity_LBEa, 'Color', colorcode{3});
plot(tme_LBtuna, normintensity_LBtuna, 'Color', colorcode{4}); 
plot(tme_LBvana, normintensity_LBvana, 'Color', colorcode{5}); 
plot(tme_LBspnta, normintensity_LBspnta, 'Color', colorcode{6}); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'treated_normintensity.png');

figure, hold on
plot(tme_LBMga, Cnew_LBMga,'Color', colorcode{1});
plot(tme_LBEa, Cnew_LBEa, 'Color', colorcode{3});
plot(tme_LBtuna, Cnew_LBtuna, 'Color', colorcode{4}); 
plot(tme_LBvana, Cnew_LBvana, 'Color', colorcode{5}); 
plot(tme_LBspnta, Cnew_LBspnta, 'Color', colorcode{6}); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'treated_corrected.png');

%% find where the mean traces plateau
[avg_LBMga, imend_LBMga]=traceAvg(Cnew_LBMga(:, 1:end-3));
[avg_LBEa, imend_LBEa]=traceAvg(Cnew_LBEa(:, 1:end-3));
[avg_LBtuna, imend_LBtuna]=traceAvg(Cnew_LBtuna);
[avg_LBvana, imend_LBvana]=traceAvg(Cnew_LBvana);
[avg_LBspnta, imend_LBspnta]=traceAvg(Cnew_LBspnta);

%% plot the mean traces
figure, hold on
plot(tme_LBMga(1:imend_LBMga), avg_LBMga, 'Color', colorcode{1})
plot(tme_LBEa(1:imend_LBEa), avg_LBEa, 'Color', colorcode{3})
plot(tme_LBtuna(1:imend_LBtuna), avg_LBtuna, 'Color', colorcode{4})
plot(tme_LBvana(1:imend_LBvana), avg_LBvana, 'Color', colorcode{6})
plot(tme_LBspnta(1:imend_LBspnta), avg_LBspnta, 'Color', colorcode{8})
ylim([0 1.1])

%% compare the spnt to the 10 min control
figure, hold on
plot(tme_LB10a(1:imend_LB10a), avg_LB10a, 'Color', colorcode{6})
plot(tme_LBspnta(1:imend_LBspnta), avg_LBspnta, 'Color', colorcode{8})
ylim([0 1.1])

%% plot the PBS experiments
cd([dirsave '/PBS'])

figure, hold on
plot(time_PBS2a, intensity_PBS2a, 'Color', colorcode{1});
plot(time_PBS2b, intensity_PBS2b, 'Color', colorcode{3});
plot(time_PBS20a, intensity_PBS20a, 'Color', colorcode{4}); 
plot(time_PBS20b, intensity_PBS20b, 'Color', colorcode{5}); 
plot(time_PBS60a, intensity_PBS60a, 'Color', colorcode{6}); 
%plot(time_PBS60b, intensity_PBS60b, 'Color', colorcode{7}); 
plot(time_PBS120a, intensity_PBS120a, 'Color', colorcode{8}); 
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
%saveas(gcf, 'PBS_intensity.png');

figure, hold on
plot(time_PBS2a, adjintensity_PBS2a, 'Color', colorcode{1});
plot(time_PBS2b, adjintensity_PBS2b, 'Color', colorcode{3});
plot(time_PBS20a, adjintensity_PBS20a, 'Color', colorcode{4}); 
plot(time_PBS20b, adjintensity_PBS20b, 'Color', colorcode{5}); 
plot(time_PBS60a, adjintensity_PBS60a, 'Color', colorcode{6}); 
%plot(time_PBS60b, adjintensity_PBS60b, 'Color', colorcode{7}); 
plot(time_PBS120a, adjintensity_PBS120a, 'Color', colorcode{8}); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
%saveas(gcf, 'PBS_adjintensity.png');

figure, hold on
plot(tme_PBS2a, normintensity_PBS2a, 'Color', colorcode{1});
plot(tme_PBS2b, normintensity_PBS2b, 'Color', colorcode{3});
plot(tme_PBS20a, normintensity_PBS20a, 'Color', colorcode{4}); 
plot(tme_PBS20b, normintensity_PBS20b, 'Color', colorcode{5}); 
plot(tme_PBS60a, normintensity_PBS60a, 'Color', colorcode{6}); 
%plot(tme_PBS60b, normintensity_PBS60b, 'Color', colorcode{7}); 
plot(tme_PBS120a, normintensity_PBS120a, 'Color', colorcode{8}); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'PBS_normintensity.png');

figure, hold on
plot(tme_PBS2a,Cnew_PBS2a, 'Color', colorcode{1});
plot(tme_PBS2b, Cnew_PBS2b, 'Color', colorcode{3});
plot(tme_PBS20a, Cnew_PBS20a, 'Color', colorcode{4}); 
plot(tme_PBS20b, Cnew_PBS20b, 'Color', colorcode{5}); 
plot(tme_PBS60a, Cnew_PBS60a, 'Color', colorcode{6}); 
%plot(tme_PBS60b, Cnew_PBS60b, 'Color', colorcode{7}); 
plot(tme_PBS120a, Cnew_PBS120a, 'Color', colorcode{8}); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
%saveas(gcf, 'PBS_corrected.png');

%% find where the mean traces plateau
[avg_PBS2a, imend_PBS2a]=traceAvg(Cnew_PBS2a);
[avg_PBS2b, imend_PBS2b]=traceAvg(Cnew_PBS2b);
[avg_PBS20a, imend_PBS20a]=traceAvg(Cnew_PBS20a);
[avg_PBS20b, imend_PBS20b]=traceAvg(Cnew_PBS20b);
[avg_PBS60a, imend_PBS60a]=traceAvg(Cnew_PBS60a);
[avg_PBS120a, imend_PBS120a]=traceAvg(Cnew_PBS120a);

%% plot the mean traces
figure, hold on
plot(tme_PBS2a(1:imend_PBS2a), avg_PBS2a, 'Color', colorcode{1})
plot(tme_PBS2b(1:imend_PBS2b), avg_PBS2b, 'Color', colorcode{3})
plot(tme_PBS20a(1:imend_PBS20a), avg_PBS20a, 'Color', colorcode{4})
plot(tme_PBS20b(1:imend_PBS20b), avg_PBS20b, 'Color', colorcode{5})
plot(tme_PBS60a(1:imend_PBS60a), avg_PBS60a, 'Color', colorcode{6})
plot(tme_PBS120a(1:imend_PBS120a), avg_PBS120a, 'Color', colorcode{8})
ylim([0 1.1])

%% predict what the 20% intensity trace will look like for the 1 minute control vs 20 minute control
alpha=152.6455;
intercept=-1.0976;
tau=91.8894;

[Pnew1]=tracePredict([1:1:120], tau, alpha, intercept);
[Pnew2]=tracePredict([1:20:120], tau, alpha, intercept);

%compare
figure,hold on, plot([1:1:120], Pnew1, '-r'), plot([1:20:120], Pnew2, '-b')
%% verify the function
alpha=32.2114;
intercept=0.1614;

[Pnew3]=tracePredict(tme_LB1a, tau, alpha, intercept);

%compare to actual plot
figure, hold on
plot(tme_LB1a, mean(normintensity_LB1a, 1, 'omitnan'), '-r')
plot(tme_LB1a, Pnew3, '-b')
%% Functions
%to aggregate data and normalize
function [intensity, adjintensity, normintensity, lCell, time, tme, imstart]=dataNormalize(datadir, imstart)
        
    %pre-allocate variables
    intensity=[];
    lCell=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, 'icell_intensity', 'time', 'lcell')
        
        [nrow, ~]=size(icell_intensity);
        
        for n=1:nrow
            if ~isnan(icell_intensity(n, imstart))
                intensity=[intensity; icell_intensity(n, :)];
                lCell=[lCell; lcell(n, :)];
            end
        end

        if i==1
            tme=time; %pre-set new time vector
        end

    end
      
    %remove cells without a value for imstart
%     idx=find(~isnan(intensity(:, imstart)));
%     intensity=intensity(idx,:);
%     lCell=lCell(idx,:);    
       
    %subtract final fluor. value from the trace
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<0)=NaN;      

%     %interpolate the fluor values during detergent perfusion
%     for n=1:height(adjintensity)
%         idx=find(adjintensity(n,:)>adjintensity(n, imstart));
%             if ~isempty(idx)
%                 [nrow, ncol]=size(adjintensity);
%                 x=setdiff(1:ncol, idx); %x=time
%                 v=adjintensity(n, x); %v=intensity
%                 vq=interp1(x, v, idx); %vq=interpolation at query time pts
%                 adjintensity(n,idx)=vq;
%             end
%     end
        
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:end)./adjintensity(:, imstart);
    [nrow, ncol]=size(normintensity);
    
    %interpolate the fluor values during detergent perfusion (& other noisy
    %spike in fluor.)
    for n=1:height(normintensity)
        
        dv=diff(normintensity(n,:), 1, 2);
        idx1=find(dv>0)+1;
        idx2=find(normintensity(n,:)>1);
        idx=union(idx1, idx2);
        
            if ~isempty(idx)
                x=setdiff(1:ncol, idx); %x=time
                v=normintensity(n, x); %v=intensity
                vq=interp1(x, v, idx); %vq=interpolation at query time pts
                normintensity(n,idx)=vq;
            end
    end
    
    %adjust the time vector
    tme=tme(imstart:end)-tme(imstart);     
        
end

%to aggregate control data and normalize
function [intensity, bgintensity, adjintensity, normintensity, lCell, time, tme, imstart]=controlNormalize(datadir)

    %pre-allocate variables
    intensity=[];
    lCell=[];
    bgintensity=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, 'icell_intensity', 'bg_intensity', 'time', 'lcell')
        
        intensity=[intensity; icell_intensity];
        lCell=[lCell; lcell];
        bgintensity=[bgintensity; bg_intensity];

        if i==1
            tme=time; %pre-set new time vector
        end

    end
        

    %find the initial post-lysis frame by identifying the first
    %time point where dt matches the final dt
    dt=round(diff(time));
    if dt(1)==dt(end)
%         dl=diff(lCell, 1);
%         lvg=mean(dl, 1, 'omitnan');
%         imstart=find(lvg<0, 1, 'first')+2;
        imstart=1;
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        intensity=intensity(idx,:);
        lCell=lCell(idx,:);    

        %subtract final fluor. value from the trace and set limit of
        %detection
        adjintensity=intensity-bgintensity; 
        adjintensity(adjintensity<200)=NaN; 
    
    else
        if dt(end)<1 %discrepancies in dt might be found in the meta data, this is the ad hoc fix
            imstart=find(dt==dt(end), 1, 'first')+2;
        else 
            imstart=find(dt==dt(end), 1, 'first');
        end
    
  
    %remove cells without a value for imstart
    idx=find(~isnan(intensity(:, imstart)));
    intensity=intensity(idx,:);
    lCell=lCell(idx,:);    
        
    %subtract final fluor. value from the trace and set limit of
    %detection
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<200)=NaN; 
    
    end
        
    %find when all values are below the limit of detection
    [~, ncol]=size(adjintensity);
    nsum=sum(isnan(adjintensity));
    if max(nsum)==ncol
        imend=min(find(nsum==ncol));
    else 
        imend=ncol;
    end
    
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:imend)./adjintensity(:, imstart);

    %adjust the time vector
    tme=tme(imstart:imend)-tme(imstart);     
        
end

%do not fit to y=(1-beta)*exp(-t./tau)+beta (beta=normalized limit of detection)
%fit to y=Ae^(-t/tau)
function [tau, yhat]=tauCalc(tme, normintensity, tau0)
    
    [nrow, ~]=size(normintensity);
    tau=nan(nrow,1);
    yhat=nan(size(normintensity));
    modelfun = @(tau, x)exp(-x./tau);
    for i=1:nrow
        tau(i, 1)=nlinfit(tme, normintensity(i,:), modelfun, tau0);
         yhat(i, :)=modelfun(tau(i), tme);
    end
    
end

%to correct for photobleaching
function [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept)
        
        %pre-allocate variables
        %assume that the initial 'measured' fluorescence values and corrected
        %fluor. values will be equal. I prefer to pre-allocate with nan in case 
        %some values are missing in the raw data
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:,1);
             
        dCB=nan(height(normintensity), length(tme)-1); %change in fluor. due to photobleaching
        dCT=nan(height(normintensity), length(tme)-1); %total change in fluor.
        dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability

        unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
        unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

        Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
        Cbl_exp(:, 1)=0;
        
        %calculate dt (the dt between frames may vary in a single run).
        %Note: this is also the frame rate
        dt=round(diff(tme), 2);
        
        %this formula comes from the slope and intercept calculated for the 1.2, 2,
        %and 3 second tau vs frame rate controls 
       dC=@(C, alpha, dt, b)(C/(alpha*dt+b))*dt;
       %dC=@(C, alpha, intercept, dt, b)(C/(alpha*dt+b));
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                dCB(n,i) = dC(normintensity(n,i), alpha, dt(i), intercept); %this is the amount of photobleaching that occured in our measured value

                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)*unb_frac(n,i); %Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
        end      

end

%to find the time point at which the corrected trace plateaus and calculate
%the average
function [avg, imend]=traceAvg(Cnew)
    cmean=mean(Cnew, 1, 'omitnan');
    imend=find(cmean-cmean(end)<0.01, 1, 'first');
    if isempty(imend)
        lidx=find(cmean>=0, 1, 'last');
        imend=find(cmean-cmean(lidx)<0.01, 1, 'first');
    end
    avg=cmean(:, 1:imend);
end

%to calculate strain
function [strain]=strainCalc(lcell, imstart)

    lcell=lcell(:, 1:imstart+3);
%     dl=diff(lcell, 1, 2);
% 
%     strain=dl./lcell(:, 1:end-1);

    [nrow, ncol]=size(lcell);
    strain=nan(nrow, ncol);

    for i=1:ncol-1
        dl=lcell(:,i)-lcell(:, i+1);
        strain(:, i)=dl./lcell(:, i+1);
    end

end

%to calculation growth rate
function [growthRate]=gRateCalc(time, lcell, imstart)

    lcell=lcell(:, 1:imstart+3);
    time=time(1:imstart+3);

    dl=diff(lcell, 1, 2);
    dt=diff(time, 1, 2);

    sl=(lcell(:, 1:end-1)+(dl+lcell(:, 1:end-1)))/2;

    growthRate=(dl./dt)./sl;

end

%to predict a normintensity trace
function [Pnew]=tracePredict(time, tau, alpha, intercept)
    
    diffusion=@(time, tau)exp(-time/tau);
    trueTrace=diffusion(time, tau);
    Pnew=trueTrace;
    
    dt=diff(time, 1, 2);
    dC=@(C, alpha, dt, b)(C/(alpha*dt+b))*dt;
    unb_frac=1;
    
    for i=1:length(time)-1
        dCP=trueTrace(1, i+1)-trueTrace(1,i);
        dCP=dCP*unb_frac;
     
        dCB = dC(trueTrace(1,i), alpha, dt(i), intercept);
        Pnew(1, i+1) = Pnew(1, i) - dCB;
        
        Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
        unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
    
    end

end