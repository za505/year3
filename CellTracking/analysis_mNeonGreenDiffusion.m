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
LB5a=dir(['02122022_Exp1' '*dm.mat']); %LB, fr = 5 min
LB20a=dir(['02092022_Exp1' '*dm.mat']); %LB, fr = 20 min

PBS2a=dir(['11192021_Exp2' '*dm.mat']); %PBS 2 min, rep 1
PBS2b=dir(['12082021_Exp3' '*dm.mat']); %PBS 2 min, rep 2
PBS20a=dir(['11192021_Exp1' '*dm.mat']); %PBS 20 min, rep 1
PBS20b=dir(['11302021_Exp1' '*dm.mat']); %PBS 20 min, rep 2
PBS60a=dir(['10232021_Exp2' '*dm.mat']); %PBS 60 min, rep 1
PBS60b=dir(['10262021_Exp2' '*dm.mat']); %PBS 60 min, rep 2
PBS120a=dir(['01142022_Exp1' '*dm.mat']); %PBS 120 min

LBMga=dir(['01172022_Exp1' '*dm.mat']); %LB, 20 mM Mg2+
LBEa=dir(['01172022_Exp2' '*dm.mat']); %LB, 10 mM EDTA
LBtuna=dir(['01242022_Exp1' '*dm.mat']); %LB, 0.5 ug/mL tunicamycin
LBvana=dir(['01262022_Exp1' '*dm.mat']); %LB, 1 ug/mL vancomycin

LB1s=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
LB2s=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
LB3s=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s

%% calculate normalized fluorescence traces
[intensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lcell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize(LB1a, 6); 
[intensity_LB1b, adjintensity_LB1b, normintensity_LB1b, lcell_LB1b, time_LB1b, tme_LB1b, imstart_LB1b]=dataNormalize(LB1b, 5); 
[intensity_LB5a, adjintensity_LB5a, normintensity_LB5a, lcell_LB5a, time_LB5a, tme_LB5a, imstart_LB5a]=dataNormalize(LB5a, 5); 
[intensity_LB20a, adjintensity_LB20a, normintensity_LB20a, lcell_LB20a, time_LB20a, tme_LB20a, imstart_LB20a]=dataNormalize(LB20a, 6); 

[intensity_PBS2a, adjintensity_PBS2a, normintensity_PBS2a, lcell_PBS2a, time_PBS2a, tme_PBS2a, imstart_PBS2a]=dataNormalize(PBS2a, 10); 
[intensity_PBS2b, adjintensity_PBS2b, normintensity_PBS2b, lcell_PBS2b, time_PBS2b, tme_PBS2b, imstart_PBS2b]=dataNormalize(PBS2b, 10); 
[intensity_PBS20a, adjintensity_PBS20a, normintensity_PBS20a, lcell_PBS20a, time_PBS20a, tme_PBS20a, imstart_PBS20a]=dataNormalize(PBS20a, 4); 
[intensity_PBS20b, adjintensity_PBS20b, normintensity_PBS20b, lcell_PBS20b, time_PBS20b, tme_PBS20b, imstart_PBS20b]=dataNormalize(PBS20b, 4); 
[intensity_PBS60a, adjintensity_PBS60a, normintensity_PBS60a, lcell_PBS60a, time_PBS60a, tme_PBS60a, imstart_PBS60a]=dataNormalize(PBS60a, 3); 
[intensity_PBS60b, adjintensity_PBS60b, normintensity_PBS60b, lcell_PBS60b, time_PBS60b, tme_PBS60b, imstart_PBS60b]=dataNormalize(PBS60b, 5); 
[intensity_PBS120a, adjintensity_PBS120a, normintensity_PBS120a, lcell_PBS120a, time_PBS120a, tme_PBS120a, imstart_PBS120a]=dataNormalize(PBS120a, 9); 

[intensity_LBMga, adjintensity_LBMga, normintensity_LBMga, lcell_LBMga, time_LBMga, tme_LBMga, imstart_LBMga]=dataNormalize(LBMga, 6); 
[intensity_LBEa, adjintensity_LBEa, normintensity_LBEa, lcell_LBEa, time_LBEa, tme_LBEa, imstart_LBEa]=dataNormalize(LBEa, 6); 
[intensity_LBtuna, adjintensity_LBtuna, normintensity_LBtuna, lcell_LBtuna, time_LBtuna, tme_LBtuna, imstart_LBtuna]=dataNormalize(LBtuna, 3); 
[intensity_LBvana, adjintensity_LBvana, normintensity_LBvana, lcell_LBvana, time_LBvana, tme_LBvana, imstart_LBvana]=dataNormalize(LBvana, 5); %the true pre-lysis frame is t=4, but there are no fluor readings 

[intensity_LB1s, bgintensity_LB1s, adjintensity_LB1s, normintensity_LB1s, lcell_LB1s, time_LB1s, tme_LB1s, imstart_LB1s]=controlNormalize(LB1s); 
[intensity_LB2s, bgintensity_LB2s, adjintensity_LB2s, normintensity_LB2s, lcell_LB2s, time_LB2s, tme_LB2s, imstart_LB2s]=controlNormalize(LB2s); 
[intensity_LB3s, bgintensity_LB3s, adjintensity_LB3s, normintensity_LB3s, lcell_LB3s, time_LB3s, tme_LB3s, imstart_LB3s]=controlNormalize(LB3s); 

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
saveas(gcf, 'control_intensity.png');

% plot the adjusted intensity
figure, hold on
plot(time_LB1s, adjintensity_LB1s, '-r');
plot(time_LB2s, adjintensity_LB2s, '-b');
plot(time_LB3s, adjintensity_LB3s, '-g'); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
saveas(gcf, 'control_adjintensity.png');

% plot the normalized intensity
figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r');
plot(tme_LB2s, normintensity_LB2s, '-b');
plot(tme_LB3s, normintensity_LB3s, '-g'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'control_normintensity.png');

%% calculate alpha
[tau1, yhat_LB1s]=tauCalc(tme_LB1s, normintensity_LB1s);
[tau2, yhat_LB2s]=tauCalc(tme_LB2s, normintensity_LB2s);
[tau3, yhat_LB3s]=tauCalc(tme_LB3s, normintensity_LB3s);

%% plots for tau
figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r')
plot(tme_LB1s, yhat_LB1s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'tauFit_LB1s.png')

figure, hold on
plot(tme_LB2s, normintensity_LB2s, '-b')
plot(tme_LB2s, yhat_LB2s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'tauFit_LB2s.png')

figure, hold on
plot(tme_LB3s, normintensity_LB3s, '-g')
plot(tme_LB3s, yhat_LB3s, '--k')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'tauFit_LB3s.png')

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
saveas(gcf, 'tau_vs_frameRate.png')

%% correct for photobleaching
alpha=32.2114;
intercept=0.1614;

[Cnew_LB1a, dCB_LB1a, dCT_LB1a, dCP_LB1a, CblExp_LB1a, unbFrac_LB1a]=photoCorrect(tme_LB1a, normintensity_LB1a, alpha, intercept);
[Cnew_LB1b, dCB_LB1b, dCT_LB1b, dCP_LB1b, CblExp_LB1b, unbFrac_LB1b]=photoCorrect(tme_LB1b, normintensity_LB1b, alpha, intercept);
[Cnew_LB5a, dCB_LB5a, dCT_LB5a, dCP_LB5a, CblExp_LB5a, unbFrac_LB5a]=photoCorrect(tme_LB5a, normintensity_LB5a, alpha, intercept);
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

[Cnew_LB1s, dCB_LB1s, dCT_LB1s, dCP_LB1s, CblExp_LB1s, unbFrac_LB1s]=photoCorrect(tme_LB1s, normintensity_LB1s, alpha, intercept);
[Cnew_LB2s, dCB_LB2s, dCT_LB2s, dCP_LB2s, CblExp_LB2s, unbFrac_LB2s]=photoCorrect(tme_LB2s, normintensity_LB2s, alpha, intercept);
[Cnew_LB3s, dCB_LB3s, dCT_LB3s, dCP_LB3s, CblExp_LB3s, unbFrac_LB3s]=photoCorrect(tme_LB3s, normintensity_LB3s, alpha, intercept);

%% plot the corrected traces for the control
figure, hold on
plot(tme_LB1s, Cnew_LB1s, '-r');
plot(tme_LB2s, Cnew_LB2s, '-b');
plot(tme_LB3s, Cnew_LB3s, '-g'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'control_corrected.png')

%% plot the LB frame rate experiments
cd([dirsave '/LB'])

figure, hold on
plot(time_LB1a, intensity_LB1a, 'Color', colorcode{1});
plot(time_LB1b, intensity_LB1b, 'Color', colorcode{3});
plot(time_LB5a, intensity_LB5a, 'Color', colorcode{4}); 
plot(time_LB20a, intensity_LB20a, 'Color', colorcode{5}); 
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
saveas(gcf, 'frameRate_intensity.png');

figure, hold on
plot(time_LB1a, adjintensity_LB1a, 'Color', colorcode{1});
plot(time_LB1b, adjintensity_LB1b, 'Color', colorcode{3});
plot(time_LB5a, adjintensity_LB5a, 'Color', colorcode{4});
plot(time_LB20a, adjintensity_LB20a,  'Color', colorcode{5}); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
saveas(gcf, 'frameRate_adjintensity.png');

figure, hold on
plot(tme_LB1a, normintensity_LB1a, 'Color', colorcode{1});
plot(tme_LB1b, normintensity_LB1b, 'Color', colorcode{3});
plot(tme_LB5a, normintensity_LB5a, 'Color', colorcode{4});
plot(tme_LB20a, normintensity_LB20a, 'Color', colorcode{5}); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'frameRate_normintensity.png');

figure, hold on
plot(tme_LB1a, Cnew_LB1a, 'Color', colorcode{1});
plot(tme_LB1b, Cnew_LB1b, 'Color', colorcode{3});
plot(tme_LB5a, Cnew_LB5a, 'Color', colorcode{4});
plot(tme_LB20a, Cnew_LB20a, 'Color', colorcode{5}); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'frameRate_corrected.png');

%% plot the LB treatment experiments
figure, hold on
plot(time_LBMga, intensity_LBMga, '-r');
plot(time_LBEa, intensity_LBEa, '-b');
plot(time_LBtuna, intensity_LBtuna, '-g'); 
plot(time_LBvana, intensity_LBvana, '-m'); 
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
saveas(gcf, 'treated_intensity.png');

figure, hold on
plot(time_LBMga, adjintensity_LBMga, '-r');
plot(time_LBEa, adjintensity_LBEa, '-b');
plot(time_LBtuna, adjintensity_LBtuna, '-g'); 
plot(time_LBvana, adjintensity_LBvana, '-m'); 
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
saveas(gcf, 'treated_adjintensity.png');

figure, hold on
plot(tme_LBMga, normintensity_LBMga, '-r');
plot(tme_LBEa, normintensity_LBEa, '-b');
plot(tme_LBtuna, normintensity_LBtuna, '-g'); 
plot(tme_LBvana, normintensity_LBvana, '-m'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'treated_normintensity.png');

figure, hold on
plot(tme_LBMga, Cnew_LBMga, '-r');
plot(tme_LBEa, Cnew_LBEa, '-b');
plot(tme_LBtuna, Cnew_LBtuna, '-g'); 
plot(tme_LBvana, Cnew_LBvana, '-m'); 
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
saveas(gcf, 'treated_corrected.png');

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
saveas(gcf, 'PBS_intensity.png');

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
saveas(gcf, 'PBS_adjintensity.png');

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
saveas(gcf, 'PBS_normintensity.png');

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
saveas(gcf, 'PBS_corrected.png');

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

    %interpolate the fluor values during detergent perfusion
    for n=1:height(adjintensity)
        idx=find(adjintensity(n,:)>adjintensity(n, imstart));
            if ~isempty(idx)
                [nrow, ncol]=size(adjintensity);
                x=setdiff(1:ncol, idx); %x=time
                v=adjintensity(n, x); %v=intensity
                vq=interp1(x, v, idx); %vq=interpolation at query time pts
                adjintensity(n,idx)=vq;
            end
    end
        
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:end)./adjintensity(:, imstart);

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
        dl=diff(lCell, 1);
        lvg=mean(dl, 1, 'omitnan');
        imstart=find(lvg<0, 1, 'first')+2;
    else
        if dt(end)<1 %discrepancies in dt might be found in the meta data, this is the ad hoc fix
            imstart=find(dt==dt(end), 1, 'first')+2;
        else 
            imstart=find(dt==dt(end), 1, 'first');
        end
    end
  
    %remove cells without a value for imstart
    idx=find(~isnan(intensity(:, imstart)));
    intensity=intensity(idx,:);
    lCell=lCell(idx,:);    
        
    %subtract final fluor. value from the trace and set limit of
    %detection
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<200)=NaN; 

        
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
function [tau, yhat]=tauCalc(tme, normintensity)
    
    [nrow, ~]=size(normintensity);
    tau=nan(nrow,1);
    yhat=nan(size(normintensity));
    modelfun = @(tau, x)exp(-x./tau);
    for i=1:nrow
        tau(i, 1)=nlinfit(tme, normintensity(i,:), modelfun, 1);
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
