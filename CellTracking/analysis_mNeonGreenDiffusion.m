%Author: Zarina Akbary
%Date: 12/17/2021
%Purpose: To combine and analyze decayMeasure (dm.mat) data from mNeonGreen diffusion
%experiments. 

clear, close all

%Inputs
%dirsave = directory where all the dm.mat files are located and where plots
%are saved

%Outputs
%intensity = m x n matrix of raw intensity values; concatenated
%icell_intensity variable from dm.mat files
%time = 1 x n matrix of time (minutes) values that correspond to
%intensity traces
%normintensity = m x n matrix of intensity values where 1) the background
%was subtracted 2) the intensity values during detergent perfusion were
%interpolated and 3) the intensity values were normalized by the initial
%pre-lysis value
%time = 1 x n matrix of time (minutes) values that correspond to
%normalized intensity traces
%Cnew = m x n matrix of normalized intensity values corrected for
%photobleaching
%beta = m x 1 matrix of final fluorescent values for corrected traces
%dCB = m x n matrix of change in corrected fluorescence attributable to
%photobleaching between two adjacent time points
%dCT = m x n matrix of change in corrected fluorescence between two time points
%dCP = m x n matrix of change in corrected fluorescence attributable to
%permeability between two adjacent time points
%CblExp = m x n matrix of the concentration of bleached fluorophores over
%time
%unbFrac = m x n matrix of the fraction of unbleached fluor over time
%midx = index of Cnew traces that correspond to normintensity traces; indexed
%Cnew traces have a beta value <=0.9 (greater than 0.9 are outliers with
%noisy corrected traces) and are not NaN
%tau = m x 1 matrix of time constants calculated by fitting Cnew to the
%exponential decay function y = (1-beta)*e^-t/tau+beta
%yhat = m x n matrix of predicted trace value
%initial = m x n matrix of initial phase intensities
%final = m x n matrix of final phase intensities
%gamma = m x 1 matrix of final/initial phase intensity ratios

%% Store data directory files

% colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#E872D0'};
% colorcode2={'#A1D5F7', '#FFBB9E', '#FFE3A3', '#EAA6F7', '#CBF09C', '#C7EEFF', '#FFB3C1', '#FACAF0'};

% colorcode={[0 0.45 0.74], [0.85 0.33 0.1], [0.49 0.18 0.56], [0.93 0.69 0.13], [0.47 0.67 0.19], [0.3 0.75 0.93], [0.64 0.08 0.18], [0.91 0.45 0.82]};
% colorcode2={[0.63 0.84 0.97], [1 0.73 0.62], [0.92 0.65 0.97], [0.97 0.91 0.69], [0.8 0.94 0.61], [0.82 0.95 1], [1 0.7 0.76], [0.98 0.79 0.94]};
%blue, orange, purple, yellow, green, light blue, red, pink

%color codes must be in RGB [0-1] format to be used in ciplot
colorcode={[204 0 0], [204 102 0], [204 204 0], [102 204 0], [0 204 204], [0 0 204], [102 0 204], [204 0 204], [204 0 102], [255 102 102], [255 178 102], [102 255 102], [102 255 255], [102 178 255], [178 102 255], [255 102 255],[255 102 178]};
colorcode2={[255 51 51], [255 153 51], [255 255 51], [153 255 51], [51 255 255], [51 51 255], [153 51 255], [255 51 255], [255 51 153], [255 204 204], [255 229 204], [204 255 204], [204 255 255], [204 229 255], [229 204 255], [255 204 255],[255 204 229]};

colorcode=cellfun(@(x)(x./255), colorcode, 'UniformOutput', false);
colorcode2=cellfun(@(x)(x./255), colorcode2, 'UniformOutput', false);

%direct the code to the location of the dm.mat files 
dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/01252022_analysis';
cd([dirsave '/MatFiles'])

LBa=dir(['10232021_Exp1' '*dm.mat']); %LB, rep 1
LBb=dir(['10262021_Exp1' '*dm.mat']); %LB, rep 2
LB5a=dir(['02012022_Exp1' '*dm.mat']); %LB, 5 min frame rate
LB20a=dir(['02042022_Exp1' '*dm.mat']); %LB, 20 min frame rate
PBS120=dir(['01142022_Exp1' '*dm.mat']); %PBS 120 min
PBS60a=dir(['10232021_Exp2' '*dm.mat']); %PBS 60 min, rep 1
PBS60b=dir(['10262021_Exp2' '*dm.mat']); %PBS 60 min, rep 2
PBS20a=dir(['11192021_Exp1' '*dm.mat']); %PBS 20 min, rep 1
PBS20b=dir(['11302021_Exp1' '*dm.mat']); %PBS 20 min, rep 2
PBS2a=dir(['11192021_Exp2' '*dm.mat']); %PBS 2 min, rep 1
PBS2b=dir(['12082021_Exp3' '*dm.mat']); %PBS 2 min, rep 2

LBMga=dir(['01172022_Exp1' '*dm.mat']); %LB, 20 mM Mg2+
LBEa=dir(['01172022_Exp2' '*dm.mat']); %LB, 10 mM EDTA
LBtuna=dir(['01242022_Exp1' '*dm.mat']); %LB, 0.5 ug/mL tunicamycin
LBdeta=dir(['01192022_Exp2' '*dm.mat']); %LB, 5% NLS
LBvanco=dir(['01262022_Exp1' '*dm.mat']); %LB, 1 ug/mL vancomycin

%direct the code to the location of the pc.mat files 
LBc=dir(['10232021_Exp1' '*pc.mat']); %LB, rep 1
LBd=dir(['10262021_Exp1' '*pc.mat']); %LB, rep 2
PBS120c=dir(['01142022_Exp1' '*pc.mat']); %PBS 120 min
PBS60c=dir(['10232021_Exp2' '*pc.mat']); %PBS 60 min, rep 1
PBS60d=dir(['10262021_Exp2' '*pc.mat']); %PBS 60 min, rep 2
PBS20c=dir(['11192021_Exp1' '*pc.mat']); %PBS 20 min, rep 1
PBS20d=dir(['11302021_Exp1' '*pc.mat']); %PBS 20 min, rep 2
PBS2c=dir(['11192021_Exp2' '*pc.mat']); %PBS 2 min, rep 1
PBS2d=dir(['12082021_Exp3' '*pc.mat']); %PBS 2 min, rep 2

LBMgc=dir(['01172022_Exp1' '*pc.mat']); %LB, 20 mM Mg2+
LBEc=dir(['01172022_Exp2' '*pc.mat']); %LB, 10 mM EDTA
LBtunc=dir(['01242022_Exp1' '*pc.mat']); %LB, 0.5 ug/mL tunicamycin
LBdetc=dir(['01192022_Exp2' '*pc.mat']); %LB, 5% NLS
LBvanc=dir(['01262022_Exp1' '*pc.mat']); %LB + 1 ug/mL vancomycin, since the cell
%tracking is done using WT, the alignment with phase is a little off, but not by much. 
%Strangely, there are huge spikes after 40 minutes (around the time there's a sudden expansion in cell length

%% Calculate phase contrast ratio
% [initial_LBc, final_LBc, gamma_LBc]=phaseRatio(LBc, 5);
% [initial_LBd, final_LBd, gamma_LBd]=phaseRatio(LBd, 4);
% [initial_PBS120, final_PBS120, gamma_PBS120]=phaseRatio(PBS120, 4);
% [initial_PBS60c, final_PBS60c, gamma_PBS60c]=phaseRatio(PBS60c, 2);
% [initial_PBS60d, final_PBS60d, gamma_PBS60d]=phaseRatio(PBS60d, 4);
% [initial_PBS20c, final_PBS20c, gamma_PBS20c]=phaseRatio(PBS20c, 3);
% [initial_PBS20d, final_PBS20d, gamma_PBS20d]=phaseRatio(PBS20d, 3);
% [initial_PBS2c, final_PBS2c, gamma_PBS2c]=phaseRatio(PBS2c, 9);
% [initial_PBS2d, final_PBS2d, gamma_PBS2d]=phaseRatio(PBS2d, 9);
% 
% [initial_LBMgc, final_LBMgc, gamma_LBMgc]=phaseRatio(LBMgc, 5);
% [initial_LBEc, final_LBEc, gamma_LBEc]=phaseRatio(LBEc, 5);
% [initial_LBtunc, final_LBtunc, gamma_LBtunc]=phaseRatio(LBtunc, 2);
% [initial_LBdetc, final_LBdetc, gamma_LBdetc]=phaseRatio(LBdetc, 5);
% [initial_LBvanc, final_LBvanc, gamma_LBvanc]=phaseRatio(LBvanc, 5, 49); %the numbers get too low at the end so truncating it is the best way to get data that makes sense

%% Calculate normalized fluorescence traces
[normintensity_LBa, intensity_LBa, time_LBa, tme_LBa, imend_LBa]=dataNormalize(LBa, 5);
[normintensity_LBb, intensity_LBb, time_LBb, tme_LBb, imend_LBb]=dataNormalize(LBb, 4);
[normintensity_LB5a, intensity_LB5a, time_LB5a, tme_LB5a, imend_LB5a]=dataNormalize(LB5a, 7);
[normintensity_LB20a, intensity_LB20a, time_LB20a, tme_LB20a, imend_LB20a]=dataNormalize(LB20a, 9);
[normintensity_PBS120, intensity_PBS120, time_PBS120, tme_PBS120, imend_PBS120]=dataNormalize(PBS120, 4);
[normintensity_PBS60a, intensity_PBS60a, time_PBS60a, tme_PBS60a, imend_PBS60a]=dataNormalize(PBS60a, 2);
[normintensity_PBS60b, intensity_PBS60b, time_PBS60b, tme_PBS60b, imend_PBS60b]=dataNormalize(PBS60b, 4);
[normintensity_PBS20a, intensity_PBS20a, time_PBS20a, tme_PBS20a, imend_PBS20a]=dataNormalize(PBS20a, 3);
[normintensity_PBS20b, intensity_PBS20b, time_PBS20b, tme_PBS20b, imend_PBS20b]=dataNormalize(PBS20b, 3);
[normintensity_PBS2a, intensity_PBS2a, time_PBS2a, tme_PBS2a, imend_PBS2a]=dataNormalize(PBS2a, 9);
[normintensity_PBS2b, intensity_PBS2b, time_PBS2b, tme_PBS2b, imend_PBS2b]=dataNormalize(PBS2b, 9);

[normintensity_LBMga, intensity_LBMga, time_LBMga, tme_LBMga, imend_LBMga]=dataNormalize(LBMga, 5);
[normintensity_LBEa, intensity_LBEa, time_LBEa, tme_LBEa, imend_LBEa]=dataNormalize(LBEa, 5);
[normintensity_LBtuna, intensity_LBtuna, time_LBtuna, tme_LBtuna, imend_LBtuna]=dataNormalize(LBtuna, 2);
[normintensity_LBdeta, intensity_LBdeta, time_LBdeta, tme_LBdeta, imend_LBdeta]=dataNormalize(LBdeta, 5);
[normintensity_LBvanco, intensity_LBvanco, time_LBvanco, tme_LBvanco, imend_LBvanco]=dataNormalize(LBvanco, 5);

%% combine datasets
time_LB=time_LBb;
time_PBS60=time_PBS60a;
time_PBS20=time_PBS20a;
time_PBS2=time_PBS2a;

tme_LB=tme_LBa;
tme_PBS60=tme_PBS60b;
tme_PBS20=tme_PBS20a;
tme_PBS2=tme_PBS2b;

imend_LB=imend_LBa;
imend_PBS60=imend_PBS60b;
imend_PBS20=imend_PBS20a;
imend_PBS2=imend_PBS2b;

intensity_LB=[intensity_LBa(:, 1:length(time_LB)); intensity_LBb];
intensity_PBS60=[intensity_PBS60a; intensity_PBS60b(:, 1:length(time_PBS60))];
intensity_PBS20=[intensity_PBS20a; intensity_PBS20b(:, 1:length(time_PBS20))];
intensity_PBS2=[intensity_PBS2a; intensity_PBS2b(:, 1:length(time_PBS2))];

normintensity_LB=[normintensity_LBa; normintensity_LBb(:, 1:length(tme_LBa))];
normintensity_PBS60=[normintensity_PBS60a(:, 1:length(tme_PBS60)); normintensity_PBS60b];
normintensity_PBS20=[normintensity_PBS20a; normintensity_PBS20b(:, 1:length(tme_PBS20))];
normintensity_PBS2=[normintensity_PBS2a(:, 1:length(tme_PBS2)); normintensity_PBS2b];

% gamma_LB=[gamma_LBc;gamma_LBd];
% gamma_PBS60=[gamma_PBS60c;gamma_PBS60d;];
% gamma_PBS20=[gamma_PBS20c;gamma_PBS20d];
% gamma_PBS2=[gamma_PBS2c;gamma_PBS20d];

%% Correct for photobleaching
parameter=28.9210;
[Cnew_LB, beta_LB, dCB_LB, dCT_LB, dCP_LB, CblExp_LB, unbFrac_LB, midx_LB]=photoCorrect(tme_LB, normintensity_LB, parameter);
[Cnew_PBS120, beta_PBS120, dCB_PBS120, dCT_PBS120, dCP_PBS120, CblExp_PBS120, unbFrac_PBS120, midx_PBS120]=photoCorrect(tme_PBS120, normintensity_PBS120, parameter);
[Cnew_PBS60, beta_PBS60, dCB_PBS60, dCT_PBS60, dCP_PBS60, CblExp_PBS60, unbFrac_PBS60, midx_PBS60]=photoCorrect(tme_PBS60, normintensity_PBS60, parameter);
[Cnew_PBS20, beta_PBS20, dCB_PBS20, dCT_PBS20, dCP_PBS20, CblExp_PBS20, unbFrac_PBS20, midx_PBS20]=photoCorrect(tme_PBS20, normintensity_PBS20, parameter);
[Cnew_PBS2, beta_PBS2, dCB_PBS2, dCT_PBS2, dCP_PBS2, CblExp_PBS2, unbFrac_PBS2, midx_PBS2]=photoCorrect(tme_PBS2, normintensity_PBS2, parameter);

[Cnew_LBMga, beta_LBMga, dCB_LBMga, dCT_LBMga, dCP_LBMga, CblExp_LBMga, unbFrac_LBMga, midx_LBMga]=photoCorrect(tme_LBMga, normintensity_LBMga, parameter);
[Cnew_LBEa, beta_LBEa, dCB_LBEa, dCT_LBEa, dCP_LBEa, CblExp_LBEa, unbFrac_LBEa, midx_LBEa]=photoCorrect(tme_LBEa, normintensity_LBEa, parameter);
[Cnew_LBtuna, beta_LBtuna, dCB_LBtuna, dCT_LBtuna, dCP_LBtuna, CblExp_LBtuna, unbFrac_LBtuna, midx_LBtuna]=photoCorrect(tme_LBtuna, normintensity_LBtuna, parameter);
[Cnew_LBdeta, beta_LBdeta, dCB_LBdeta, dCT_LBdeta, dCP_LBdeta, CblExp_LBdeta, unbFrac_LBdeta, midx_LBdeta]=photoCorrect(tme_LBdeta, normintensity_LBdeta, parameter);
[Cnew_LBvanco, beta_LBvanco, dCB_LBvanco, dCT_LBvanco, dCP_LBvanco, CblExp_LBvanco, unbFrac_LBvanco, midx_LBvanco]=photoCorrect(tme_LBvanco, normintensity_LBvanco, parameter);

[Cnew_LB5a, beta_LB5a, dCB_LB5a, dCT_LB5a, dCP_LB5a, CblExp_LB5a, unbFrac_LB5a, midx_LB5a]=photoCorrect(tme_LB5a, normintensity_LB5a, parameter);
[Cnew_LB20a, beta_LB20a, dCB_LB20a, dCT_LB20a, dCP_LB20a, CblExp_LB20a, unbFrac_LB20a, midx_LB20a]=photoCorrect(tme_LB20a, normintensity_LB20a, parameter);
%% fit corrected plot to exponential decay function to calculate tau
% tau0_LB=30;
% tau0_PBS120=1;
% tau0_PBS60=1;
% tau0_PBS20=10;
% tau0_PBS2=25;
% tau0_LBMga=30;
% tau0_LBEa=10;
% tau0_LBtuna=10;
% tau0_LBdeta=20;
% tau0_LBvanco=1;
% 
% [tau_LB, yhat_LB]=expFit(tme_LB, Cnew_LB, tau0_LB);
% [tau_LB5a, yhat_LB5a]=expFit(tme_LB5a, Cnew_LB5a);
% [tau_PBS120, yhat_PBS120]=expFit(tme_PBS120, Cnew_PBS120, tau0_PBS120);
% [tau_PBS60, yhat_PBS60]=expFit(tme_PBS60, Cnew_PBS60, tau0_PBS60);
% [tau_PBS20, yhat_PBS20]=expFit(tme_PBS20, Cnew_PBS20, tau0_PBS20);
% [tau_PBS2, yhat_PBS2]=expFit(tme_PBS2, Cnew_PBS2, tau0_PBS2);
% [tau_LBMga, yhat_LBMga]=expFit(tme_LBMga, Cnew_LBMga, tau0_LBMga);
% [tau_LBEa, yhat_LBEa]=expFit(tme_LBEa, Cnew_LBEa, tau0_LBEa);
% [tau_LBtuna, yhat_LBtuna]=expFit(tme_LBtuna, Cnew_LBtuna, tau0_LBtuna);
% [tau_LBdeta, yhat_LBdeta]=expFit(tme_LBdeta, Cnew_LBdeta, tau0_LBdeta);
% [tau_LBvanco, yhat_LBvanco]=expFit(tme_LBvanco, Cnew_LBvanco, tau0_LBvanco);

% %% Plot the corrected fluorescent traces
% cd(dirsave)
% alpha=0.3;
% 
% %plot the PBS incubation data
% figure, hold on
% ciplot(mean(Cnew_LB(:, 1:30), 1, 'omitnan')-std(Cnew_LB(:, 1:30), 0, 1, 'omitnan'), mean(Cnew_LB(:, 1:30), 1, 'omitnan')+std(Cnew_LB(:, 1:30), 0, 1, 'omitnan'), tme_LB(:, 1:30), colorcode2{1}, alpha)
% plot(tme_LB(:, 1:30), mean(Cnew_LB(:, 1:30), 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_PBS2(:, 1:30), 1, 'omitnan')-std(Cnew_PBS2(:, 1:30), 0, 1, 'omitnan'), mean(Cnew_PBS2(:, 1:30), 1, 'omitnan')+std(Cnew_PBS2(:, 1:30), 0, 1, 'omitnan'), tme_PBS2(:, 1:30), colorcode2{2}, alpha)
% plot(tme_PBS2(:, 1:30), mean(Cnew_PBS2(:, 1:30), 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_PBS20(:, 1:30), 1, 'omitnan')-std(Cnew_PBS20(:, 1:30), 0, 1, 'omitnan'), mean(Cnew_PBS20(:, 1:30), 1, 'omitnan')+std(Cnew_PBS20(:, 1:30), 0, 1, 'omitnan'), tme_PBS20(:, 1:30), colorcode2{3}, alpha)
% plot(tme_PBS20(:, 1:30), mean(Cnew_PBS20(:, 1:30), 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_PBS60(:, 1:30), 1, 'omitnan')-std(Cnew_PBS60(:, 1:30), 0, 1, 'omitnan'), mean(Cnew_PBS60(:, 1:30), 1, 'omitnan')+std(Cnew_PBS60(:, 1:30), 0, 1, 'omitnan'), tme_PBS60(:, 1:30), colorcode2{4}, alpha)
% plot(tme_PBS60(:, 1:30), mean(Cnew_PBS60(:, 1:30), 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_PBS120(:, 1:20), 1, 'omitnan')-std(Cnew_PBS120(:, 1:20), 0, 1, 'omitnan'), mean(Cnew_PBS120(:, 1:20), 1, 'omitnan')+std(Cnew_PBS120(:, 1:20), 0, 1, 'omitnan'), tme_PBS120(:, 1:20), colorcode2{5}, alpha)
% plot(tme_PBS120(:, 1:20), mean(Cnew_PBS120(:, 1:20), 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% 
% legend({'LB, std', 'LB, mean', 'PBS 2 min, std', 'PBS 2 min, mean', 'PBS 20 min, std', 'PBS 20 min, mean', 'PBS 1 hour, std', 'PBS 1 hour, mean', 'PBS 2 hour, std', 'PBS 2 hour, mean'}, 'Location', 'southeast')
% xlabel('Time (minutes)')
% ylabel('Normalized Fluorescence (A.U.)')
% ylim([0 1.05])
% xlim([0 45])
% saveas(gcf, 'pbsCorrectedIntensity.png')
% saveas(gcf, 'pbsCorrectedIntensity.fig')
% 
% %plot the LB + treatment data
% figure, hold on
% ciplot(mean(Cnew_LB(:, 1:45), 1, 'omitnan')-std(Cnew_LB(:, 1:45), 0, 1, 'omitnan'), mean(Cnew_LB(:, 1:45), 1, 'omitnan')+std(Cnew_LB(:, 1:45), 0, 1, 'omitnan'), tme_LB(:, 1:45), colorcode2{1}, alpha)
% plot(tme_LB(:, 1:45), mean(Cnew_LB(:, 1:45), 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_LBMga(:, 1:45), 1, 'omitnan')-std(Cnew_LBMga(:, 1:45), 0, 1, 'omitnan'), mean(Cnew_LBMga(:, 1:45), 1, 'omitnan')+std(Cnew_LBMga(:, 1:45), 0, 1, 'omitnan'), tme_LBMga(:, 1:45), colorcode2{6}, alpha)
% plot(tme_LBMga(:, 1:45), mean(Cnew_LBMga(:, 1:45), 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_LBEa(:, 1:45), 1, 'omitnan')-std(Cnew_LBEa(:, 1:45), 0, 1, 'omitnan'), mean(Cnew_LBEa(:, 1:45), 1, 'omitnan')+std(Cnew_LBEa(:, 1:45), 0, 1, 'omitnan'), tme_LBEa(:, 1:45), colorcode2{7}, alpha)
% plot(tme_LBEa(:, 1:45), mean(Cnew_LBEa(:, 1:45), 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_LBtuna(:, 1:45), 1, 'omitnan')-std(Cnew_LBtuna(:, 1:45), 0, 1, 'omitnan'), mean(Cnew_LBtuna(:, 1:45), 1, 'omitnan')+std(Cnew_LBtuna(:, 1:45), 0, 1, 'omitnan'), tme_LBtuna(:, 1:45), colorcode2{8}, alpha)
% plot(tme_LBtuna(:, 1:45), mean(Cnew_LBtuna(:, 1:45), 1, 'omitnan'), 'Color', colorcode{8}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_LBdeta(:, 1:45), 1, 'omitnan')-std(Cnew_LBdeta(:, 1:45), 0, 1, 'omitnan'), mean(Cnew_LBdeta(:, 1:45), 1, 'omitnan')+std(Cnew_LBdeta(:, 1:45), 0, 1, 'omitnan'), tme_LBdeta(:, 1:45), colorcode2{9}, alpha)
% plot(tme_LBdeta(:, 1:45), mean(Cnew_LBdeta(:, 1:45), 1, 'omitnan'), 'Color', colorcode{9}, 'LineWidth', 1)
% 
% ciplot(mean(Cnew_LBvanco(:, 1:45), 1, 'omitnan')-std(Cnew_LBvanco(:, 1:45), 0, 1, 'omitnan'), mean(Cnew_LBvanco(:, 1:45), 1, 'omitnan')+std(Cnew_LBvanco(:, 1:45), 0, 1, 'omitnan'), tme_LBvanco(:, 1:45), colorcode2{10}, alpha)
% plot(tme_LBvanco(:, 1:45), mean(Cnew_LBvanco(:, 1:45), 1, 'omitnan'), 'Color', colorcode{10}, 'LineWidth', 1)
% 
% legend({'LB, std', 'LB, mean', 'LB Mg^{2+}, std', 'LB Mg^{2+}, mean', 'LB EDTA, std', 'LB EDTA, mean', 'LB tunicamycin, std', 'LB tunicamycin, mean', 'LB 5% NLS, std', 'LB 5% NLS, mean', 'LB vancomycin, std', 'LB vancomycin, mean'}, 'Location', 'southeast')
% xlabel('Time (minutes)')
% ylabel('Normalized Fluorescence (A.U.)')
% ylim([0 1.05])
% xlim([0 60])
% saveas(gcf, 'LBcorrectedIntensity.png')
% saveas(gcf, 'LBcorrectedIntensity.fig')
% 
% %% plot beta values vs PBS incubation
% linearCoef1 = polyfit([zeros(1, length(beta_LB)), repelem(2, length(beta_PBS2)), repelem(20, length(beta_PBS20)), repelem(60, length(beta_PBS60)), repelem(120, length(beta_PBS120))],[beta_LB', beta_PBS2', beta_PBS20', beta_PBS60', beta_PBS120'],1);
% linearFit1= polyval(linearCoef1,[0 2 20 60 120]);
% 
% linearCoef1 = polyfit([0 2 20 60 120], [mean(beta_LB, 'omitnan'), mean(beta_PBS2, 'omitnan'), mean(beta_PBS20, 'omitnan'), mean(beta_PBS60, 'omitnan'), mean(beta_PBS120, 'omitnan')], 1);
% linearFit1 = polyval(linearCoef1, [0 2 20 60 120]);
% 
% %plot beta vs PBS incubation
% figure, hold on
% scatter(zeros(1, length(beta_LB)), beta_LB', 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter(repelem(120, length(beta_PBS120)), beta_PBS120', 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter(repelem(60, length(beta_PBS60)), beta_PBS60', 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter(repelem(20, length(beta_PBS20)), beta_PBS20', 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter(repelem(2, length(beta_PBS2)), beta_PBS2', 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter([0 2 20 60 120], [mean(beta_LB, 'omitnan'), mean(beta_PBS2, 'omitnan'), mean(beta_PBS20, 'omitnan'), mean(beta_PBS60, 'omitnan'), mean(beta_PBS120, 'omitnan')], 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
% 
% errorbar(0, mean(beta_LB, 'omitnan'), std(beta_LB, 'omitnan'), 'Color', 'black')
% errorbar(2, mean(beta_PBS2, 'omitnan'), std(beta_PBS2, 'omitnan'), 'Color', 'black')
% errorbar(20, mean(beta_PBS20, 'omitnan'), std(beta_PBS20, 'omitnan'), 'Color', 'black')
% errorbar(60, mean(beta_PBS60, 'omitnan'), std(beta_PBS60, 'omitnan'), 'Color', 'black')
% errorbar(120, mean(beta_PBS120, 'omitnan'), std(beta_PBS120, 'omitnan'), 'Color', 'black')
% 
% plot([0 2 20 60 120], linearFit1, '--b')
% 
% xlim([-2 122])
% ylim([0 1.1])
% xlabel('Time in PBS (minutes)')
% ylabel('Beta (A.U.)')
% 
% saveas(gcf, 'beta_vs_PBS.png')
% saveas(gcf, 'beta_vs_PBS.fig')
% 
% %% plot tau as a function of PBS incubation
% linearCoef2 = polyfit([0 2 20 60 120], [mean(tau_LB, 'omitnan'), mean(tau_PBS2, 'omitnan'), mean(tau_PBS20, 'omitnan'), mean(tau_PBS60, 'omitnan'), mean(tau_PBS120, 'omitnan')], 1);
% linearFit2 = polyval(linearCoef2, [0 2 20 60 120]);
% 
% figure, hold on
% scatter(zeros(1, length(tau_LB)), tau_LB, 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter(repelem(120, length(tau_PBS120)), tau_PBS120, 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter(repelem(60, length(tau_PBS60)), tau_PBS60, 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter(repelem(20, length(tau_PBS20)), tau_PBS20, 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter(repelem(2, length(tau_PBS2)), tau_PBS2, 'MarkerFaceColor', colorcode{14}, 'MarkerEdgeColor', colorcode{14})
% scatter([0 2 20 60 120], [mean(tau_LB, 'omitnan'), mean(tau_PBS2, 'omitnan'), mean(tau_PBS20, 'omitnan'), mean(tau_PBS60, 'omitnan'), mean(tau_PBS120, 'omitnan')], 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
% 
% errorbar(0, mean(tau_LB, 'omitnan'), std(tau_LB, 'omitnan'), 'Color', 'black')
% errorbar(2, mean(tau_PBS2, 'omitnan'), std(tau_PBS2, 'omitnan'), 'Color', 'black')
% errorbar(20, mean(tau_PBS20, 'omitnan'), std(tau_PBS20, 'omitnan'), 'Color', 'black')
% errorbar(60, mean(tau_PBS60, 'omitnan'), std(tau_PBS60, 'omitnan'), 'Color', 'black')
% errorbar(120, mean(tau_PBS120, 'omitnan'), std(tau_PBS120, 'omitnan'), 'Color', 'black')
% 
% plot([0 2 20 60 120], linearFit2, '--b')
% 
% xlim([-2 122])
% %ylim([0 1.1])
% xlabel('Time in PBS (minutes)')
% ylabel('Tau (minutes)')
% 
% [DF_LB]=tauProb(tau_LB);
% [DF_PBS2]=tauProb(tau_PBS2);
% [DF_PBS20]=tauProb(tau_PBS20);
% [DF_PBS60]=tauProb(tau_PBS60);
% [DF_PBS120]=tauProb(tau_PBS120);
% 
% figure, hold on
% fnplt(DF_LB, colorcode{1}, 2)
% fnplt(DF_PBS2, colorcode{2}, 2)
% fnplt(DF_PBS20, colorcode{3}, 2)
% fnplt(DF_PBS60, colorcode{4}, 2)
% fnplt(DF_PBS120, colorcode{5}, 2)
% 
% % figure, hold on
% % histogram(tau_LB, length(tau_LB), 'Normalization', 'probability')
% % histogram(tau_PBS2, length(tau_PBS2), 'Normalization', 'probability')
% % histogram(tau_PBS20, length(tau_PBS20), 'Normalization', 'probability')
% % histogram(tau_PBS60, length(tau_PBS60), 'Normalization', 'probability')
% % histogram(tau_PBS120, length(tau_PBS120), 'Normalization', 'probability')
% 
% saveas(gcf, 'tau_vs_PBS.png')
% saveas(gcf, 'tau_vs_PBS.fig')

% %% plot tau as a function of beta
% figure, hold on
% scatter(beta_LB, tau_LB, 'MarkerFaceColor', colorcode{1}, 'MarkerEdgeColor', colorcode{1})
% scatter(beta_PBS2, tau_PBS2, 'MarkerFaceColor', colorcode{2}, 'MarkerEdgeColor', colorcode{2})
% scatter(beta_PBS20, tau_PBS20, 'MarkerFaceColor', colorcode{3}, 'MarkerEdgeColor', colorcode{3})
% scatter(beta_PBS60, tau_PBS60, 'MarkerFaceColor', colorcode{4}, 'MarkerEdgeColor', colorcode{4})
% scatter(beta_PBS120, tau_PBS120, 'MarkerFaceColor', colorcode{5}, 'MarkerEdgeColor', colorcode{5})
% scatter(beta_LBMga, tau_LBMga, 'MarkerFaceColor', colorcode{6}, 'MarkerEdgeColor', colorcode{6})
% scatter(beta_LBEa, tau_LBEa, 'MarkerFaceColor', colorcode{7}, 'MarkerEdgeColor', colorcode{7})
% scatter(beta_LBtuna, tau_LBtuna, 'MarkerFaceColor', colorcode{8}, 'MarkerEdgeColor', colorcode{8})
% scatter(beta_LBdeta, tau_LBdeta, 'MarkerFaceColor', colorcode{9}, 'MarkerEdgeColor', colorcode{9})
% scatter(beta_LBvanco, tau_LBvanco, 'MarkerFaceColor', colorcode{10}, 'MarkerEdgeColor', colorcode{10})
% xlabel('Beta (A.U.)')
% ylabel('Tau (minutes)')
% xlim([0 1])
% legend({'LB', 'PBS 2 min', 'PBS 20 min', 'PBS 1 hour', 'PBS 2 hour', 'LB + Mg^{2+}', 'LB + EDTA', 'LB + tunicamycin', 'LB + 5% NLS', 'LB + vancomycin'}, 'Location', 'northwest')
% saveas(gcf, 'tau_vs_beta.png')
% saveas(gcf, 'tau_vs_beta.fig')
% 
% %% plot beta as a function of gamma
% figure, hold on
% scatter(gamma_LB(midx_LB,:), beta_LB, 'MarkerFaceColor', colorcode{1}, 'MarkerEdgeColor', colorcode{1})
% scatter(gamma_PBS2(midx_PBS2,:), beta_PBS2, 'MarkerFaceColor', colorcode{2}, 'MarkerEdgeColor', colorcode{2})
% scatter(gamma_PBS20(midx_PBS20,:), beta_PBS20, 'MarkerFaceColor', colorcode{3}, 'MarkerEdgeColor', colorcode{3})
% scatter(gamma_PBS60(midx_PBS60,:), beta_PBS60, 'MarkerFaceColor', colorcode{4}, 'MarkerEdgeColor', colorcode{4})
% scatter(gamma_PBS120(midx_PBS120,:), beta_PBS120, 'MarkerFaceColor', colorcode{5}, 'MarkerEdgeColor', colorcode{5})
% scatter(gamma_LBMgc(midx_LBMga,:), beta_LBMga, 'MarkerFaceColor', colorcode{6}, 'MarkerEdgeColor', colorcode{6})
% scatter(gamma_LBEc(midx_LBEa,:), beta_LBEa, 'MarkerFaceColor', colorcode{7}, 'MarkerEdgeColor', colorcode{7})
% scatter(gamma_LBtunc(midx_LBtuna,:), beta_LBtuna, 'MarkerFaceColor', colorcode{8}, 'MarkerEdgeColor', colorcode{8})
% scatter(gamma_LBdetc(midx_LBdeta,:), beta_LBdeta, 'MarkerFaceColor', colorcode{9}, 'MarkerEdgeColor', colorcode{9})
% scatter(gamma_LBvanc(midx_LBvanco,:), beta_LBvanco, 'MarkerFaceColor', colorcode{10}, 'MarkerEdgeColor', colorcode{10})
% xlabel('Gamma (A.U.)')
% ylabel('Beta (A.U.)')
% %xlim([0 1])
% legend({'LB', 'PBS 2 min', 'PBS 20 min', 'PBS 1 hour', 'PBS 2 hour', 'LB + Mg^{2+}', 'LB + EDTA', 'LB + tunicamycin', 'LB + 5% NLS', 'LB + vancomycin'}, 'Location', 'northwest')
% saveas(gcf, 'beta_vs_gamma.png')
% saveas(gcf, 'beta_vs_gamma.fig')
% 
% %% plot tau as a function of gamma
% figure, hold on
% scatter(gamma_LB(midx_LB,:), tau_LB, 'MarkerFaceColor', colorcode{1}, 'MarkerEdgeColor', colorcode{1})
% scatter(gamma_PBS2(midx_PBS2,:), tau_PBS2, 'MarkerFaceColor', colorcode{2}, 'MarkerEdgeColor', colorcode{2})
% scatter(gamma_PBS20(midx_PBS20,:), tau_PBS20, 'MarkerFaceColor', colorcode{3}, 'MarkerEdgeColor', colorcode{3})
% scatter(gamma_PBS60(midx_PBS60,:), tau_PBS60, 'MarkerFaceColor', colorcode{4}, 'MarkerEdgeColor', colorcode{4})
% scatter(gamma_PBS120(midx_PBS120,:), tau_PBS120, 'MarkerFaceColor', colorcode{5}, 'MarkerEdgeColor', colorcode{5})
% scatter(gamma_LBMgc(midx_LBMga,:), tau_LBMga, 'MarkerFaceColor', colorcode{6}, 'MarkerEdgeColor', colorcode{6})
% scatter(gamma_LBEc(midx_LBEa,:), tau_LBEa, 'MarkerFaceColor', colorcode{7}, 'MarkerEdgeColor', colorcode{7})
% scatter(gamma_LBtunc(midx_LBtuna,:), tau_LBtuna, 'MarkerFaceColor', colorcode{8}, 'MarkerEdgeColor', colorcode{8})
% scatter(gamma_LBdetc(midx_LBdeta,:), tau_LBdeta, 'MarkerFaceColor', colorcode{9}, 'MarkerEdgeColor', colorcode{9})
% scatter(gamma_LBvanc(midx_LBvanco,:), tau_LBvanco, 'MarkerFaceColor', colorcode{10}, 'MarkerEdgeColor', colorcode{10})
% xlabel('Gamma (A.U.)')
% ylabel('Tau (minutes)')
% %xlim([0 1])
% legend({'LB', 'PBS 2 min', 'PBS 20 min', 'PBS 1 hour', 'PBS 2 hour', 'LB + Mg^{2+}', 'LB + EDTA', 'LB + tunicamycin', 'LB + 5% NLS', 'LB + vancomycin'}, 'Location', 'northwest')
% saveas(gcf, 'tau_vs_gamma.png')
% saveas(gcf, 'tau_vs_gamma.fig')

%% Functions
function [normintensity, intensity, time, tme, imend]=dataNormalize(datadir, imstart)
        
        %pre-allocate variables
        intensity=[];
      
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'time')
            
            for n=1:height(icell_intensity)

                %make sure there are fluor readings during the initial frame, otherwise the adjust. and norm. will look off
                if ~isnan(icell_intensity(n, imstart))
                    intensity=[intensity; icell_intensity(n, :)];
                end
            end
            
            if i==1
                tme=time; %pre-set new time vector
            end

        end
        
        %set the limit of detection 
        lmt=1300;
        adjintensity=intensity;
        adjintensity(adjintensity<=lmt)=NaN; 
        idx=sum(isnan(adjintensity));
        imend=find(idx>=(height(adjintensity)*0.5)); %the time point where half the values are nan
        imend=imend(min(find(imend>imstart))); 
            
        %adjust the time vector
        tme=tme(imstart:end)-tme(imstart);
                                                                      
        %interpolate the fluor values during detergent perfusion
        if imstart<9
            omit=[2:5]+imstart;
            idx=setdiff(1:length(time), omit);
            v=adjintensity(:, idx);
            vq=nan(height(v), length(omit));
            for n=1:height(v)
                vq(n, :)=interp1(idx, v(n, :), omit);
            end
            adjintensity(:, omit)=vq;
        end

        %adjust the background
        adjintensity = adjintensity(:, 1:imend);
        adjintensity = adjintensity-lmt;
        
        %normalize to the initial pre-lysis frame
        normintensity=adjintensity(:, imstart:end)./adjintensity(:,imstart);
        normintensity(normintensity<0)=0;
        
end

function [initial, final, gamma]=phaseRatio(datadir, imstart, imend)
        
        %pre-allocate variables
        intensity=[];
      
        %go through the data for each position
        for i=1:length(datadir)
            
            %load decayMeasure .mat file
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity')

            for n=1:height(icell_intensity)
                
                intensity=[intensity; icell_intensity(n, :)];

            end
            %the tricky thing about this is that I don't know if the ratio
            %belongs to a corresponding cell in Cnew because some may have
            %gotten filtered out
        end
        
        initial=intensity(:, imstart);
        final=intensity(:, imend);
        gamma=final./initial;
        
end

function [Cnew, beta, dCB, dCT, dCP, Cbl_exp, unb_frac, midx]=photoCorrect(tme, normintensity, alpha)
        
        %Correct for photobleaching
        %calculate dt (the dt between frames may vary in a single run)
        dt=diff(tme);
        
        %this formula comes from the slope and intercept calculated for the 1.2, 2,
        %and 3 second tau vs frame rate controls 
       
        parameter = (alpha.*dt) + 0.2482; 
        
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
                
                dCB(n,i) = normintensity(n,i)/parameter(i); %this is the amount of photobleaching that occured in our measured value

                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)*unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
        end      
        
        beta=1;
        midx=1;
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

function [tau, yhat]=expFit(tme, Cnew, tau0)
        
        %pre-allocate variables
        tau=[];
        yhat=[];
        
        if nargin < 3
            for i=1:height(Cnew)
                linearCoef = polyfit(tme, Cnew(i,:), 1);
                linearFit = polyval(linearCoef, tme);  

                tau=[tau, linearCoef(1)];
                yhat=[yhat; linearFit]; 
            end
        else         
            for i=1:height(Cnew)

                modelfun=@(tau,t)exp(-t./tau);
                tau_temp=nlinfit(tme, Cnew(i,:), modelfun, tau0);
                y_hat=modelfun(tau_temp, tme);   

                tau=[tau, tau_temp];
                yhat=[yhat; y_hat]; 

            end           
        end
end

function [DF] = tauProb(tau)

    h = histogram(tau, length(tau));
    heights=h.BinCounts/h.NumBins;
    centers=tau;
    n = length(centers);
    w = centers(2)-centers(1);
    t = linspace(centers(1) - w/2, centers(end) + w/2, n + 1);
    dt = diff(t);
    Fvals = cumsum([0, heights.*dt]);
    F = spline(t, [0, Fvals, 0]);
    DF = fnder(F);
    
end
