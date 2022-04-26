%sanity_checks.m
%Author: Zarina Akbary
%Date: 04/19/2022
%Purpose: this is basically a sandbox for me to look through my data and
%keep track of how I do it. 

clear, close all

okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

%% User Input
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/aggregate/';
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/figures/';

control1_basenames = {'04112022_Exp1'};
control1_labels = {'Frame Rate = 1 min'};

control2_basenames = {'12082021_Exp2'};
control2_labels = {'Frame Rate = 3 s'};

%% quick check
cd(dirsave)

[control1_cellTrace, control1_bgTrace, control1_adjTrace, control1_times, control1_lcell, control1_frameRate] = loadData(dirpath, control1_basenames);
[control1_lysis, control1_elongation] = lysisCalc(control1_times, control1_lcell);
[control1_normTrace] = dataNormalize(control1_times, control1_adjTrace, control1_lysis);
[control1_correctedTrace]=photoCorrect(control1_normTrace(:, 1), control1_normTrace(:, 2), 0.004094, 0.028591, 4);

[control2_cellTrace, control2_bgTrace, control2_adjTrace, control2_times, control2_lcell, control2_frameRate] = loadData(dirpath, control2_basenames);
[control2_lysis, control2_elongation] = lysisCalc(control2_times, control2_lcell);
[control2_normTrace] = dataNormalize(control2_times, control2_adjTrace, control2_lysis);
[control2_correctedTrace]=photoCorrect(control2_normTrace(:, 1), control2_normTrace(:, 2), 0.004094, 0.028591, 1);

alpha0 = 1;
rho0 = 0;

[control1_alpha, control1_alpha_yhat] = alphaCalc(control1_normTrace(:, 1), control1_normTrace(:, 2), alpha0);
[control1_rho, control1_rho_yhat] = rhoCalc(control1_normTrace(:, 1), control1_normTrace(:, 2), rho0, 4);

%plot the fit of alpha
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1 %:2
    
%     subplot(1, 2, i)
    
    x = control1_normTrace{i,1};
    
    y1 = control1_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control1_alpha_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([0 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control1_labels{i})
end
% saveas(gcf, 'figure02.png')
% saveas(gcf, 'figure02.fig')

%plot the fit of rho
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1 %:2
    
    %subplot(1, 2, i)
    
    [~, ncol] = size(control1_normTrace{i, 1});
    x = [1:ncol]*4;
    
    y1 = control1_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control1_rho_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([0 1.1])
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence')
    title(control1_labels{i})
end
% saveas(gcf, 'figure01.png')
% saveas(gcf, 'figure01.fig')

%corrected trace
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1 %:2
    
    %subplot(1, 2, i)
    
    x = control1_normTrace{i,1};
    x2 = control2_normTrace{i,1};
   
    y1 = control1_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control1_correctedTrace{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    y3 = control2_correctedTrace{i};
    ybar3 = mean(y3, 1, 'omitnan');
    ybar3 = movingaverage(ybar3, 3);
    err3 = std(y3, 0, 1, 'omitnan');
    
    errorbar(x, ybar2, err2, 'Color', okabeIto{i}), hold on
    %errorbar(x, ybar2, err2, '--k')
    errorbar(x2, ybar3, err3, 'Color', okabeIto{i+1})
    %errorbar(x, ybar2, err2, '--k')
    ylim([0 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    %title(control1_labels{i})
end
%% load data
[exp_cellTrace, exp_bgTrace, exp_adjTrace, exp_times, exp_lcell, exp_frameRate] = loadData(dirpath, exp_basenames);
[control2_cellTrace, control2_bgTrace, control2_adjTrace, control2_times, control2_lcell, control2_frameRate] = loadData(dirpath, control2_basenames);
[control3_cellTrace, control3_bgTrace, control3_adjTrace, control3_times, control3_lcell, control3_frameRate] = loadData(dirpath, control3_basenames);
[pbs_cellTrace, pbs_bgTrace, pbs_adjTrace, pbs_times, pbs_lcell, pbs_frameRate] = loadData(dirpath, pbs_basenames);

%% calculate lysis frames and elongation rate
[exp_lysis, exp_elongation] = lysisCalc(exp_times, exp_lcell);
[control2_lysis, control2_elongation] = lysisCalc(control2_times, control2_lcell);
[control3_lysis, control3_elongation] = lysisCalc(control3_times, control3_lcell);
[pbs_lysis, pbs_elongation] = lysisCalc(pbs_times, pbs_lcell);

%% normalize traces
[exp_normTrace] = dataNormalize(exp_times, exp_adjTrace, exp_lysis);
[control2_normTrace] = dataNormalize(control2_times, control2_adjTrace, control2_lysis);
[control3_normTrace] = dataNormalize(control3_times, control3_adjTrace, control3_lysis);
[pbs_normTrace] = dataNormalize(pbs_times, pbs_adjTrace, pbs_lysis);

%% calculate alpha and rho
alpha0 = 1;
rho0 = 0;

[exp_alpha, exp_alpha_yhat] = alphaCalc(exp_normTrace(:, 1), exp_normTrace(:, 2), alpha0);
[control2_alpha, control2_alpha_yhat] = alphaCalc(control2_normTrace(:, 1), control2_normTrace(:, 2), alpha0);
[control3_alpha, control3_alpha_yhat] = alphaCalc(control3_normTrace(:, 1), control3_normTrace(:, 2), alpha0);
[pbs_alpha, pbs_alpha_yhat] = alphaCalc(pbs_normTrace(:, 1), pbs_normTrace(:, 2), alpha0);

[exp_rho, exp_rho_yhat] = rhoCalc(exp_normTrace(:, 1), exp_normTrace(:, 2), rho0);
[control2_rho, control2_rho_yhat] = rhoCalc(control2_normTrace(:, 1), control2_normTrace(:, 2), rho0);
[control3_rho, control3_rho_yhat] = rhoCalc(control3_normTrace(:, 1), control3_normTrace(:, 2), rho0);
[pbs_rho, pbs_rho_yhat] = rhoCalc(pbs_normTrace(:, 1), pbs_normTrace(:, 2), rho0);

%% fit rho and alpha to linear equations (without the 1 minute data)
[control2_alphaCoef, control2_alphaFit] = alphaLinear(control2_frameRate(1:5, 1), control2_alpha(1:5, 1));
[control3_alphaCoef, control3_alphaFit] = alphaLinear(control3_frameRate(1:2, 1), control3_alpha(1:2, 1));

[control2_rhoCoef, control2_rhoFit] = rhoLinear(control2_frameRate(1:5, 1), control2_rho(1:5, 1));
[control3_rhoCoef, control3_rhoFit] = rhoLinear(control3_frameRate(1:2, 1), control3_rho(1:2, 1));

%% correct for photobleaching
[exp_correctedTrace]=photoCorrect(exp_normTrace(:, 1), exp_normTrace(:, 2), control2_alphaCoef(1), control2_alphaCoef(2));
[control2_correctedTrace]=photoCorrect(control2_normTrace(:, 1), control2_normTrace(:, 2), control2_alphaCoef(1), control2_alphaCoef(2));
%[control3_correctedTrace]=photoCorrect(control3_normTrace(:, 1), control3_normTrace(:, 2), control3_alphaCoef(1), control3_alphaCoef(2));
[control3_correctedTrace]=photoCorrect(control3_normTrace(:, 1), control3_normTrace(:, 2), 15, control3_alphaCoef(2));
[pbs_correctedTrace]=photoCorrect(pbs_normTrace(:, 1), pbs_normTrace(:, 2), 15, control3_alphaCoef(2));

%% calculate C effective
[exp_ceff]=effectiveLoss(exp_normTrace(:, 1), exp_correctedTrace);
[control2_ceff]=effectiveLoss(control2_normTrace(:, 1), control2_correctedTrace);
[control3_ceff]=effectiveLoss(control3_normTrace(:, 1), control3_correctedTrace);
[pbs_ceff]=effectiveLoss(pbs_normTrace(:, 1), pbs_correctedTrace);

%% calculate tau
tau0 = 50;

[exp_tau, exp_tau_yhat] = tauCalc(exp_normTrace(:, 1), exp_correctedTrace, tau0);
[control2_tau, control2_tau_yhat] = tauCalc(control2_normTrace(:, 1), control2_correctedTrace, tau0);
[control3_tau, control3_tau_yhat] = tauCalc(control3_normTrace(:, 1), control3_correctedTrace, tau0);
[pbs_tau, pbs_tau_yhat] = tauCalc(pbs_normTrace(:, 1), pbs_correctedTrace, tau0);

%% the story starts with the controls, so plot the raw fluor. traces of those first
cd(dirsave)

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:8
    subplot(4, 2, i)
    x = control2_times{i};
    y1 = control2_cellTrace{i};
    y2 = control2_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(control2_labels{i})
end
saveas(gcf, 'figure02.png')
saveas(gcf, 'figure02.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:2
    subplot(1, 2, i)
    x = control3_times{i};
    y1 = control3_cellTrace{i};
    y2 = control3_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(control3_labels{i})
end
saveas(gcf, 'figure03.png')
saveas(gcf, 'figure03.fig')

%% after this, show the fits to alpha and rho, since that'll also show the normalization
cd(dirsave)

%plot the fit of alpha
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:8
    subplot(4, 2, i)
    x = control2_normTrace{i,1};
    
    y1 = control2_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control2_alpha_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([0 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure04.png')
saveas(gcf, 'figure04.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:2
    
    subplot(1, 2, i)
    
    x = control3_normTrace{i,1};
    
    y1 = control3_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control3_alpha_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([0 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end
saveas(gcf, 'figure05.png')
saveas(gcf, 'figure05.fig')

%plot the fit of rho
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:8
    subplot(4, 2, i)
    [~, ncol] = size(control2_normTrace{i, 1});
    x = 1:ncol;
    
    y1 = control2_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control2_rho_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([0 1.1])
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure06.png')
saveas(gcf, 'figure06.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:2
    
    subplot(1, 2, i)
    
    [~, ncol] = size(control3_normTrace{i, 1});
    x = 1:ncol;
    
    y1 = control3_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control3_rho_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([0 1.1])
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end
saveas(gcf, 'figure07.png')
saveas(gcf, 'figure07.fig')

%% plot normalized traces over frame number

p = 0;
%plot the fit of rho
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=[2, 3, 7]
    
    p = p + 1;
    
    subplot(1, 3, p)
    [~, ncol] = size(control2_normTrace{i, 1});
    x = 1:ncol;
    
    y1 = control2_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
%     y2 = control2_rho_yhat{i};
%     ybar2 = mean(y2, 1, 'omitnan');
%     ybar2 = movingaverage(ybar2, 3);
%     err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    %errorbar(x, ybar2, err2, '--k')
    ylim([0 1.1])
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end

%% plot the corrected trace
p = 0;
%plot the fit of rho
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=[2, 3, 7]
    
    p = p + 1;
    
    subplot(1, 3, p)
%     [~, ncol] = size(control2_normTrace{i, 1});
%     x = 1:ncol;
    
    x = control2_normTrace{i, 1};
    y1 = control2_correctedTrace{i, 1};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
%     y2 = control2_rho_yhat{i};
%     ybar2 = mean(y2, 1, 'omitnan');
%     ybar2 = movingaverage(ybar2, 3);
%     err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    %errorbar(x, ybar2, err2, '--k')
    ylim([0 1.3])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
%% plot alpha vs frame rate & rho vs frame rate for the controls
dx1 = [];
dx2 = [];

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:5   
    x = control2_frameRate{i}(end);
    dx1 = [dx1, x];
 
    y = control2_alpha{i}';
    ybar = mean(y, 2, 'omitnan');
    err = std(y, 0, 2, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')    
end    

for i=1:2
    x = control3_frameRate{i}(end);
    dx2 = [dx2, x];
 
    y = control3_alpha{i}';
    ybar = mean(y, 2, 'omitnan');
    err = std(y, 0, 2, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')  
end
plot([0, dx1], control2_alphaFit, '--k', 'LineWidth', 1.5)
plot([0, dx2], control3_alphaFit, '--k', 'LineWidth', 1.5)
xlabel('Frame Rate (minutes)')
ylabel('\alpha (min)')
xlim([-0.002 x+x/10])
%legend([p1(1), p1(2), p2(1), p3(1)], [control2_labels{1}, control2_labels{2}, 'mean', 'linear fit'])
caption1 = sprintf('%f * frame rate + %f', control2_alphaCoef(1), control2_alphaCoef(2));
caption2 = sprintf('%f * frame rate + %f', control3_alphaCoef(1), control3_alphaCoef(2));
text(0.35, 10, ['\alpha = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
text(0.3, 2, ['\alpha = ' caption2], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
saveas(gcf, 'figure08.png')
saveas(gcf, 'figure08.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:5    
    x = control2_frameRate{i}(end);
    y = control2_rho{i}';
    ybar = mean(y, 2, 'omitnan');
    err = std(y, 0, 2, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')
end   

for i=1:2
    x = control3_frameRate{i}(end);
    y = control3_rho{i}';
    ybar = mean(y, 2, 'omitnan');
    err = std(y, 0, 2, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')  
end
plot([0, dx1], control2_rhoFit, '--k', 'LineWidth', 1.5)
plot([0, dx2], control3_rhoFit, '--k', 'LineWidth', 1.5)
xlabel('Frame Rate (minutes)')
ylabel('\rho')
ylim([0 .11])
xlim([-0.002 0.55])
%legend([p1(1), p1(2), p2(1), p3(1)], [control2_labels{1}, control2_labels{2}, 'mean', 'linear fit'])
caption1 = sprintf('%f * frame rate + %f', control2_rhoCoef(1), control2_rhoCoef(2));
caption2 = sprintf('%f * frame rate + %f', control3_rhoCoef(1), control3_rhoCoef(2));
text(0.35, 0.03, ['\rho = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
text(0.3, 0.06, ['\rho = ' caption2], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
saveas(gcf, 'figure09.png')
saveas(gcf, 'figure09.fig')

%% now, plot the correction for the GFP-only controls
cd(dirsave)

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(control2_labels)
    
    subplot(3, 4, i)
    x = control2_normTrace{i, 1};
    y = control2_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
%% now, plot the correction for the controls
cd(dirsave)

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:2
    
    p = p+1;
    subplot(2, 2, p)
    x = control2_times{i};
    y1 = control2_cellTrace{i};
    y2 = control2_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(control2_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(control2_labels{i})
    
    p = p+1;
    subplot(2, 2, p)
    x = control2_normTrace{i, 1};
    y = control2_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure10.png')
saveas(gcf, 'figure10.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=3:5
    
    p = p+1;
    subplot(3, 2, p)
    x = control2_times{i};
    y1 = control2_cellTrace{i};
    y2 = control2_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(control2_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(control2_labels{i})
    
    p = p+1;
    subplot(3, 2, p)
    x = control2_normTrace{i, 1};
    y = control2_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure11.png')
saveas(gcf, 'figure11.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=6:7
    
    p = p+1;
    subplot(2, 2, p)
    x = control2_times{i};
    y1 = control2_cellTrace{i};
    y2 = control2_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(control2_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(control2_labels{i})
    
    p = p+1;
    subplot(2, 2, p)
    x = control2_normTrace{i, 1};
    y = control2_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure12.png')
saveas(gcf, 'figure12.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=8:9
    
    p = p+1;
    subplot(2, 2, p)
    x = control2_times{i};
    y1 = control2_cellTrace{i};
    y2 = control2_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(control2_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(control2_labels{i})
    
    p = p+1;
    subplot(2, 2, p)
    x = control2_normTrace{i, 1};
    y = control2_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure13.png')
saveas(gcf, 'figure13.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=10:12
    
    p = p+1;
    subplot(3, 2, p)
    x = control2_times{i};
    y1 = control2_cellTrace{i};
    y2 = control2_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(control2_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(control2_labels{i})
    
    p = p+1;
    subplot(3, 2, p)
    x = control2_normTrace{i, 1};
    y = control2_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure14.png')
saveas(gcf, 'figure14.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:2
    p = p + 1;
    subplot(2, 2, p)
    x = control3_times{i};
    y1 = control3_cellTrace{i};
    y2 = control3_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    xline(control3_lysis{i}(1), '--k')
    ylim([0 Inf])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence')
    title(control3_labels{i})
    
    p = p + 1;
    subplot(2, 2, p)
    x = control3_normTrace{i,1};
    y = control3_correctedTrace{i,1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])    
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end
saveas(gcf, 'figure15.png')
saveas(gcf, 'figure15.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=3:5
    p = p + 1;
    subplot(3, 2, p)
    x = control3_times{i};
    y1 = control3_cellTrace{i};
    y2 = control3_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    xline(control3_lysis{i}(1), '--k')
    ylim([0 Inf])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence')
    title(control3_labels{i})
    
    p = p + 1;
    subplot(3, 2, p)
    x = control3_normTrace{i,1};
    y = control3_correctedTrace{i,1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])    
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end
saveas(gcf, 'figure16.png')
saveas(gcf, 'figure16.fig')

close all
%% plot the corrections for the experimental traces
cd(dirsave)

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:3
    
    p = p+1;
    subplot(3, 2, p)
    x = exp_times{i};
    y1 = exp_cellTrace{i};
    y2 = exp_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(exp_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(exp_labels{i})
    
    p = p+1;
    subplot(3, 2, p)
    x = exp_normTrace{i, 1};
    y = exp_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(exp_labels{i})
end
saveas(gcf, 'figure17.png')
saveas(gcf, 'figure17.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=4:6
    
    p = p+1;
    subplot(3, 2, p)
    x = exp_times{i};
    y1 = exp_cellTrace{i};
    y2 = exp_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(exp_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(exp_labels{i})
    
    p = p+1;
    subplot(3, 2, p)
    x = exp_normTrace{i, 1};
    y = exp_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(exp_labels{i})
end
saveas(gcf, 'figure18.png')
saveas(gcf, 'figure18.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:3
    
    p = p+1;
    subplot(3, 2, p)
    x = pbs_times{i};
    y1 = pbs_cellTrace{i};
    y2 = pbs_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(pbs_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(pbs_labels{i})
    
    p = p+1;
    subplot(3, 2, p)
    x = pbs_normTrace{i, 1};
    y = pbs_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(pbs_labels{i})
end
saveas(gcf, 'figure19.png')
saveas(gcf, 'figure19.fig')

p = 0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=4:6
    
    p = p+1;
    subplot(3, 2, p)
    x = pbs_times{i};
    y1 = pbs_cellTrace{i};
    y2 = pbs_bgTrace{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xline(pbs_lysis{i}(1), '--k')
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(pbs_labels{i})
    
    p = p+1;
    subplot(3, 2, p)
    x = pbs_normTrace{i, 1};
    y = pbs_correctedTrace{i, 1};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(pbs_labels{i})
end
saveas(gcf, 'figure20.png')
saveas(gcf, 'figure20.fig')

close all

%% now compare the pbs traces
cd(dirsave)

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(pbs_correctedTrace)
    subplot(2, 3, i)
    x = pbs_normTrace{i,1};
    y = pbs_correctedTrace{i};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 1.5]) 
    xlim([-0.5 Inf])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(pbs_labels{i})
end
saveas(gcf, 'figure21.png')
saveas(gcf, 'figure21.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(pbs_correctedTrace)
    subplot(2, 3, i)
    x = pbs_normTrace{i,1};
    y = pbs_correctedTrace{i};
    
    dx = diff(x, 1, 2);
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl./y(:, 1:end-1);
    
    ybar = mean(el, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(el, 0, 1, 'omitnan');
    
    errorbar(x(1:end-1), -ybar, err, 'Color', okabeIto{i})
    xlim([-0.5 Inf])
    %ylim([-0.15 0.06]) 
    ylim([-0.06 0.15]) 
    xlabel('Time (minutes)')
    ylabel('C_{eff}')
    title(pbs_labels{i})
end
saveas(gcf, 'figure22.png')
saveas(gcf, 'figure22.fig')


%% now compare the experimental traces
cd(dirsave)

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_correctedTrace)
    subplot(2, 3, i)
    x = exp_normTrace{i,1};
    y = exp_correctedTrace{i};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 1.5]) 
    xlim([-0.5 Inf])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(exp_labels{i})
end
saveas(gcf, 'figure23.png')
saveas(gcf, 'figure23.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_correctedTrace)
    subplot(2, 3, i)
    x = exp_normTrace{i,1};
    y = exp_correctedTrace{i};
    
    dx = diff(x, 1, 2);
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl./y(:, 1:end-1);
    
    ybar = mean(el, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(el, 0, 1, 'omitnan');
    
    errorbar(x(1:end-1), -ybar, err, 'Color', okabeIto{i})
    xlim([-0.5 Inf])
    %ylim([-0.05 0.03])
    ylim([-0.03 0.05])
    xlabel('Time (minutes)')
    ylabel('C_{eff}')
    title(exp_labels{i})
end
saveas(gcf, 'figure24.png')
saveas(gcf, 'figure24.fig')

%% plot the maximum Ceff and the time
all_labels = [control2_labels, control3_labels, pbs_labels, exp_labels];
control2_maxe = cell(length(control2_labels), 1);
control3_maxe = cell(length(control3_labels), 1);
pbs_maxe = cell(length(pbs_labels), 1);
exp_maxe = cell(length(exp_labels), 1);

for i=1:height(control2_ceff)
    y = control2_ceff{i};
    [control2_maxe{i}(:, 1), control2_maxe{i}(:, 2)] = max(y, [], 2);
end

for i=1:height(control3_ceff)
    y = control3_ceff{i};
    [control3_maxe{i}(:,1), control3_maxe{i}(:, 2)] = max(y, [], 2);
end

for i=1:height(pbs_ceff)
    y = control2_ceff{i};
    [pbs_maxe{i}(:, 1), pbs_maxe{i}(:, 2)] = max(y, [], 2);
end

for i=1:height(exp_ceff)
    y = control2_ceff{i};
    [exp_maxe{i}(:, 1), exp_maxe{i}(:, 2)] = max(y, [], 2);
end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_ceff)
    x = control2_maxe{i}(:, 2);
    y = control2_maxe{i}(:, 1);
    scatter(x, y, 'o', 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
end
for i=1:height(control3_ceff)
    x = control3_maxe{i}(:, 2);
    y = control3_maxe{i}(:, 1);
    scatter(x, y, '^', 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
end
for i=1:height(pbs_ceff)
    x = pbs_maxe{i}(:, 2);
    y = pbs_maxe{i}(:, 1);
    scatter(x, y, 's', 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
end
for i=1:height(exp_ceff)
    x = exp_maxe{i}(:, 2);
    y = exp_maxe{i}(:, 1);
    scatter(x, y, 'v', 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
end
legend(all_labels)
%% see how well tau fits for the pbs and experimental traces
%plot the fit of tau
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:5
    subplot(2, 3, i)
    x = pbs_normTrace{i,1};
    
    y1 = pbs_correctedTrace{i};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = pbs_tau_yhat{i};
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    plot(x, y2, '--k')
    ylim([0 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(pbs_labels{i})
end
saveas(gcf, 'figure25.png')
saveas(gcf, 'figure25.fig')

%plot the fit of tau
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:5
    subplot(2, 3, i)
    x = exp_normTrace{i,1};
    
    y1 = exp_correctedTrace{i};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = exp_tau_yhat{i};
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    plot(x, y2, '--k')
    ylim([0 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(exp_labels{i})
end
saveas(gcf, 'figure26.png')
saveas(gcf, 'figure26.fig')
%% plot smoothed mean length traces
cd(dirsave)

% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
% for i=1:height(exp_lcell)
%     subplot(2, 3, i)
%     x = exp_times{i};
%     y = exp_lcell{i};
%     ybar = mean(y, 1, 'omitnan');
%     ybar = movingaverage(ybar, 3);
%     err = std(y, 0, 1, 'omitnan');
%     
%     dx = diff(x, 1, 2);
%     dy = diff(y, 1, 2);
%     dl = dy./dx;
%     el = dl./y(:, 1:end-1);
%     exp_elongation{i} = el;
%     
%     dy = diff(ybar, 1, 2);
%     dl = dy./dx;
%     el = dl./ybar(:, 1:end-1);
%     
%     ilyse = find(el < -0.01, 1, 'first');
%     flyse = find(el(ilyse:end) > -0.005, 1, 'first') + (ilyse-1);
%     
%     exp_lysis{i}(1, 1) = ilyse;
%     exp_lysis{i}(1, 2) = flyse;
%     idx1 = exp_lysis{i}(1,1);
%     idx2 = exp_lysis{i}(1,2);
%     
%     errorbar(x, ybar, err, 'Color', okabeIto{i})
%     ylim([0 15]) 
%     xlim([-0.5 Inf])
%     xlabel('Time (minutes)')
%     ylabel('Length (\mum)')
%     title(exp_labels{i})
%     xline(x(idx1), '--k')
%     xline(x(idx2), '--k')
% end
% saveas(gcf, 'figure05.png')
% saveas(gcf, 'figure05.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_lcell)
    subplot(3, 3, i)
    x = control2_times{i};
    y = control2_lcell{i};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    dx = diff(x, 1, 2);
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl./y(:, 1:end-1);
    control2_elongation{i} = el;
    
    dy = diff(ybar, 1, 2);
    dl = dy./dx;
    el = dl./ybar(:, 1:end-1);
    
    ilyse = find(el < -0.01, 1, 'first');
    flyse = find(el(ilyse:end) > -0.005, 1, 'first') + (ilyse-1);
    
    control2_lysis{i}(1, 1) = ilyse;
    control2_lysis{i}(1, 2) = flyse;
    idx1 = control2_lysis{i}(1,1);
    idx2 = control2_lysis{i}(1,2);
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 15])    
    xlim([-0.5 Inf])
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    title(control2_labels{i})
    xline(x(idx1), '--k')
    xline(x(idx2), '--k')
end
% saveas(gcf, 'figure06.png')
% saveas(gcf, 'figure06.fig')
% 
% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
% for i=1:height(control3_lcell)
%     
%     if i==5
%         subplot(2, 3, [i, i+1])
%     else
%         subplot(2, 3, i)
%     end
%     
%     x = control3_times{i};
%     y = control3_lcell{i};
%     ybar = mean(y, 1, 'omitnan');
%     ybar = movingaverage(ybar, 3);
%     err = std(y, 0, 1, 'omitnan');
%     
%     dx = diff(x, 1, 2);
%     dy = diff(y, 1, 2);
%     dl = dy./dx;
%     el = dl./y(:, 1:end-1);
%     control3_elongation{i} = el;
%     
%     dy = diff(ybar, 1, 2);
%     dl = dy./dx;
%     el = dl./ybar(:, 1:end-1);
%     
%     ilyse = find(el < -0.01, 1, 'first');
%     flyse = find(el(ilyse:end) > -0.005, 1, 'first') + (ilyse-1);
%     
%     control3_lysis{i}(1, 1) = ilyse;
%     control3_lysis{i}(1, 2) = flyse;
%     idx1 = control3_lysis{i}(1,1);
%     idx2 = control3_lysis{i}(1,2);
%     
%     errorbar(x, ybar, err, 'Color', okabeIto{i})
%     ylim([0 20])
%     xlim([-0.5 Inf])
%     xlabel('Time (minutes)')
%     ylabel('Length (\mum)')
%     title(control3_labels{i})
%     xline(x(idx1), '--k')
%     xline(x(idx2), '--k')
% end
% saveas(gcf, 'figure07.png')
% saveas(gcf, 'figure07.fig')

%% determining change in fluor
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_adjTrace)
    subplot(2, 3, i)
    x = exp_times{i};
    y = exp_adjTrace{i};
    
    dx = diff(x, 1, 2);
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl./y(:, 1:end-1);
    ybar = mean(el, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(el, 0, 1, 'omitnan');
%     exp_elongation{i} = el;
    
    errorbar(x(1:end-1), ybar, err, 'Color', okabeIto{i})
%     ylim([0 15]) 
%     xlim([-0.5 Inf])
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    title(exp_labels{i})
%     xline(x(idx1), '--k')
%     xline(x(idx2), '--k')
end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_adjTrace)
    subplot(3, 3, i)
    x = control2_times{i};
    y = control2_adjTrace{i};
    
    dx = diff(x, 1, 2);
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl./y(:, 1:end-1);
    ybar = mean(el, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(el, 0, 1, 'omitnan');
%     control2_elongation{i} = el;
    
    errorbar(x(1:end-1), ybar, err, 'Color', okabeIto{i})
%     ylim([0 15]) 
%     xlim([-0.5 Inf])
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    title(control2_labels{i})
%     xline(x(idx1), '--k')
%     xline(x(idx2), '--k')
end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control3_adjTrace)
    
    if i==5
        subplot(2, 3, [i, i+1])
    else
        subplot(2, 3, i)
    end
    
    x = control3_times{i};
    y = control3_adjTrace{i};
    
    dx = diff(x, 1, 2);
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl./y(:, 1:end-1);
    ybar = mean(el, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(el, 0, 1, 'omitnan');
%     control3_elongation{i} = el;
    
    errorbar(x(1:end-1), ybar, err, 'Color', okabeIto{i})
%     ylim([0 15]) 
%     xlim([-0.5 Inf])
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    title(control3_labels{i})
%     xline(x(idx1), '--k')
%     xline(x(idx2), '--k')
end

%% plot the fit of alpha for all the traces
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_normTrace)
    subplot(2, 3, i)
    x = exp_normTrace{i,1};
    
    y1 = exp_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = exp_alpha_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([-Inf 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(exp_labels{i})
end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_normTrace)
    subplot(3, 3, i)
    x = control2_normTrace{i,1};
    
    y1 = control2_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control2_alpha_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([-Inf 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end


figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control3_normTrace)
    
    if i==5
        subplot(2, 3, [i, i+1])
    else
        subplot(2, 3, i)
    end
    
    x = control3_normTrace{i,1};
    
    y1 = control3_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = control3_alpha_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
    ylim([-Inf 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end

%% plot corrected traces
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:2
    subplot(1, 2, i)
    x = control2_normTrace{i, 1};
    y = control2_correctedTrace{i, 1};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 Inf])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure16.png')
saveas(gcf, 'figure16.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:2
    subplot(1, 2, i)
    x = control3_normTrace{i, 1};
    y = control3_correctedTrace{i, 1};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 Inf])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end
saveas(gcf, 'figure17.png')
saveas(gcf, 'figure17.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_normTrace)
    subplot(3, 3, i)
    x = control2_normTrace{i, 1};
    y = control2_correctedTrace{i, 1};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 Inf])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
saveas(gcf, 'figure18.png')
saveas(gcf, 'figure18.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control3_normTrace)
    
    if i==5
        subplot(2, 3, [i, i+1])
    else
        subplot(2, 3, i)
    end
    
    x = control3_normTrace{i, 1};
    y = control3_correctedTrace{i, 1};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 Inf])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end
saveas(gcf, 'figure19.png')
saveas(gcf, 'figure19.fig')

%% plot normalized traces
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_normTrace)
    subplot(2, 3, i)
    x = exp_normTrace{i, 1};
    y = exp_normTrace{i, 2};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([-Inf 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(exp_labels{i})
end
% saveas(gcf, 'figure06.png')
% saveas(gcf, 'figure06.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_normTrace)
    subplot(3, 3, i)
    x = control2_normTrace{i, 1};
    y = control2_normTrace{i, 2};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([-Inf 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
% saveas(gcf, 'figure04.png')
% saveas(gcf, 'figure04.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control3_normTrace)
    
    if i==5
        subplot(2, 3, [i, i+1])
    else
        subplot(2, 3, i)
    end
    
    x = control3_normTrace{i, 1};
    y = control3_normTrace{i, 2};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([-Inf 1.1])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end
% saveas(gcf, 'figure05.png')
% saveas(gcf, 'figure05.fig')

%% plot normalized traces over frame number
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_normTrace)
    subplot(2, 3, i)
    [~, ncol] = size(exp_normTrace{i, 1});
    x = 1:ncol;
    y = exp_normTrace{i, 2};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([-Inf 1.1])
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence')
    title(exp_labels{i})
end
% saveas(gcf, 'figure06.png')
% saveas(gcf, 'figure06.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_normTrace)
    subplot(3, 3, i)
    [~, ncol] = size(control2_normTrace{i, 1});
    x = 1:ncol;
    y = control2_normTrace{i, 2};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([-Inf 1.1])
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence')
    title(control2_labels{i})
end
% saveas(gcf, 'figure04.png')
% saveas(gcf, 'figure04.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control3_normTrace)
    
    if i==5
        subplot(2, 3, [i, i+1])
    else
        subplot(2, 3, i)
    end
    
    [~, ncol] = size(control3_normTrace{i, 1});
    x = 1:ncol;
    y = control3_normTrace{i, 2};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([-Inf 1.1])
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end
% saveas(gcf, 'figure05.png')
% saveas(gcf, 'figure05.fig')

%% Functions
function [cellTrace, bgTrace, adjTrace, times, lCell, frameRate] = loadData(dirpath, basenames)

    cellTrace = cell(length(basenames), 1);
    bgTrace = cell(length(basenames), 1);
    adjTrace = cell(length(basenames), 1);
    times = cell(length(basenames), 1);
    lCell = cell(length(basenames), 1);
    frameRate = cell(length(basenames), 1);

    for i=1:length(basenames)

        cd(dirpath);
        basename = basenames{i};
        datadir=dir([basename '*']);
        load(datadir(1).name, '-regexp', 'intensity$', 'time', 'lcell');

        cellTrace{i} = intensity;
        bgTrace{i} = bgintensity;
        adjTrace{i} = adjintensity;
        times{i} = time;
        lCell{i} = lcell;
        frameRate{i} = [0, diff(time, 1, 2)];
        clear intensity adjintensity bgintensity time

    end

end

function y = movingaverage(x,s)

        [mx,lx]=size(x);
        y=zeros(mx,lx);

        if s>lx
            s=round(lx/2);
        end

        for i=1:lx
                s2=floor(s/2);
            if i-s2<=0
                y(:,i)=mean(x(:,1:i+s2), 2, 'omitnan');
            elseif i+s2>lx
                y(:,i)=mean(x(:,i-s2:end), 2, 'omitnan');
            else
                y(:,i)=mean(x(:,i-s2:i+s2), 2, 'omitnan');
            end
        end
    
end

function [lysis, elongation] = lysisCalc(times, lcell)

    lysis = cell(size(lcell));
    elongation = cell(size(lcell));

    for i=1:height(lcell)
        x = times{i};
        y = lcell{i};
        ybar = mean(y, 1, 'omitnan');
        ybar = movingaverage(ybar, 3);

        %calculate the elongation rate of each length trace
        dx = diff(x, 1, 2);
        dy = diff(y, 1, 2);
        dl = dy./dx;
        el = dl./y(:, 1:end-1);
        elongation{i} = el;

        %calculate the elongation rate of the smoothed averaged length trace
        dy = diff(ybar, 1, 2);
        dl = dy./dx;
        el = dl./ybar(:, 1:end-1);

        %find the initial and final frame of detergent perfusion
        ilyse = find(el < -0.01, 1, 'first');
        flyse = find(el(ilyse:end) > -0.005, 1, 'first') + (ilyse-1);

        lysis{i}(1, 1) = ilyse;
        lysis{i}(1, 2) = flyse;
    end

end

function [normTrace] = dataNormalize(times, adjTrace, lysis)
    
    %pre-allocate variables
    [nrow, ~] = size(adjTrace);
    normTrace = cell(nrow, 2);
    
    for i=1:height(adjTrace)  
        
        adjintensity = adjTrace{i};
        tme = times{i};
        
        %interpolate values obtained during detergent perfusion
        [nrow, ncol] = size(adjintensity);
        prelysis = lysis{i}(1,1)-1;
        xq = lysis{i}(1,1):lysis{i}(1,2)+2;
        x = setdiff(1:ncol, xq);
        for j=1:nrow
            vq = interp1(tme(x), adjintensity(j, x), tme(xq));
            adjintensity(j, xq) = vq;
        end
        
        %subtract final fluor value (if that value is reached)    
        dx = diff(tme, 1, 2);
        dy = diff(adjintensity, 1, 2);
        df = dy./dx;
        ef = df./adjintensity(:, 1:end-1);
        ybar = mean(ef, 1, 'omitnan');
        ybar = movingaverage(ybar, 3);
        
%         if ybar(end) > -0.005
%             
%             [~, minx] = min(ybar);
%             fidx = find(ybar(minx:end) > -0.005, 1, 'first') + (minx-1);
%             
%             if isempty(fidx)
%                 fidx = ncol;
%             end
%             
%             %adjintensity = adjintensity - adjintensity(:, fidx);
%             
%         else
%             [~, minx] = min(ybar);
%             fidx = find(ybar(minx:end) > -0.005, 1, 'first') + (minx-1);
%             
%             if isempty(fidx)
%                 fidx = ncol;
%             %else
%                 %adjintensity = adjintensity - adjintensity(:, fidx);
%             end
%         end
        
        %adjust the time vector
        normTrace{i,1} = tme(prelysis:end) - tme(prelysis);
        
        %determine the frame number
%         frame = normTrace{i,1}(2:end)./diff(normTrace{i,1}, 1, 2);
%         normTrace{i, 2} = [1, frame+1];
        normTrace{i, 2} = adjintensity(:, prelysis:end) ./ adjintensity(:, prelysis);
        
    end
    
end

function [alpha, alpha_yhat] = alphaCalc(times, normTrace, alpha0)
    
    %pre-allocate variables
    alpha = cell(height(normTrace), 1);
    alpha_yhat = cell(height(normTrace), 1);
    
    for i=1:height(normTrace)

        normintensity = normTrace{i};
        tme = times{i};
        
        %find out if the frame rate changes
        dt = diff(tme, 1, 2);
        
        %convert negative values to NaN
        normintensity(normintensity<0) = NaN;
    
        %index traces that have at least 1 data point
        [nrow, ncol]=size(normintensity); 
        idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
        normintensity=normintensity(idx, :);
        
        %remove traces without a start value
        if dt(1)~=dt(end)
            imstart = find(dt == dt(end), 1, 'first');
        else 
            imstart = 1;
        end
        
        idx = find(~isnan(normintensity(:, imstart))); 
        normintensity = normintensity(idx, :);
        
        %pre-allocate output variables
        [nrow, ~]=size(normintensity);
        alphaTemp=nan(nrow,1);    
        yhatTemp=nan(size(normintensity)); 
    
        %define exponential function
        %modelfun = @(param, t)param(1)*exp(-t./param(2));
        
        if dt(1)~=dt(end)
            idx = find(dt == dt(end), 1, 'first');
            xq = idx:ncol;
            vq = setdiff(1:ncol, xq);
            
            %fit data to the function and evaluate
            for j=1:nrow
                A = normintensity(j,idx);
                modelfun = @(alpha, t)A*exp(-t./alpha);
                alphaTemp(j, 1)=nlinfit(tme(xq)-tme(idx), normintensity(j,xq), modelfun, alpha0);
                yhatTemp(j, xq)=modelfun(alphaTemp(j), tme(xq)-tme(idx));
                yhatTemp(j, vq)=NaN;
            end
        else
            %fit data to the function and evaluate
            modelfun = @(alpha, t)exp(-t./alpha);
            for j=1:nrow
                alphaTemp(j, 1)=nlinfit(tme, normintensity(j,:), modelfun, alpha0);
                yhatTemp(j, :)=modelfun(alphaTemp(j), tme);
            end            
        end
   
        
        alpha{i} = alphaTemp;
        alpha_yhat{i} = yhatTemp;
        
    end    
    
end

function [alphaCoef, alphaFit] = alphaLinear(frameRate, alpha)
    
    dx = [];
    rx = [];
    xv = [];
    yv = [];

    for i=1:height(alpha)
        x = frameRate{i}(end);
        y = alpha{i}';
        dx = [dx, x];

        rx = repelem(x, length(y));
        xv = [xv, rx];
        yv = [yv, y];
    end

    alphaCoef = polyfit(xv, yv, 1);
    alphaFit = polyval(alphaCoef,[0 dx]);

end

function [rho, rho_yhat] = rhoCalc(times, normTrace, rho0, positions)
    
    %pre-allocate variables
    rho = cell(height(normTrace), 1);
    rho_yhat = cell(height(normTrace), 1);
    
    for i=1:height(normTrace)
        
        normintensity = normTrace{i};
        tme = times{i};
       
        %find out if the frame rate changes
        dt = diff(tme, 1, 2);
        
        %convert negative values to NaN
        normintensity(normintensity<0) = NaN;
    
        %index traces that have at least 1 data point
        [nrow, ncol]=size(normintensity); 
        idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
        normintensity=normintensity(idx, :);
        frames = [1:ncol] * positions;
        
        if dt(1)~=dt(end)
            imstart = find(dt == dt(end), 1, 'first');
        else 
            imstart = 1;
        end
        
        idx = find(~isnan(normintensity(:, imstart))); 
        normintensity = normintensity(idx, :);
        
        %pre-allocate output variables
        [nrow, ~]=size(normintensity);
        rhoTemp=nan(nrow,1);    
        yhatTemp=nan(size(normintensity));
    
        %define exponential function
        %modelfun = @(param, x)param(1)*exp(-x*param(2));
        
        if dt(1)~=dt(end)
            idx = find(dt == dt(end), 1, 'first');
            xq = idx:ncol;
            vq = setdiff(1:ncol, xq);
            frames = (xq-(idx-1))*4;
            
            %fit data to the function and evaluate
            for j=1:nrow
                A = normintensity(j,idx);
                modelfun = @(rho, x)A*exp(-x*rho);
                rhoTemp(j, 1)=nlinfit(frames, normintensity(j,xq), modelfun, rho0);
                yhatTemp(j, xq)=modelfun(rhoTemp(j), frames);
                yhatTemp(j, vq) = NaN;
            end
        else
            %fit data to the function and evaluate
            modelfun = @(rho, x)exp(-x*rho);
            for j=1:nrow
                rhoTemp(j, 1)=nlinfit(frames, normintensity(j,:), modelfun, rho0);
                yhatTemp(j, :)=modelfun(rhoTemp(j, 1), frames);
            end          
        end
                 
        rho{i} = rhoTemp;
        rho_yhat{i} = yhatTemp;  
        
    end  
    
end

function [rhoCoef, rhoFit] = rhoLinear(frameRate, rho)
    
    dx = [];
    rx = [];
    xv = [];
    yv = [];

    for i=1:height(rho)
        x = frameRate{i}(end);
        y = rho{i}';
        dx = [dx, x];

        rx = repelem(x, length(y));
        xv = [xv, rx];
        yv = [yv, y];
    end

    rhoCoef = polyfit(xv, yv, 1);
    rhoFit = polyval(rhoCoef,[0 dx]);

end

function [correctedTrace]=photoCorrect(times, normTrace, omega, intercept, positions)
        
        correctedTrace = cell(height(normTrace), 6);
        
        for h=1:height(normTrace)
            
            normintensity = normTrace{h};
            tme = times{h};
            
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
       dC=@(C, omega, dt, b, positions)(C*((omega*dt+b)^positions));
       %dC=@(C, alpha, dt, b)((C*exp(-dt/(alpha*dt+b)))*(1/(alpha*dt+b))*dt);
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                dCB(n,i) = dC(normintensity(n,i), omega, dt(i), intercept, positions); %this is the amount of photobleaching that occured in our measured value
   
                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)*unb_frac(n,i); %Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
        end
        
        correctedTrace{h, 1} = Cnew;
        correctedTrace{h, 2} = dCB;
        correctedTrace{h, 3} = dCT;
        correctedTrace{h, 4} = dCP;
        correctedTrace{h, 5} = NaN; %Cbl_exp;
        correctedTrace{h, 6} = NaN; %unb_frac;
        
        end
        
end

function [Ceff] = effectiveLoss(times, correctedTrace)

    [nrow, ~] = size(correctedTrace);
    Ceff= cell(nrow, 1);
    
    for i=1:nrow
        x = times{i};
        y = correctedTrace{i};
        
        dx = diff(x, 1, 2);
        dy = diff(y, 1, 2);
        dl = dy./dx;
        Ceff{i} = dl./y(:, 1:end-1);      
    end

end

function [tau, tau_yhat] = tauCalc(times, correctedTrace, tau0)
    
    tau = cell(height(correctedTrace), 1);
    tau_yhat = cell(height(correctedTrace), 1);
    
    for i = 1:height(correctedTrace)
        x = times{i};
        y = correctedTrace{i};
          
        ybar = mean(y, 1, 'omitnan');
        A = ybar(1);
        B = ybar(end);
        
        %modelfun = @(tau, t) A * exp(-t./tau) + B;
        modelfun = @(tau, t) A * exp(log(A/B)*(1-exp(t./tau)));
        
        tau{i} = nlinfit(x, ybar, modelfun, tau0);
        tau_yhat{i} = modelfun(tau{i}, x);
    end
    
end