%sanity_checks.m
%Author: Zarina Akbary
%Date: 04/19/2022
%Purpose: this is basically a sandbox for me to look through my data and
%keep track of how I do it. 

clear, close all

okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

%% load data
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/aggregate/';
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/figures/';

exp_basenames = {'04042022_Exp1', '04112022_Exp1', '04112022_Exp3', '04132022_Exp1', '04132022_Exp2', '04142022_Exp2'};
exp_labels = {'untreated', 'untreated', 'glycerol (5 minutes)', 'PBS (5 minutes)', 'spent LB (5 minutes)', 'glucose (5 minutes)'};

control2_basenames = {'04052022_Exp2', '04052022_Exp1', '04042022_Exp1', '04112022_Exp1', '04052022_Exp3', '02122022_Exp1', '02122022_Exp2', '02092022_Exp1', '01282022_Exp1'};
control2_labels = {'Frame Rate = 10 s', 'Frame Rate = 20 s', 'Frame Rate = 1 min', 'Frame Rate = 1 min', 'Frame Rate = 2 min', 'Frame Rate = 5 min', 'Frame Rate = 10 min', 'Frame Rate = 20 min', 'Frame Rate = 20 min'};

control3_basenames = {'10302021_Exp2', '10302021_Exp1', '10232021_Exp1', '10262021_Exp1', '10282021_Exp1'};
control3_labels = {'Frame Rate = 15 s', 'Frame Rate = 30 s', 'Frame Rate = 1 min', 'Frame Rate = 1 min', 'Frame Rate = 5 min'};

[exp_cellTrace, exp_bgTrace, exp_adjTrace, exp_times, exp_lcell, exp_frameRate] = loadData(dirpath, exp_basenames);
[control2_cellTrace, control2_bgTrace, control2_adjTrace, control2_times, control2_lcell, control2_frameRate] = loadData(dirpath, control2_basenames);
[control3_cellTrace, control3_bgTrace, control3_adjTrace, control3_times, control3_lcell, control3_frameRate] = loadData(dirpath, control3_basenames);

%% normalize traces
[exp_normTrace] = dataNormalize(exp_times, exp_adjTrace, exp_lysis);
[control2_normTrace] = dataNormalize(control2_times, control2_adjTrace, control2_lysis);
[control3_normTrace] = dataNormalize(control3_times, control3_adjTrace, control3_lysis);

%% calculate alpha and rho
alpha0 = 1;
rho0 = 0;

[exp_alpha, exp_alpha_yhat] = alphaCalc(exp_normTrace(:, 1), exp_normTrace(:, 2), alpha0);
[control2_alpha, control2_alpha_yhat] = alphaCalc(control2_normTrace(:, 1), control2_normTrace(:, 2), alpha0);
[control3_alpha, control3_alpha_yhat] = alphaCalc(control3_normTrace(:, 1), control3_normTrace(:, 2), alpha0);

[exp_rho, exp_rho_yhat] = rhoCalc(exp_normTrace(:, 1), exp_normTrace(:, 2), rho0);
[control2_rho, control2_rho_yhat] = rhoCalc(control2_normTrace(:, 1), control2_normTrace(:, 2), rho0);
[control3_rho, control3_rho_yhat] = rhoCalc(control3_normTrace(:, 1), control3_normTrace(:, 2), rho0);

%% fit rho and alpha to linear equations (without the 1 minute data)
dx = [control2_frameRate{1}(end), control2_frameRate{2}(end)];
x = [repelem(dx(1), length(control2_rho{1})), repelem(dx(2),length(control2_rho{2}))];

y = [control2_rho{1}', control2_rho{2}'];
control2_rhoCoef = polyfit(x, y, 1);
control2_rhoFit = polyval(control2_rhoCoef,[0 dx]);

y = [control2_alpha{1}', control2_alpha{2}'];
control2_alphaCoef = polyfit(x, y, 1);
control2_alphaFit = polyval(control2_alphaCoef,[0 dx]);

dx = [control3_frameRate{1}(end), control3_frameRate{2}(end)];
x = [repelem(dx(1), length(control3_rho{1})), repelem(dx(2),length(control3_rho{2}))];

y = [control3_rho{1}', control3_rho{2}'];
control3_rhoCoef = polyfit(x, y, 1);
control3_rhoFit = polyval(control3_rhoCoef,[0 dx]);

y = [control3_alpha{1}', control3_alpha{2}'];
control3_alphaCoef = polyfit(x, y, 1);
control3_alphaFit = polyval(control3_alphaCoef,[0 dx]);

%% correct for photobleaching
[exp_correctedTrace]=photoCorrect(exp_normTrace(:, 1), exp_normTrace(:, 2), control2_alphaCoef(1), control2_alphaCoef(2));
[control2_correctedTrace]=photoCorrect(control2_normTrace(:, 1), control2_normTrace(:, 2), control2_alphaCoef(1), control2_alphaCoef(2));
[control3_correctedTrace]=photoCorrect(control3_normTrace(:, 1), control3_normTrace(:, 2), control3_alphaCoef(1), control3_alphaCoef(2));

%% plot the fit of alpha
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
% saveas(gcf, 'figure06.png')
% saveas(gcf, 'figure06.fig')

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
% saveas(gcf, 'figure04.png')
% saveas(gcf, 'figure04.fig')

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

%% plot the fit of rho
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_normTrace)
    subplot(2, 3, i)
    [~, ncol] = size(exp_normTrace{i, 1});
    x = 1:ncol;
    
    y1 = exp_normTrace{i, 2};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar1 = movingaverage(ybar1, 3);
    err1 = std(y1, 0, 1, 'omitnan');
    
    y2 = exp_rho_yhat{i};
    ybar2 = mean(y2, 1, 'omitnan');
    ybar2 = movingaverage(ybar2, 3);
    err2 = std(y2, 0, 1, 'omitnan');
    
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, '--k')
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
    ylim([-Inf 1.1])
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence')
    title(control3_labels{i})
end

%% plot alpha vs frame rate for controls 2 & 3
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_alpha)
    
    x = control2_frameRate{i}(end);
    
    y = control2_alpha{i};
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')
    xlabel('Frame Rate')
    ylabel('\alpha (min)')
end
legend(control2_labels)
    
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control3_alpha)-1
    
    x = control3_frameRate{i}(end);
    
    y = control3_alpha{i};
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')
    xlabel('Frame Rate')
    ylabel('\alpha (min)')
end
legend(control3_labels)

%% plot rho vs frame rate for controls 2 & 3
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:2
    
    x = control2_frameRate{i}(end);
    
    y = control2_rho{i};
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')
end
xlabel('Frame Rate')
ylabel('\rho')
%legend(control2_labels)
    
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control3_rho)-3
    
    x = control3_frameRate{i}(end);
    
    y = control3_rho{i};
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')
end
xlabel('Frame Rate')
ylabel('\rho')
%legend(control3_labels)


%% plot corrected traces
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
% saveas(gcf, 'figure05.png')
% saveas(gcf, 'figure05.fig')

%% plot smoothed mean length traces and use it to calculate the detergent frames
cd(dirsave)

exp_lysis = cell(size(exp_lcell));
control2_lysis = cell(size(control2_lcell));
control3_lysis = cell(size(control3_lcell));

exp_elongation = cell(size(exp_lcell));
control2_elgonation = cell(size(control2_lcell));
control3_elongation = cell(size(control3_lcell));

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_lcell)
    subplot(2, 3, i)
    x = exp_times{i};
    y = exp_lcell{i};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    dx = diff(x, 1, 2);
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl./y(:, 1:end-1);
    exp_elongation{i} = el;
    
    dy = diff(ybar, 1, 2);
    dl = dy./dx;
    el = dl./ybar(:, 1:end-1);
    
    ilyse = find(el < -0.01, 1, 'first');
    flyse = find(el(ilyse:end) > -0.005, 1, 'first') + (ilyse-1);
    
    exp_lysis{i}(1, 1) = ilyse;
    exp_lysis{i}(1, 2) = flyse;
    idx1 = exp_lysis{i}(1,1);
    idx2 = exp_lysis{i}(1,2);
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 15]) 
    xlim([-0.5 Inf])
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    title(exp_labels{i})
    xline(x(idx1), '--k')
    xline(x(idx2), '--k')
end
saveas(gcf, 'figure06.png')
saveas(gcf, 'figure06.fig')

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
saveas(gcf, 'figure04.png')
saveas(gcf, 'figure04.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control3_lcell)
    
    if i==5
        subplot(2, 3, [i, i+1])
    else
        subplot(2, 3, i)
    end
    
    x = control3_times{i};
    y = control3_lcell{i};
    ybar = mean(y, 1, 'omitnan');
    ybar = movingaverage(ybar, 3);
    err = std(y, 0, 1, 'omitnan');
    
    dx = diff(x, 1, 2);
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl./y(:, 1:end-1);
    control3_elongation{i} = el;
    
    dy = diff(ybar, 1, 2);
    dl = dy./dx;
    el = dl./ybar(:, 1:end-1);
    
    ilyse = find(el < -0.01, 1, 'first');
    flyse = find(el(ilyse:end) > -0.005, 1, 'first') + (ilyse-1);
    
    control3_lysis{i}(1, 1) = ilyse;
    control3_lysis{i}(1, 2) = flyse;
    idx1 = control3_lysis{i}(1,1);
    idx2 = control3_lysis{i}(1,2);
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 20])
    xlim([-0.5 Inf])
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    title(control3_labels{i})
    xline(x(idx1), '--k')
    xline(x(idx2), '--k')
end
saveas(gcf, 'figure05.png')
saveas(gcf, 'figure05.fig')


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
%% plot smooth length traces vs time
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_lcellSmooth)
    subplot(2, 3, i)
    x = exp_times{i};
    y = exp_lcellSmooth{i};    
    plot(x, y, 'Color', okabeIto{i})
    title(exp_labels{i})
    ylim([0 17])
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_lcellSmooth)
    subplot(2, 3, i)
    x = exp_times{i};
    y = exp_lcellSmooth{i};
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    
    idx1 = exp_lysis{i}(1,1);
    idx2 = exp_lysis{i}(1,2);
    
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    ylim([0 15])    
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    title(exp_labels{i})
    xline(x(idx1), '--k')
    xline(x(idx2), '--k')
end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(exp_lcellSmooth)
    subplot(2, 3, i)
    x = exp_times{i};
    dx = diff(x, 1, 2);
    y = exp_lcellSmooth{i};
    dy = diff(y, 1, 2);
    dl = dy./dx;
    el = dl ./ y(:, 1:end-1);
    ebar = mean(el, 1, 'omitnan');
    err = std(el, 0, 1, 'omitnan');
    errorbar(x(1:end-1), ebar, err, 'Color', okabeIto{i})  
    ylim([-0.06 0.08])
    xlabel('Time (minutes)')
    ylabel('$\dot{e} (min^{-1})$', 'Interpreter', 'latex')
    title(exp_labels{i})
end



%% plot fluor. trace over time
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(cellTrace)
    subplot(2, 3, i)
    x = times{i};
    y = cellTrace{i};
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    title(labels{i})
end

%% plot the normalized traces
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(cellTrace)
    subplot(2, 3, i)
    x = times{i};
    y = cellTrace{i};
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i})
    title(labels{i})
end
%% plot the effective elongation rates over time
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(normTrace)
    x = normTrace{i, 1};
    y = normTrace{i, 2};
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i})
end
legend(labels)
%% plot the elongation rates for the first five minutes of perfusion
phase1 = 4;
phase2 = 9;

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(lcell)
%     x = times{i}(phase1);
    y = lcell{i}(:, phase1:phase1+1);
    dy = diff(y, 1, 2);
    ydot = dy ./ lcell{i}(:, phase1);
    
    h = histogram(ydot, 'Normalization', 'probability', 'FaceColor', okabeIto{i}, 'NumBins', height(ydot), 'BinWidth', 0.0015, 'EdgeColor', okabeIto{i})
    
%     ybar = [0, mean(ydot, 1, 'omitnan')];
%     err = [0, std(ydot, 0, 1, 'omitnan')];
%     errorbar(x, ybar, err, 'Color', okabeIto{i})
end
legend(labels)

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(lcell)
    if ismember(i, [1:2])
        y = lcell{i}(:, phase1:phase1+1);
        dy = diff(y, 1, 2);
        ydot = dy ./ lcell{i}(:, phase1);
    else
        y = lcell{i}(:, phase2:phase2+1);
        dy = diff(y, 1, 2);
        ydot = dy ./ lcell{i}(:, phase2);        
    end
    
    h = histogram(ydot, 'Normalization', 'probability', 'FaceColor', okabeIto{i}, 'NumBins', height(ydot), 'BinWidth', 0.0015, 'EdgeColor', okabeIto{i})
    
%     ybar = [0, mean(ydot, 1, 'omitnan')];
%     err = [0, std(ydot, 0, 1, 'omitnan')];
%     errorbar(x, ybar, err, 'Color', okabeIto{i})
end
legend(labels)
%% focus on the 1-minute frame rate traces
control_normTrace = cell(4, 3);

phase1 = 5; %end of initial perfusion
omit = 6:10; %detergent perfusion
phase3 = 11; %final perfusion

control_cellTrace = [control2_cellTrace(3:4); control3_cellTrace(3:4)];
control_bgTrace = [control2_bgTrace(3:4); control3_bgTrace(3:4)];
control_adjTrace = [control2_adjTrace(3:4); control3_adjTrace(3:4)];
control_times = [control2_times(3:4); control3_times(3:4)];
control_frameRate = [control2_frameRate(3:4); control3_frameRate(3:4)];
control_labels = [control2_labels(3:4), control3_labels(3:4)];

for i=1:height(control_adjTrace)
    [~, ncol] = size(control_times{i});
    idx = setdiff(phase1:ncol, omit);   
    x = control_times{i}(idx);
    control_normTrace{i, 1} = x - x(1);  
    y = control_adjTrace{i}(:, idx) - control_adjTrace{i}(:, end);
    control_normTrace{i, 2} = y./y(:, 1);   
    control_normTrace{i, 3} = control_frameRate{i}(idx);
end

%% plot the normalized traces
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control_normTrace)
    x = control_normTrace{i, 1};
    y = control_normTrace{i, 2};
    ybar = mean(control_normTrace{i, 2}, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i})
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
legend(control_labels)

%% calculate rho and alpha
alpha0 = 1;
rho0 = 0;

[alpha, alpha_yhat] = alphaCalc(control_normTrace(:, 1), control_normTrace(:, 2), alpha0);
[rho, rho_yhat, frames] = rhoCalc(control_normTrace(:, 1), control_normTrace(:, 3), control_normTrace(:, 2), rho0);

%% plot rho and alpha
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(alpha)
    subplot(2, 2, i)
    x = control_normTrace{i, 1};
    y1 = control_normTrace{i, 2};
    y2 = alpha_yhat{i};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar2 = mean(y2, 1, 'omitnan');
    err1 = std(y1, 0, 1, 'omitnan');
    err2 = std(y2, 0, 1, 'omitnan');
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, 'LineStyle', '--', 'Color', 'black')
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title(control_labels{i})
end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(rho)
    subplot(2, 2, i)
    x = frames{i};
    y1 = control_normTrace{i, 2};
    y2 = rho_yhat{i};
    ybar1 = mean(y1, 1, 'omitnan');
    ybar2 = mean(y2, 1, 'omitnan');
    err1 = std(y1, 0, 1, 'omitnan');
    err2 = std(y2, 0, 1, 'omitnan');
    errorbar(x, ybar1, err1, 'Color', okabeIto{i}), hold on
    errorbar(x, ybar2, err2, 'LineStyle', '--', 'Color', 'black')
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title(control_labels{i})
end

%%
%compare the background fluorescence
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=3:4
    x = control2_times{i};
    y = mean(control2_bgTrace{i}, 1, 'omitnan');
    err = std(control2_bgTrace{i}, 0, 1, 'omitnan');
    errorbar(x, y, err, 'Color', okabeIto{i+1})
    xlabel('Time (minutes)')
    ylabel('Background Fluorescence (A.U.)')
end

for i=3:4
    x = control3_times{i};
    y = mean(control3_bgTrace{i}, 1, 'omitnan');
    err = std(control3_bgTrace{i}, 0, 1, 'omitnan');
    errorbar(x, y, err, 'Color', okabeIto{i+3})
    xlabel('Time (minutes)')
    ylabel('Background Fluorescence (A.U.)')
end
ylim([0 Inf])
legend([control2_labels(3:4), control3_labels(3:4)])

%compare the cellular traces
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=3:4
    x = control2_times{i};
    y = mean(control2_cellTrace{i}, 1, 'omitnan');
    err = std(control2_cellTrace{i}, 0, 1, 'omitnan');
    errorbar(x, y, err, 'Color', okabeIto{i+1})
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
end

for i=3:4
    x = control3_times{i};
    y = mean(control3_cellTrace{i}, 1, 'omitnan');
    err = std(control3_cellTrace{i}, 0, 1, 'omitnan');
    errorbar(x, y, err, 'Color', okabeIto{i+3})
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
end
ylim([0 Inf])
legend([control2_labels(3:4), control3_labels(3:4)])

%compare cellular traces minus background
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=3:4
    x = control2_times{i};
    y = mean(control2_adjTrace{i}, 1, 'omitnan');
    err = std(control2_adjTrace{i}, 0, 1, 'omitnan');
    errorbar(x, y, err, 'Color', okabeIto{i+1})
    xlabel('Time (minutes)')
    ylabel('Cellular Fluorescence - Background (A.U.)')
end

for i=3:4
    x = control3_times{i};
    y = mean(control3_adjTrace{i}, 1, 'omitnan');
    err = std(control3_adjTrace{i}, 0, 1, 'omitnan');
    errorbar(x, y, err, 'Color', okabeIto{i+3})
    xlabel('Time (minutes)')
    ylabel('Cellular Fluorescence - Background (A.U.)')
end
ylim([0 Inf])
legend([control2_labels(3:4), control3_labels(3:4)])

%% subtract the final fluor. value and plot phase 1 vs time and phase 3 vs time
phase1 = 5; %end of initial perfusion
omit = 6:10; %detergent perfusion
phase3 = 11; %final perfusion

%Phase 1: compare cellular traces minus background minus final fluor. value
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=3:4
    x = control2_times{i}(1:phase1);
    y = control2_adjTrace{i}(:, 1:phase1) - control2_adjTrace{i}(:, end);
    y = y./y(:, 1);
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i+1})
end

for i=3:4    
    x = control3_times{i}(1:phase1);
    y = control3_adjTrace{i}(:, 1:phase1) - control3_adjTrace{i}(:, end);
    y = y./y(:, 1);
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i+3})
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
ylim([0 Inf])
legend([control2_labels(3:4), control3_labels(3:4)])

%Phase 3: compare cellular traces minus background minus final fluor. value
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=3:4
    [~, ncol] = size(control2_times{i});
    idx = setdiff(phase1:ncol, omit);   
    x = control2_times{i}(idx);
    y = control2_adjTrace{i}(:, idx) - control2_adjTrace{i}(:, end);
    y = y./y(:, 1);
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i+1})
end

for i=3:4    
    [~, ncol] = size(control3_times{i});
    idx = setdiff(phase1:ncol, omit);   
    x = control3_times{i}(idx);
    y = control3_adjTrace{i}(:, idx) - control3_adjTrace{i}(:, end);
    y = y./y(:, 1);
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i+3})
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
ylim([0 Inf])
legend([control2_labels(3:4), control3_labels(3:4)])

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=3:4  
    x = control2_times{i}(phase3:end);
    y = control2_adjTrace{i}(:, phase3:end) - control2_adjTrace{i}(:, end);
    y = y./y(:, 1);
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i+1})
end

for i=3:4    
    x = control3_times{i}(phase3:end);
    y = control3_adjTrace{i}(:, phase3:end) - control3_adjTrace{i}(:, end);
    y = y./y(:, 1);
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'Color', okabeIto{i+3})
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
ylim([0 Inf])
legend([control2_labels(3:4), control3_labels(3:4)])

% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
% for i=1:length(control3_basenames)
%     x = control3_times{i};
%     y = mean(control3_bgTrace{i}, 1, 'omitnan');
%     err = std(control3_bgTrace{i}, 0, 1, 'omitnan');
%     errorbar(x, y, err, 'Color', okabeIto{i})
%     xlabel('Time (minutes)')
%     ylabel('Fluorescence (A.U.)')
% end
% ylim([0 Inf])
% legend(control3_labels)
%% plot length over time
% control2_lcell_mean = cellfun(@(x)mean(x, 1, 'omitnan'), control2_lcell, 'UniformOutput', false);
% control2_lcell_std = cellfun(@(x)std(x, 0, 1, 'omitnan'), control2_lcell, 'UniformOutput', false);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(control2_basenames)
    subplot(3, 3, i)
    x = control2_times{i};
    y = mean(control2_lcell{i}, 1, 'omitnan');
    err = std(control2_lcell{i}, 0, 1, 'omitnan');
    errorbar(x, y, err, '-', 'Color', okabeIto{i})
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    ylim([0 Inf])
    title(control2_labels{i})
end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(control2_basenames)
    subplot(3, 3, i)
    x = control2_times{i};
    dy = diff(control2_lcell{i}, 1, 2);
    y_dot = dy./control2_lcell{i}(:, 2:end);
    y = [0 mean(y_dot, 1, 'omitnan')];
    err = [0 std(y_dot, 0, 1, 'omitnan')];
    errorbar(x, y, err, '-', 'Color', okabeIto{i})
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    %ylim([0 Inf])
    title(control2_labels{i})
end

%%
%since these controls are all taken in one continuous movie, the first 8-10
%minutes can be ignored when subtracting the background value and
%normalizing. How necessary is it to normlize to the pre-lysis vs
%post-lysis value? How much of a difference does it make (mathematically)?
cutoff = 4; %point in the movie at which to start the analyze
omit = [5:7];

control2_sliceTime = dataSlice(control2_times, cutoff);
control2_sliceIntensity = dataSlice(control2_cellTrace, cutoff, omit);
control2_sliceBackground = dataSlice(control2_bgTrace, cutoff, omit);
control2_sliceAdjintensity = dataSlice(control2_adjTrace, cutoff, omit);

control3_sliceTime = dataSlice(control3_times, cutoff);
control3_sliceIntensity = dataSlice(control3_cellTrace, cutoff, omit);
control3_sliceBackground = dataSlice(control3_bgTrace, cutoff, omit);
control3_sliceAdjintensity = dataSlice(control3_adjTrace, cutoff, omit);

%% plot the 'sliced' traces
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(control2_basenames)
    subplot(3, 3, i)
    x = control2_sliceTime{i};
    y = control2_sliceAdjintensity{i};
    plot(x, y, 'Color', okabeIto{i})
%     [idx, ~] = size(control2_sliceAdjintensity{i});
%     for j=1:idx
%         scatter(x, y(j, :), 15, okabeIto{i}), hold on
%     end
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    ylim([0 Inf])
    title(control2_labels{i})
end

cd(dirsave)
saveas(gcf, 'figure06.png')
saveas(gcf, 'figure06.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(control3_basenames)
    subplot(1, 5, i)
    x = control3_sliceTime{i};
    y = control3_sliceAdjintensity{i};
    plot(x, y, 'Color', okabeIto{i})
%     [idx, ~] = size(control3_sliceAdjintensity{i});
%     for j=1:idx
%         scatter(x, y(j, :), 15, okabeIto{i}), hold on
%     end
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    ylim([0 Inf])
    title(control3_labels{i})
end

cd(dirsave)
saveas(gcf, 'figure07.png')
saveas(gcf, 'figure07.fig')


%% determine alpha and rho
alpha0=1;
rho0=1;

[control2_alpha, control2_alphaYhat] = alphaCalc(control2_sliceTime, control2_sliceAdjintensity, alpha0);
[control3_alpha, control3_alphaYhat] = alphaCalc(control3_sliceTime, control3_sliceAdjintensity, alpha0);

[control2_rho, control2_rhoYhat, control2_frames] = rhoCalc(control2_sliceTime, control2_sliceAdjintensity, rho0);
[control3_rho, control3_rhoYhat, control3_frames] = rhoCalc(control3_sliceTime, control3_sliceAdjintensity, rho0);

%% what do the sliced traces look like?
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(control2_basenames)
    subplot(3, 3, i)
    plot(control2_times{i}, control2_adjTrace{i}, 'Color', okabeIto{i}), hold on
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    ylim([0 Inf])
    title(control2_labels{i})
end

cd(dirsave)
saveas(gcf, 'figure04.png')
saveas(gcf, 'figure04.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(control3_basenames)
    subplot(1, 5, i)
    plot(control3_times{i}, control3_adjTrace{i}, 'Color', okabeIto{i}), hold on
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    ylim([0 Inf])
    title(control3_labels{i})
end

cd(dirsave)
saveas(gcf, 'figure05.png')
saveas(gcf, 'figure05.fig')

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

function [slicedData] = dataSlice(data, cutoff, omit)
    
    slicedData = cell(size(data));
    
    if nargin < 3
        for i=1:height(data)
            dataTemp = data{i};
            slicedData{i} = dataTemp(:, cutoff:end)-dataTemp(:, cutoff);
        end       
    else
        if isempty(omit)==0
            for i=1:height(data)
                dataTemp = data{i};
                [nrow, ncol] = size(dataTemp);
                x = setdiff(1:ncol, omit);
                for j=1:nrow
                    v = dataTemp(j, x);
                    vq = interp1(x, v, omit);
                    dataTemp(j, omit) = vq;
                end
                slicedData{i} = dataTemp(:, cutoff:end);
            end
        else
            for i=1:height(data)
                slicedData{i} = data{i}(:, cutoff:end);
            end
        end
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
        xq = lysis{i}(1,1):lysis{i}(1,2)-1;
        x = setdiff(1:ncol, xq);
        for j=1:nrow
            vq = interp1(tme(x), adjintensity(j, x), tme(xq))
            adjintensity(j, xq) = vq;
        end
        
        %subtract final fluor value
        adjintensity = adjintensity - adjintensity(:, end);
        
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

function [rho, rho_yhat] = rhoCalc(times, normTrace, rho0)
    
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
        frames = 1:ncol;
        
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
            
            %fit data to the function and evaluate
            for j=1:nrow
                A = normintensity(j,idx);
                modelfun = @(rho, x)A*exp(-x*rho);
                rhoTemp(j, 1)=nlinfit(xq-(idx-1), normintensity(j,xq), modelfun, rho0);
                yhatTemp(j, xq)=modelfun(rhoTemp(j), xq-(idx-1));
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

%to correct for photobleaching
function [correctedTrace]=photoCorrect(times, normTrace, beta, intercept)
        
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
       dC=@(C, beta, dt, b)(C/(beta*dt+b))*dt;
       %dC=@(C, alpha, dt, b)((C*exp(-dt/(alpha*dt+b)))*(1/(alpha*dt+b))*dt);
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                dCB(n,i) = dC(normintensity(n,i), beta, dt(i), intercept); %this is the amount of photobleaching that occured in our measured value
   
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
        correctedTrace{h, 5} = Cbl_exp;
        correctedTrace{h, 6} = unb_frac;
        
        end
        
end
