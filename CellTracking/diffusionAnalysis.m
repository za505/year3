%diffusionAnalysis.m
%Author: Zarina Akbary
%Date: 01 May 2022
%purpose: to analyze mNeonGreen.m diffusion out of the cell

clear, close all

%color palette
okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

%% User Input
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/aggregate/';
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/figures/CZI/';

basenames = {'04242022_Exp2', '04052022_Exp2', '04052022_Exp1'};
labels = {'Frame Rate = 3 s', 'Frame Rate = 10 s', 'Frame Rate = 20 s'};

%% Load the data and normalize
[cellTrace, bgTrace, adjTrace, times, lCell, frameRate] = loadData(dirpath, basenames);
[normTime, normTrace] = dataNormalize(times, adjTrace);

%% calculate alpha
% alpha0 = 1;
% 
% [alpha, alpha_yhat] = alphaCalc(normTime, normTrace, alpha0);
% [alphaCoef, alphaFit] = alphaLinear(frameRate, alpha);

%% calculate rho
param0 = [0, 0.01];

[rho, B, rho_yhat] = rhoCalc(normTrace,  param0);
[rhoCoef, rhoFit] = rhoLinear(frameRate, rho);

%% plot the fit of alpha
% p = 0;
% n1 = 1;
% n2 = length(basenames);

% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
% for i=1:length(basenames)
%     
%     p = p + 1;
%     subplot(n1, n2, p)
%     x = normTrace{i, 1};
%     y1 = normTrace{i, 2};
%     y2 = alpha_yhat{i};
%     
%     plot(x, y1, 'Color', okabeIto{i}), hold on
%     plot(x, y2, '--k')
%     ylim([0 Inf])    
%     xlabel('Time (minutes)')
%     ylabel('Fluorescence (A.U.)')
%     title([labels{i}])
% end

%% plot the fit of rho
p = 0;
n1 = 1;
n2 = length(basenames);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=1:length(basenames)
    
    p = p + 1;
    subplot(n1, n2, p)
    x = [0:length(normTime{i})-1];
    y1 = normTrace{i};
    y2 = rho_yhat{i};
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, '--k')
    ylim([0 Inf])    
    xlabel('Frame Number')
    ylabel('Normalized Fluorescence (A.U.)')
    title([labels{i}])
end

%% correct the control traces
rho1= 0.03;
[correctedTrace]=photoCorrect(normTime, normTrace, rho1);

%% plot the control traces
p = 0;
n1 = length(basenames);
n2 = 2; %3;

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(basenames)
    
%     p = p + 1;
%     subplot(n1, n2, p)
%     x = times{i};
%     y1 = cellTrace{i};
%     y2 = bgTrace{i};
%     
%     plot(x, y1, 'Color', okabeIto{i}), hold on
%     plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
%     ylim([0 Inf])    
%     xlabel('Time (minutes)')
%     ylabel('Fluorescence (A.U.)')
%     title(['Raw ' labels{i}])
    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i};
    y = normTrace{i};
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 Inf])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Normalized ' labels{i}])
    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i};
    y = correctedTrace{i,1};
    [nr, nc] = size(y);
    ysmooth = nan(size(y));
    for j=1:nr
        ysmooth(j, :) = movingaverage(y(j, :), 3);
    end
    
    plot(x, ysmooth, 'Color', okabeIto{i})
    ylim([0 1.5])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Corrected ' labels{i}])
end

%% plot rho vs frame rate
dx = [];
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=1:length(basenames)
    x = frameRate{i}(end);
    dx = [dx x];
    y = rho{i}';
    ybar = mean(y, 2, 'omitnan');
    err = std(y, 0, 2, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')
end  
plot([0, dx], rhoFit, '--k', 'LineWidth', 1.5)
xlabel('Frame Rate (minutes)')
ylabel('\rho')

%% plot B as function of frame rate
dx = [];
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(basenames)  
    x = frameRate{i}(end);
    dx = [dx x];
    y = B{i}';
    ybar = mean(y, 2, 'omitnan');
    err = std(y, 0, 2, 'omitnan');
    
    scatter(repelem(x, length(y)), y, 'MarkerFaceColor', okabeIto{i}, 'MarkerEdgeColor', okabeIto{i}), hold on
    scatter(x, ybar, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
    errorbar(x, ybar, err, 'Color', 'black')
end  
xlim([0 0.4])
xlabel('Frame Rate (minutes)')
ylabel('B')

%% the experiments
clear, close all

%color palette
okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/aggregate/';
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/figures/CZI/';

basenames = {'04242022_Exp1', '04252022_Exp1', '04262022_Exp1', '05052022_Exp2'};
labels = {'Untreated', 'PBS', 'Tunicamycin', 'Spent Media'};

rho1 = 0.03;
%imstart = [6; 11; 11];
%imstart = [7; 12; 12; 7];
%imstart = [5; 9; 10; 5];
imstart = [6; 11; 11; 5];
[cellTrace, bgTrace, adjTrace, times, lCell, frameRate] = loadData(dirpath, basenames);
[normTime, normTrace] = dataNormalize(times, adjTrace, imstart);
[correctedTrace]=photoCorrect(normTime, normTrace, rho1);

%% plot to compare lengths
p = 0;
n1 = 1;
n2 = 2; %length(basenames);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=[2,3] %1:length(basenames)
    
    p = p + 1;
    subplot(n1, n2, p)
    x = times{i};
    y = lCell{i,1};
    [nr, nc] = size(y);
    ysmooth = nan(size(y));
    for j=1:nr
        ysmooth(j, :) = movingaverage(y(j, :), 3);
    end
    
    plot(x, ysmooth, 'Color', okabeIto{i})  
    xlabel('Time (minutes)')
    ylabel('Length (\mum)')
    ylim([0 Inf])
    title([labels{i}])
end

%% calculate effective elongation rate
dt1=times{1,1}(2:end)-times{1,1}(1:end-1);
dl1=(lCell{1,1}(:,2:end)-lCell{1,1}(:,1:end-1))./((lCell{1,1}(:,1:end-1)+lCell{1,1}(:,2:end))/2);
v1=dl1./dt1;

dt2=times{4,1}(2:end)-times{4,1}(1:end-1);
dl2=(lCell{4,1}(:,2:end)-lCell{4,1}(:,1:end-1))./((lCell{4,1}(:,1:end-1)+lCell{4,1}(:,2:end))/2);
v2=dl2./dt2;

ef = cell(size(lCell));

for i = 1:length(basenames)
    dl = diff(lCell{i}, 1, 2);
    dt = diff(times{i}./60, 1, 2);
    ef{i} = (dl./dt) ./ lCell{i}(:, 1:end-1);
end

figure, hold on
for i = 1:length(basenames)
    mef = mean(ef{i}, 1, 'omitnan');
    sef = movingaverage(mef, 3);
    plot(times{i}(1:end-1), sef, 'Color', okabeIto{i})
      %plot(times{i}(1:20), ef{i}(:, 1:20), 'Color', okabeIto{i})
end
xline(7, '--k')
xline(12, '--k')

[nrow1, ~] = size(lCell{1});
[nrow2, ~] = size(lCell{2});
[nrow3, ~] = size(lCell{3});
[nrow4, ~] = size(lCell{4});

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
plot([repelem(1, nrow1),repelem(3, nrow2), repelem(4, nrow3)], [ef{1}(:, 4)', ef{2}(:, 4)', ef{3}(:, 4)'], 'o', 'MarkerFaceColor', okabeIto{1}, 'MarkerEdgeColor', okabeIto{8}, 'MarkerSize', 20)
plot(repelem(2, nrow4), ef{4}(:, 5)', 'o', 'MarkerFaceColor', okabeIto{4}, 'MarkerEdgeColor', okabeIto{8}, 'MarkerSize', 20)
plot(repelem(3, nrow2), ef{2}(:, 6)', 'o', 'MarkerFaceColor', okabeIto{2}, 'MarkerEdgeColor', okabeIto{8}, 'MarkerSize', 20)
plot(repelem(4, nrow3), ef{3}(:, 8)', 'o', 'MarkerFaceColor', okabeIto{3}, 'MarkerEdgeColor', okabeIto{8}, 'MarkerSize', 20)
xlim([0 5])
xticks([0 1 2 3 4 5])
xticklabels({' ', 'untreated', 'spent media', 'PBS', 'tunicamcyin', ' '})
ylabel('Effective Elongation Rate (h^{-1})')
xlabel('Treatment')
legend({'steady-state LB', 'steady-state spent media', 'PBS shift', 'tunicamycin shift'})

%% plot the experimental traces
cd(dirsave)

p = 0;
n1 = length(basenames);
n2 = 3;
lim = 60;

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:n1
    
    p = p + 1;
    subplot(n1, n2, p)
    x = times{i}(1:lim);
    y1 = cellTrace{i}(:, 1:lim);
    y2 = bgTrace{i}(:, 1:lim);
    
    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, 'LineStyle', '--', 'Color', okabeIto{i})
    ylim([0 Inf])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Raw ' labels{i}])
    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i}(1:lim);
    y = normTrace{i}(:, 1:lim);
    
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.2])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Normalized ' labels{i}])
    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i}(1:lim);
    y = correctedTrace{i,1}(:, 1:lim);
    [nr, nc] = size(y);
    ysmooth = nan(size(y));
    for j=1:nr
        ysmooth(j, :) = movingaverage(y(j, :), 3);
    end
    
    %plot(x, ysmooth, 'Color', okabeIto{i})
    plot(x, y, 'Color', okabeIto{i})
    ylim([0 1.5])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(['Corrected ' labels{i}])
end
% saveas(gcf, 'experimentalTraces.png')
% saveas(gcf, 'experimentalTraces.fig')
% saveas(gcf, 'untreatedOnly.png')
% saveas(gcf, 'untreatedOnly.fig')

%% Untreated vs Autolysis
x1 = normTime{1}(1:40);
y1 = correctedTrace{1,1}(:, 1:40);
ybar1 = mean(y1, 1, 'omitnan');
yerr1 = std(y1, 0, 1, 'omitnan');
ysmooth1 = movingaverage(ybar1, 3);

x2 = normTime{2}(1:40);
y2 = correctedTrace{2,1}(:, 1:40);
ybar2 = mean(y2, 1, 'omitnan');
yerr2 = std(y2, 0, 1, 'omitnan');
ysmooth2 = movingaverage(ybar2, 3);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
errorbar(x1, ysmooth1, yerr1, 'Color', okabeIto{1}, 'LineWidth', 1.5)
errorbar(x2, ysmooth2, yerr2, 'Color', okabeIto{2}, 'LineWidth', 1.5)
ylim([0 1.2])    
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
legend(labels(1:2))

%% Untreated vs Tunicamycin
x1 = normTime{1}(1:60);
y1 = correctedTrace{1,1}(:, 1:60);
ybar1 = mean(y1, 1, 'omitnan');
yerr1 = std(y1, 0, 1, 'omitnan');
ysmooth1 = movingaverage(ybar1, 3);

x3 = normTime{3}(1:60);
y3 = correctedTrace{3,1}(:, 1:60);
ybar3 = mean(y3, 1, 'omitnan');
yerr3 = std(y3, 0, 1, 'omitnan');
ysmooth3 = movingaverage(ybar3, 3);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
errorbar(x1, ysmooth1, yerr1, 'Color', okabeIto{1}, 'LineWidth', 1.5)
errorbar(x3, ysmooth3, yerr3, 'Color', okabeIto{3}, 'LineWidth', 1.5)
ylim([0 1.2])    
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
legend(labels([1,3]))

%% Untreated vs Spent Media
x1 = normTime{1}(1:60);
y1 = correctedTrace{1,1}(:, 1:60);
ybar1 = mean(y1, 1, 'omitnan');
yerr1 = std(y1, 0, 1, 'omitnan');
ysmooth1 = movingaverage(ybar1, 3);

x3 = normTime{4}(1:60);
y3 = correctedTrace{4,1}(:, 1:60);
ybar3 = mean(y3, 1, 'omitnan');
yerr3 = std(y3, 0, 1, 'omitnan');
ysmooth3 = movingaverage(ybar3, 3);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
errorbar(x1, ysmooth1, yerr1, 'Color', okabeIto{1}, 'LineWidth', 1.5)
errorbar(x3, ysmooth3, yerr3, 'Color', okabeIto{4}, 'LineWidth', 1.5)
ylim([0 1.2])    
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
legend(labels([1,4]))

%% calculate tau
tau0 = 1;
lim = 60;

[tau, tau_yhat] = tauCalc(normTime, correctedTrace, tau0, lim);

%% plot the fits to tau
cd(dirsave)

p = 0;
n1 = 1;
n2 = length(basenames);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=1:n2    
    p = p + 1;
    subplot(n1, n2, p)
    x = normTime{i}(1:lim);
    y1 = correctedTrace{i}(:, 1:lim);
    y2 = tau_yhat{i};

    plot(x, y1, 'Color', okabeIto{i}), hold on
    plot(x, y2, '--k')
    ylim([0 Inf])    
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    title([labels{i}])
end
saveas(gcf, 'tauFits.png')
saveas(gcf, 'tauFits.fig')

%% re-do correction
imstart2 = [5; 5; 5; 5];
[correctedTrace2, correctedTime]=fixNstitch(times, adjTrace, correctedTrace, rho1, imstart2);

%% re-plot corrected experimental traces
cd(dirsave)

Fmax = nan(height(correctedTrace2), 1);
foldingTrace = cell(height(correctedTrace2), 1);
flux = cell(height(correctedTrace2), 1);

j=0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=[1, 3, 4]
    i
    x = correctedTime{i}(1:40);
    y = correctedTrace2{i,1}(:, 1:40);
    if i==3
        idx = setdiff(1:height(y), [2,6]);
        y = y(idx, :);
    end
    ybar = mean(y, 1, 'omitnan');
    yerr = std(y, 0, 1, 'omitnan')/sqrt(height(correctedTrace2{i,1}));
    %ysmooth = movingaverage(ybar, 3);
    %errorbar(x, ysmooth, yerr, 'Color', okabeIto{i}, 'LineWidth', 1.5)
    Fmax(i) = max(ybar);
    modelfun = @(t, tau)(Fmax(i)-1)*(1-exp(-t./tau)) + 1;
    %tau = nlinfit(normTime{i}(1:40), correctedTrace{i}(1:40), modelfun, 3)
    if i==1
        tau=4.5;
    elseif i==2
        tau=5.5;
    else
        tau=3;
    end
        foldingTrace{i} = modelfun(normTime{i}, tau);
    
    ysub = (ybar(:, 5:40) - foldingTrace{i}(:, 1:36)) + 1;
    ycomb = [ybar(:, 1:4), ysub];
    errorbar(x, ycomb, yerr, 'Color', okabeIto{i}, 'LineWidth', 1.5)
    %errorbar(x, ybar, yerr, 'Color', okabeIto{i}, 'LineWidth', 1.5)
    %plot(correctedTime{i}(5:40), foldingTrace{i}(:, 1:36), 'Color', okabeIto{i}, 'LineStyle', '--')

end
ylim([0 1.2])    
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
legend(labels([1, 3, 4]))
% saveas(gcf, 'untreatedvsautolysis.png')
% saveas(gcf, 'untreatedvsautolysis.fig')

%% plot the derivative
cd(dirsave)

j=0;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=[1, 3, 4]
    i
    j=j+1;
    x = correctedTime{i}(1:40);
    y = correctedTrace2{i,1}(:, 1:40);
    if i==3
        idx = setdiff(1:height(y), [2,6]);
        y = y(idx, :);
    end
    dy = diff(y, 1, 2);
    dt = diff(x, 1, 2);
    df = dy./dt;
    ybar = mean(df(:, 16:end), 2, 'omitnan');
    yerr = std(ybar, 1, 'omitnan'); %/sqrt(height(correctedTrace2{i,1}));
    ybar = mean(ybar, 'omitnan');
    
    bar(j, abs(ybar), 'FaceColor', okabeIto{i}, 'EdgeColor', okabeIto{i})
    errorbar(j, abs(ybar), yerr, 'Color', okabeIto{i}, 'LineStyle', 'none')
    %errorbar(x(1:end-1), ybar, yerr, 'Color', okabeIto{i}, 'LineWidth', 1.5)
    %plot(x(1:end-1), y,'Color', okabeIto{i}, 'LineWidth', 1.5)
end  
xlabel('Treatment')
ylabel('Fluorescence Derivative')
legend(labels([1, 3, 4]))

%% calculate max flux 
[maxFlux, Fmax, foldingTrace] = fluxCalc(normTime([1, 3, 4],1), correctedTrace([1, 3, 4],1));

%% plot the max flux
cd(dirsave)

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=1:length(basenames)
    x = i;
    y = maxFlux{i,1}(:, 1);
    if i==3
       idx = setdiff(1:height(y), [2,6]);
       y = y(idx, :);
    end
    ybar = abs(mean(y, 1, 'omitnan'));
    yerr = std(y, 0, 1, 'omitnan');
    bar(x, ybar, 'FaceColor', okabeIto{i}, 'EdgeColor', okabeIto{i})
    errorbar(x, ybar, yerr, 'Color', okabeIto{i}, 'LineWidth', 1.5, 'LineStyle', 'none')
end 
xlabel('Treatment')
xticks(1:4)
xticklabels(labels)
ylabel('Maximum Flux')
legend({'Untreated', '', 'PBS', '', 'Tunicamycin', '', 'Spent Media', ''})

%% re-plot untreated vs spent
cd(dirsave)
x1 = correctedTime{1}(1:60);
y1 = correctedTrace2{1,1}(:, 1:60);
ybar1 = mean(y1, 1, 'omitnan');
yerr1 = std(y1, 0, 1, 'omitnan');
ysmooth1 = movingaverage(ybar1, 3);

x4 = correctedTime{4}(1:60);
y4 = correctedTrace2{4,1}(:,1:60);
ybar4 = mean(y4, 1, 'omitnan');
yerr4 = std(y4, 0, 1, 'omitnan');
ysmooth4 = movingaverage(ybar4, 3);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
errorbar(x1, ysmooth1, yerr1, 'Color', okabeIto{1}, 'LineWidth', 1.5)
errorbar(x4, ysmooth4, yerr4, 'Color', okabeIto{5}, 'LineWidth', 1.5)
ylim([0 1.2])    
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
legend(labels([1,4]))
saveas(gcf, 'untreatedvsspent.png')
saveas(gcf, 'untreatedvsspent.fig')

%% plot the distribution of tau
cd(dirsave)

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=1:length(basenames)
    x = [0:80];
    y = pdf('Normal', x, mean(tau{i}, 'omitnan'), std(tau{i}, 'omitnan'));
    %histogram(tau{i}, 'Normalization', 'pdf', 'BinWidth', 1, 'EdgeColor', okabeIto{i}, 'FaceColor', okabeIto{i})
    plot(x, y, 'Color', okabeIto{i}, 'LineWidth', 1.5)
end
legend(labels)
xlabel(['\tau (minutes)'])
ylabel('Probability')
saveas(gcf, 'tauPDF.png')
saveas(gcf, 'tauPDF.fig')

%% plot corrected experimental traces
p = 0;
n1 = 2;
n2 = 2;

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 30), hold on
for i=1:length(basenames)

    p = p + 1;
    subplot(n1, n2, p)
    x = correctedTime{i}(1:60);
    y = correctedTrace2{i,1}(:, 1:60);
    [nr, nc] = size(y);
    ysmooth = nan(size(y));
    for j=1:nr
        ysmooth(j, :) = movingaverage(y(j, :), 3);
    end
    
    plot(x, ysmooth, 'Color', okabeIto{i})
    ylim([0 1.1])    
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    title(labels{i})
end

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
        frameRate{i} = diff(time, 1, 2);
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

function [normTime, normTrace] = dataNormalize(times, adjTrace, imstart)
    
    %pre-allocate variables
    [nrow, ~] = size(adjTrace);
    normTime = cell(nrow, 1);
    normTrace = cell(nrow, 1);
    
    for i=1:height(adjTrace)  
        
        adjintensity = adjTrace{i};
        tme = times{i};
        
        dt = diff(tme, 1, 2);
        
        [nrow, ncol] = size(adjintensity);
        
        if nargin < 3
            imstart = find(dt == dt(end), 1, 'first');
            prelysis = imstart;
        else
            prelysis = imstart(i);
            lval = [prelysis + 1 : prelysis + 4];
            nval = setdiff(1:ncol, lval); %non-lysis frames
            pval = setdiff(prelysis:ncol, lval);
            fval = setdiff(prelysis:ncol-1, lval);
        
            xq = tme(lval); %the lysis times
            x = setdiff(tme, xq); %the non-lysis times
        
            for j=1:nrow
                vq = interp1(x, adjintensity(j, nval), xq);
                adjintensity(j, lval) = vq;
            end
        end
        
        
        %adjust the time vector
        normTime{i,1} = tme(prelysis:end) - tme(prelysis);
        
        %normalize the trace by the pre-lysis frame
        normTrace{i, 1} = adjintensity(:, prelysis:end) ./ adjintensity(:, prelysis);

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
                modelfun = @(alpha, t)A*exp(-t./alpha) + (1-A);
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

function [rho, B, rho_yhat] = rhoCalc(normTrace, param0)
    
    %pre-allocate variables
    rho = cell(height(normTrace), 1);
    B = cell(height(normTrace), 1);
    rho_yhat = cell(height(normTrace), 1);
    
    for i=1:height(normTrace)
        
        normintensity = normTrace{i};
        
        %convert negative values to NaN
        normintensity(normintensity<0) = NaN;
    
        %index traces that have at least 1 data point
        [nrow, ncol]=size(normintensity); 
        idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
        normintensity=normintensity(idx, :);
        frames = 0:ncol-1;
        
        %pre-allocate output variables
        [nrow, ~]=size(normintensity);
        rhoTemp=nan(nrow,2);    
        yhatTemp=nan(size(normintensity));
    
        %find the last time point of interest
%         di = diff(normintensity, 2, 2);
%         mdi = abs(mean(di, 1, 'omitnan'));
%         [~, fidx] = min(mdi);
        
        mint = mean(normintensity, 1, 'omitnan');
        fidx = find(mint >= 0.15, 1, 'last');
        
        %define exponential function    
        modelfun = @(param, x)(1-param(2)) * exp(-x*param(1)) + param(2);
        x = frames(1:fidx);
        
        %fit data to the function and evaluate
        for j=1:nrow
            y = normintensity(j,1:fidx);
            paramTemp(j, :)=nlinfit(x, y, modelfun, param0);
            yhatTemp(j, 1:fidx)=modelfun(paramTemp(j, :), x);
            yhatTemp(j, fidx+1:end)=NaN;
        end
                          
        rho{i} = paramTemp(:, 1);
        B{i} = paramTemp(:, 2);
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

function [correctedTrace]=photoCorrect(normTime, normTrace, rho)
        
         correctedTrace = cell(height(normTrace), 6);
        
        for h=1:height(normTrace)
            
            normintensity = normTrace{h};
            tme = normTime{h};
            
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
        
        %this formula comes from the slope and intercept calculated for the 1.2, 2,
        %and 3 second tau vs frame rate controls 
       %dC=@(C, beta, dt, b)(C/(beta*dt+b))*dt;
       %dC=@(C, alpha, dt, b)((C*exp(-dt/(alpha*dt+b)))*(1/(alpha*dt+b))*dt);
        %dC=@(C, beta, t, dt, b)(C/(beta*t+b))*dt;
       %dC=@(C, omega, dt, b)(C*(phi*dt+b))*dt;
       
        %calculate dt (the dt between frames may vary in a single run).
        %Note: this is also the frame rate
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                %dCB(n,i) = dC(normintensity(n,i), beta, tme(i), dt(i), intercept); %this is the amount of photobleaching that occured in our measured value
   
                %dCB(n,i) = dC(normintensity(n,i), phi, dt(i), intercept);
                %rho = phi + intercept;
                
                dCB(n,i) = normintensity(n,i) * rho;
                
                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                %dCP(n,i)=dCP(n,i)*unb_frac(n,i); %Correcting for the fact that a fraction of fluorophores are unbleached

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

function [correctedTrace2, correctedTime]=fixNstitch(times, adjTrace, correctedTrace, rho, imstart)
    
    correctedTrace2 = cell(size(correctedTrace));
    correctedTime = cell(size(times));
    
    for i = 1:height(adjTrace)
        prelysis = 1:imstart(i)-1;
        adjintensity = adjTrace{i};
        tme = times{i}(prelysis);
        
        [nrow, ~] = size(adjintensity);
        normintensity = adjintensity(:, prelysis) ./ adjintensity(:, 1);
        %[~, ncol2] = size(correctedTrace{i});
        
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:,1);
               
        %the correction
        for n=1:nrow
           
            for t=1:length(tme)-1
                
                dCB = normintensity(n,t) * rho;
                
                dCT = normintensity(n, t+1) - normintensity(n, t); %this is the total fluor. loss for the measured value

                dCP = dCT + dCB; %this is the amount of loss attributable to permeability

                Cnew(n,t+1)=Cnew(n,t)+dCP;

            end  
        end
        
%         %interpolate values where LB is not perfused/where a time lapse
%         %occured
%         tcol = length(prelysis) + ncol2;
%         
%         if tcol~=ncol
%             postlysis = ncol-ncol2;
%             xq = setdiff(1:ncol, prelysis);
%             xq = setdiff(xq, postlysis+1:ncol);
%             x = [tme(prelysis), times{i}(postlysis+1)];
%             x = tme(prelysis);
%             vq = nan(nrow, length(xq));
%             for n = 1:nrow
%                 y = [normintensity(n, prelysis), correctedTrace{i}(n, 1)];
%                 y = normintensity(n, prelysis);
%                 vq(n, :) = interp1(x, y, xq, 'linear', 'extrap');
%             end
%             correctedTrace2{i} = [Cnew, vq, correctedTrace{i}];
%         else
            correctedTrace2{i} = [Cnew, correctedTrace{i}];
            [~, ncol] = size(correctedTrace2{i});
            correctedTime{i} = times{i}(1:ncol);
%         end
        
    end
end

function [tau, tau_yhat] = tauCalc(times, correctedTrace, tau0, lim)
    
    tau = cell(height(correctedTrace), 1);
    tau_yhat = cell(height(correctedTrace), 1);
    
    for i = 1:height(correctedTrace)
        if nargin > 3
            x = times{i}(1:lim);
            y = correctedTrace{i}(:, 1:lim);
        else
            x = times{i};
            y = correctedTrace{i};
        end
        [nrow, ncol] = size(y);
        tauTemp = nan(nrow, 1);
        yhatTemp = nan(size(y));
        
%         ybar = mean(y, 1, 'omitnan');
%         A = ybar(1);
%         B = ybar(end);
        
        %modelfun = @(tau, t) A * exp(-t./tau) + B;
        
        for j = 1:nrow
            %aidx = find(~isnan(y(j, :)), 1, 'first');
            %A = y(j, aidx);
            %bidx = find(~isnan(y(j, :)), 1, 'last');
            [~, aidx] = max(y(j, :));
            A = y(j, aidx);
            [~, bidx] = min(y(j, :));
            B = y(j, bidx);
            modelfun = @(tau, t) (A-B) * exp(-t./tau) + B;
            %modelfun = @(tau, t) A * exp(log(A/B)*(1-exp(t./tau)));
            %modelfun = @(tau, t) exp(-t./tau);
            tauTemp(j) = nlinfit(x(aidx:bidx)-x(aidx), y(j, aidx:bidx), modelfun, tau0);
            yhatTemp(j, aidx:bidx) = modelfun(tauTemp(j, :), x(aidx:bidx)-x(aidx));
            if aidx~=1
                yhatTemp(j, 1:aidx-1) = NaN;
            end
        end
        
        tau{i} = tauTemp;
        tau_yhat{i} = yhatTemp;
    end
    
end

function [maxFlux, Fmax, foldingTrace] = fluxCalc(times, correctedTrace)
    
    %pre-allocate variables
    maxFlux = cell(height(correctedTrace), 1);
    Fmax = nan(height(correctedTrace), 1);
    foldingTrace = cell(height(correctedTrace), 1);

    for i=1:height(correctedTrace)
        t = times{i};
        y = correctedTrace{i};
        if i==2
            idx = setdiff(height(y), [2,6]);
            y = correctedTrace{i}(idx, :);
        end
        ybar = mean(y, 1, 'omitnan');
        ybar = movingaverage(ybar, 3);
        %figure,plot(ybar),pause
        Fmax(i) = max(ybar);
        y = movingaverage(y, 3);
        dy = diff(y, 1, 2);
        [minv, mini] = min(dy, [], 2, 'omitnan');

        maxFlux{i, 1} = [minv, mini];
        
        modelfun = @(t, tau)(Fmax(i)-1)*(1-exp(-t./tau)) + 1;
        %tau = nlinfit(t, ybar, modelfun, 3);
        foldingTrace{i} = modelfun(t, 3);

    end  
    
end
