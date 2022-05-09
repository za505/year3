%Author: Zarina Akbary
%Date: 04/07/2022
%Purpose: To aggregate mNeonGreen diffusion data and save figures for my
%04/07/2022 meeting with Rico

clear, close all

%% User Input
okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

color_gray = {[204 204 204], [150, 150, 150], [82, 82, 82]};
color_gray = cellfun(@(x)x./255, color_gray, 'UniformOutput', false);

transparency = 0.3; %this is the alpha argument for the ciplot function

%location and names of of the *.mat files 
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis';
figsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/04072022_analysis/';
basenames = {'10232021_Exp1', '10262021_Exp1', '04042022_Exp1', '04052022_Exp3', '02212022_Exp2', '02122022_Exp1','02122022_Exp2', '02092022_Exp1', '02212022_Exp1', '11192021_Exp2', '12082021_Exp3', '11192021_Exp1', '11302021_Exp1', '10232021_Exp2', '10262021_Exp2', '01142022_Exp1', '01172022_Exp1', '01172022_Exp2', '01242022_Exp1', '01262022_Exp1', '02122022_Exp3', '02192022_Exp1', '03012022_Exp1', '04112022_Exp3', '11202021_Exp1', '12082021_Exp1', '12082021_Exp2', '04052022_Exp1', '04052022_Exp2', '02192022_Exp2', '02192022_Exp3', '02192022_Exp4'};

%% Load the control data
basenameControls = {'11202021_Exp1', '12082021_Exp1', '12082021_Exp2', '04052022_Exp2', '04052022_Exp1'};
variNames = {'time', 'tme', 'imstart', 'lcell', 'bgintensity', 'intensity', 'adjintensity', 'normintensity', 'Cnew', 'alpha', 'alpha_yhat', 'rho', 'rho_yhat'};
alpha0 = 1;
rho0 = 0;

%pre-allocate cell array
controls = cell(length(basenameControls)+1, length(variNames));
controls(1, :) = variNames; 

for i=1:length(basenameControls)
    
    cd([dirsave '/normalizedFiles'])
    
    basename=basenameControls{i}; 
    load([basename '_norm.mat'])
    
    controls{i+1, 1} = time;
    controls{i+1, 2} = tme;
    controls{i+1, 3} = imstart;
    controls{i+1, 4} = lcell;    
    controls{i+1, 5} = bgintensity;
    controls{i+1, 6} = intensity;    
    controls{i+1, 7} = adjintensity;
    controls{i+1, 8} = normintensity;        
    
    [controls{i+1, 10}, controls{i+1, 11}] = alphaCalc(tme, normintensity, alpha0);
    [controls{i+1, 12}, controls{i+1, 13}] = rhoCalc(tme, normintensity, rho0);
    
    cd([dirsave '/correctedFiles'])
    load([basename '_corrected.mat'])
    
    controls{i+1, 9} = Cnew;
        
end

%% plot untreated: corrected intensity vs tme
labels = {'1 min untreated', 'media switch'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=[3,8] %1:height(untreated)-1
p1 = meanPlot(untreated{i+1,2}, untreated{i+1,9}, okabeIto{i}, transparency);
end
%p2 = meanPlot(untreated{3,2}, untreated{3,9}, okabeIto{2}, transparency);
%p3 = meanPlot(untreated{4,2}, untreated{4,9}, okabeIto{3}, transparency);
%p4 = meanPlot(untreated{5,2}, untreated{5,9}, okabeIto{4}, transparency);
%p5 = meanPlot(untreated{6,2}, untreated{6,9}, okabeIto{5}, transparency);
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')
legend(labels)
%lgd = legend([p1(1), p2(1), p3(1), p4(1), p5(1)], labels)
%title(lgd,'Frame Rate')
ylim([0 Inf])

%% load untreated (longer frame rate) data
basenameUntreated = {'10232021_Exp1', '10262021_Exp1', '04042022_Exp1', '04052022_Exp3','02122022_Exp1', '02122022_Exp2', '02092022_Exp1', '04112022_Exp3'};
variNames = {'time', 'tme', 'imstart', 'lcell', 'bgintensity', 'intensity', 'adjintensity', 'normintensity', 'Cnew', 'alpha', 'alpha_yhat', 'rho', 'rho_yhat'};
alpha0 = 1;
rho0 = 0;

%pre-allocate cell array
untreated = cell(length(basenameUntreated)+1, length(variNames));
untreated(1, :) = variNames; 

for i=1:length(basenameUntreated)
    
    cd([dirsave '/normalizedFiles'])
    
    basename=basenameUntreated{i}; 
    load([basename '_norm.mat'])
    
    untreated{i+1, 1} = time;
    untreated{i+1, 2} = tme;
    untreated{i+1, 3} = imstart;
    untreated{i+1, 4} = lcell;    
    untreated{i+1, 5} = bgintensity;
    untreated{i+1, 6} = intensity;    
    untreated{i+1, 7} = adjintensity;
    untreated{i+1, 8} = normintensity;        
    
    [untreated{i+1, 10}, untreated{i+1, 11}] = alphaCalc(tme, normintensity, alpha0);
    [untreated{i+1, 12}, untreated{i+1, 13}] = rhoCalc(tme, normintensity, rho0);
    
    cd([dirsave '/correctedFiles'])
    load([basename '_corrected.mat'])
    
    untreated{i+1, 9} = Cnew;
        
end

%% plot B as a function of frame rate
B = [controls{2, 5}(end), mean([untreated{2, 5}(end), untreated{3, 5}(end), untreated{4, 5}(end)]), untreated{5,5}(end), untreated{6,5}(end), untreated{7,5}(end), untreated{8,5}(end)];
dt = [1.2/60, 1, 2, 5, 10, 20]; 

%calculate the final background value, B, as a function of frame rate
% linearCoef = polyfit(dt, B, 1);
% linearFit = polyval(linearCoef,[0 dt]);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 20), hold on
scatter(dt, B, 50, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
%plot([0, dt], linearFit, '--k', 'LineWidth', 1.5)

% caption1 = sprintf('%f * frame rate + %f', linearCoef(1), linearCoef(2));
% text(dt(3), B(3)-B(3)/10, ['B = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

xlim([-.2, dt(end)+dt(end)/10])
ylabel('B (AU)')
xlabel('Frame Rate (minutes)')

%Now exclude the 10 and 20-minute time points
B = [controls{2, 5}(end), mean([untreated{2, 5}(end), untreated{3, 5}(end), untreated{4, 5}(end)]), untreated{5,5}(end), untreated{6,5}(end)];
dt = [1.2/60, 1, 2, 5]; 

%calculate the final background value, B, as a function of frame rate
linearCoef = polyfit(dt, B, 1);
linearFit = polyval(linearCoef,[0 dt]);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 20), hold on
scatter(dt, B, 50, 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black')
plot([0, dt], linearFit, '--k', 'LineWidth', 1.5)

caption1 = sprintf('%f * frame rate + %f', linearCoef(1), linearCoef(2));
text(dt(3), B(3)-B(3)/10, ['B = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

xlim([-.2, dt(end)+dt(end)/10])
ylabel('B (AU)')
xlabel('Frame Rate (minutes)')

%% plot controls: intensity & background vs time
labels = {'Cellular Trace', 'Background Trace'}
titles = {'Frame Rate = 1 second', 'Frame Rate = 2 seconds', 'Frame Rate = 3 seconds', 'Frame Rate = 10 seconds', 'Frame Rate = 20 seconds'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15)
for i=1:length(titles)
    x = controls{i+1, 1}(:, controls{i+1, 3}:end) - controls{i+1, 1}(:, controls{i+1, 3});
    y1 = controls{i+1, 6}(:, controls{i+1, 3}:end);
    y2 = controls{i+1, 5}(:, controls{i+1, 3}:end);
    
    subplot(1, 5, i), hold on
    p1 = plot(x, y1, 'Color', okabeIto{i});
    p2 = plot(x, y2, '--', 'Color', okabeIto{i}, 'LineWidth', 1.5);
    ylabel('Fluorescence (A.U.)')
    xlabel('Time (minutes)')
    title(titles{i})
    legend([p1(1), p2], labels)
end

%% Simulate averaged control traces
labels = {'Trace', 'Simulation'}
titles = {'Frame Rate = 1 second', 'Frame Rate = 2 seconds', 'Frame Rate = 3 seconds', 'Frame Rate = 10 seconds', 'Frame Rate = 20 seconds'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15)
for i=1:length(titles)
    y = mean(controls{i+1, 6}(:, controls{i+1, 3}:end), 1, 'omitnan');
    time = controls{i+1, 1}(:, controls{i+1, 3}:end) - controls{i+1, 1}(:, controls{i+1, 3});
    tau = 0;
    A = y(1);
    B = polyval(linearCoef, time(end)-time(end-1));
    beta = 33.5358;
    b = -0.0204;

    [Csim] = photobleachSim(time, tau, A, beta, b)
    
    subplot(1, 5, i), hold on
    p1 = plot(time, y, 'Color', okabeIto{i}, 'LineWidth', 1);
    p2 = plot(time, Csim, '--', 'Color', okabeIto{8}, 'LineWidth', 1);
    ylabel('Fluorescence (A.U.)')
    xlabel('Time (minutes)')
    title(titles{i})
    legend([p1(1), p2(1)], labels)
end

%% Simulate averaged control traces
labels = {'Trace', 'Simulation'}
titles = {'Frame Rate = 1 second', 'Frame Rate = 2 seconds', 'Frame Rate = 3 seconds', 'Frame Rate = 10 seconds', 'Frame Rate = 20 seconds'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15)
for i=1:length(titles)
    y = mean(controls{i+1, 6}(:, controls{i+1, 3}:end), 1, 'omitnan');
    time = controls{i+1, 1}(:, controls{i+1, 3}:end) - controls{i+1, 1}(:, controls{i+1, 3});
    tau = 0;
    A = y(1);
    B = polyval(linearCoef, time(end)-time(end-1));
    beta = 50;
    b = -0.0204;

    [Csim] = photobleachSim(time, tau, A, beta, b);
    
    subplot(1, 5, i), hold on
    p1 = plot(time, y, 'Color', okabeIto{i}, 'LineWidth', 1);
    p2 = plot(time, Csim, '--', 'Color', okabeIto{8}, 'LineWidth', 1);
    ylabel('Fluorescence (A.U.)')
    xlabel('Time (minutes)')
    title(titles{i})
    legend([p1(1), p2(1)], labels)
end

%% plot controls: normalized intensity vs tme
labels = {'1 second', '2 seconds', '3 seconds', '10 seconds', '20 seconds'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
p1 = meanPlot(controls{2,2}, controls{2,8}, okabeIto{1}, transparency);
p2 = meanPlot(controls{3,2}, controls{3,8}, okabeIto{2}, transparency);
p3 = meanPlot(controls{4,2}, controls{4,8}, okabeIto{3}, transparency);
p4 = meanPlot(controls{5,2}, controls{5,8}, okabeIto{4}, transparency);
p5 = meanPlot(controls{6,2}, controls{6,8}, okabeIto{5}, transparency);
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
lgd = legend([p1(1), p2(1), p3(1), p4(1), p5(1)], labels)
title(lgd,'Frame Rate')
ylim([0 Inf])

%% show the exponential fit for alpha = beta * f + b
labels = {'Trace', 'Fit'}
titles = {'Frame Rate = 1 second', 'Frame Rate = 2 seconds', 'Frame Rate = 3 seconds', 'Frame Rate = 10 seconds', 'Frame Rate = 20 seconds'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15)
for i=1:length(titles)
    subplot(1, 5, i), hold on
    p1 = plot(controls{i+1, 2}, controls{i+1, 8}, 'Color', okabeIto{i}, 'LineWidth', 1);
    p2 = plot(controls{i+1, 2}, controls{i+1, 11}, '--', 'Color', okabeIto{8}, 'LineWidth', 1);
    ylabel('Normalized Fluorescence (A.U.)')
    xlabel('Time (minutes)')
    ylim([0 Inf])
    title(titles{i})
    legend([p1(1), p2(1)], labels)
end

%% calculate the slope, beta, from alpha and frame rate
alpha = {controls{2, 10}; controls{3, 10}; controls{4, 10}};
dt = [1.2/60, 2/60, 3/60]; 

linearCoef1 = polyfit([repelem(dt(1), length(alpha{1})), repelem(dt(2), length(alpha{2})), repelem(dt(3), length(alpha{3}))],[alpha{1}', alpha{2}', alpha{3}'], 1);
linearFit1= polyval(linearCoef1,[0 dt]);

%% plot alpha as a function of frame rate
alpha_means = cellfun(@(x)mean(x, 1, 'omitnan'), alpha);
alpha_std = cellfun(@(x)std(x, 0, 1, 'omitnan'), alpha);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 20), hold on
p1 = scatter(repelem(dt(1), length(alpha{1})), alpha{1}, 50, 'MarkerFaceColor', '#969696', 'MarkerEdgeColor', '#969696')
scatter(repelem(dt(2), length(alpha{2})), alpha{2}, 50, 'MarkerFaceColor', '#969696', 'MarkerEdgeColor', '#969696')
scatter(repelem(dt(3), length(alpha{3})), alpha{3}, 50, 'MarkerFaceColor', '#969696', 'MarkerEdgeColor', '#969696')

p2 = scatter(dt, alpha_means', 50, 'MarkerFaceColor', 'black')
errorbar(dt(1), alpha_means(1), alpha_std(1), 'Color', 'black', 'LineWidth', 1.5)
errorbar(dt(2), alpha_means(2), alpha_std(2), 'Color', 'black', 'LineWidth', 1.5)
errorbar(dt(3), alpha_means(3), alpha_std(3), 'Color', 'black', 'LineWidth', 1.5)

p3 = plot([0, dt], linearFit1, '--k', 'LineWidth', 1.5)

xlim([0, dt(end)+dt(end)/5])
ylabel('\alpha')
xlabel('Frame Rate (minutes)')

caption1 = sprintf('%f * frame rate + %f', linearCoef1(1), linearCoef1(2));
text(dt(2), alpha_means(2)-alpha_means(2)/4, ['\alpha = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

%% plot rho as a function of frame rate
rho = {controls{2, 12}; controls{3, 12}; controls{4, 12}; controls{5, 12}; ; controls{6, 12}; untreated{4, 12}; untreated{5, 12}};
dt = [1.2/60, 2/60, 3/60, 10/60, 20/60, 1, 2]; 

rho_means = cellfun(@(x)mean(x, 1, 'omitnan'), rho);
rho_std = cellfun(@(x)std(x, 0, 1, 'omitnan'), rho);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 20), hold on
for i=1:length(dt)
p1 = scatter(repelem(dt(i), length(rho{i})), rho{i}, 50, 'MarkerFaceColor', '#969696', 'MarkerEdgeColor', '#969696')
errorbar(dt(i), rho_means(i), rho_std(i), 'Color', 'black', 'LineWidth', 1.5)
end
scatter(dt, rho_means', 50, 'MarkerFaceColor', 'black')

%p3 = plot([0, dt], linearFit1, '--k', 'LineWidth', 1.5)

xlim([0, dt(end)+dt(end)/5])
ylabel('\rho')
xlabel('Frame Rate (minutes)')

%caption1 = sprintf('%f * frame rate + %f', linearCoef1(1), linearCoef1(2));
%text(dt(2), rho_means(2)-rho_means(2)/4, ['\rho = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

%% plot controls: corrections vs tme
labels = {'Normalized Trace', 'Corrected Trace'}
titles = {'Frame Rate = 1 second', 'Frame Rate = 2 seconds', 'Frame Rate = 3 seconds', 'Frame Rate = 10 seconds', 'Frame Rate = 20 seconds'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15)
for i=1:length(titles)
    subplot(1, 5, i), hold on
    p1 = plot(controls{i+1, 2}, controls{i+1, 8}, 'Color', okabeIto{i});
    p2 = plot(controls{i+1, 2}, controls{i+1, 9}, '--', 'Color', okabeIto{8}, 'LineWidth', 1.5);
    ylabel('Normalized Fluorescence (A.U.)')
    xlabel('Time (minutes)')
    %ylim([0 Inf])
    title(titles{i})
    legend([p1(1), p2(1)], labels)
end

%% simulate untreated traces
labels = {'Trace', 'Simulation'}
titles = {'Frame Rate = 1 minute', 'Frame Rate = 1 minute', 'Frame Rate = 2 minutes', 'Frame Rate = 5 minutes', 'Frame Rate = 10 minutes', 'Frame Rate = 20 minutes'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15)
for i=1:length(titles)
    y = mean(untreated{i+1, 6}(:, untreated{i+1, 3}:end), 1, 'omitnan');
    time = untreated{i+1, 1}(:, untreated{i+1, 3}:end) - untreated{i+1, 1}(:, untreated{i+1, 3});
    tau = 100;
    A = y(1);
    B = polyval(linearCoef, time(end)-time(end-1));
    %beta = 33.5358;
    beta = 50;
    b = -0.0204;

    [Csim] = photobleachSim(time, tau, A, beta, b);
    
    subplot(2, 3, i), hold on
    p1 = plot(time, y, 'Color', okabeIto{i}, 'LineWidth', 1);
    p2 = plot(time, Csim, '--', 'Color', okabeIto{8}, 'LineWidth', 1);
    ylabel('Fluorescence (A.U.)')
    xlabel('Time (minutes)')
    title(titles{i})
    legend([p1(1), p2(1)], labels)
end

%% Functions
function [alpha, yhat] = alphaCalc(tme, normintensity, alpha0)
    
    %convert negative values to NaN
    normintensity(normintensity<0) = NaN;
    
    %index traces that have at least 1 data point
    [nrow, ncol]=size(normintensity); 
    idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
    normintensity=normintensity(idx, :);
    
    %pre-allocate output variables
    [nrow, ~]=size(normintensity);
    alpha=nan(nrow,1);    
    yhat=nan(size(normintensity));
    
    %define exponential function
    modelfun = @(alpha, t)exp(-t./alpha);
    
    %fit data to the function and evaluate
    for i=1:nrow
        alpha(i, 1)=nlinfit(tme, normintensity(i,:), modelfun, alpha0);
         yhat(i, :)=modelfun(alpha(i), tme);
    end
    
end

function [rho, yhat] = rhoCalc(tme, normintensity, rho0)
    
    %convert negative values to NaN
    normintensity(normintensity<0) = NaN;
    
    %index traces that have at least 1 data point
    [nrow, ncol]=size(normintensity); 
    idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
    normintensity=normintensity(idx, :);
    
    %pre-allocate output variables
    [nrow, ncol]=size(normintensity);
    rho=nan(nrow,1);
    yhat=nan(size(normintensity));
    
    %calculate frame number
    dt = diff(tme, 1);
    frame = [0, tme(2:end)./dt];
    
    %define the exponential function
    modelfun = @(rho, x)exp(-x*rho);
    
    %fit data to find rho and evaluate
    for i=1:nrow
        rho(i, 1)=nlinfit(frame, normintensity(i,:), modelfun, rho0);
        yhat(i, :)=modelfun(rho(i, 1), frame);
    end
    
end

function [p] = meanPlot(x, y, color, transparency)

    %modified from the ciplot function (Raymond Reynolds 24/11/06)
    if nargin<3
        color = 'b';
        transparency = 1;
    end
    
    if nargin<4
        transparency = 1;
    end

    ymean = mean(y, 1, 'omitnan');
    ystd = std(y, 0, 1, 'omitnan');
    
    upper = ymean + ystd;
    lower = ymean - ystd;
 
    p = plot(x, ymean, 'Color', color, 'LineWidth', 1), hold on
    fill([x fliplr(x)],[upper fliplr(lower)], color,'LineStyle','none','FaceAlpha', transparency)

end

function [Cu] = photobleachSim(time, tau, A, beta, b)

    %diffusion = @(time, A, B, tau)(A * exp(-time./tau) + B);
    diffusionDerv = @(time, A, tau)(-A/tau * exp(-time./tau));
    photobleaching = @(Cu, alpha) Cu/alpha;
    dt = diff(time, 1);
    
    %pre-allocate variables
    Cu = nan(1, length(time));
    
    if tau == 0 
        
        Cu(1) = A;
        
        %'add' photobleaching   
        for i=1:length(time)-1
            alpha = (beta * dt(i)) + b; 
            dCb = photobleaching(Cu(i), alpha);
            dCu = -dCb * dt(i);
            Cu(i+1) = Cu(i) + dCu;
        end
    
    else
        
        Cu(1) = A;
        
        %'add' photobleaching   
        for i=1:length(time)-1
            alpha = (beta * dt(i)) + b; 
            dCb = photobleaching(Cu(i), alpha);
            dCp = diffusionDerv(time(i), A, tau);
            dCu = (-dCb + dCp) * dt(i);
            Cu(i+1) = Cu(i) + dCu;
        end

    end

end