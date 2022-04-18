%sanity_checks.m
%Author: Zarina Akbary
%Date: 04/15/2022
%Purpose: this is basically a sandbox for me to look through my data and
%keep track of how I do it. 

clear, close all

okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

%control1 = GFP alone, perfusion with LB + IPTG
%control2= GFP alone, perfusion with LB
%control3 = GFP + CY5, perfusion with LB + IPTG
%control4 = GFP + CY5, perfusion with LB

%% load data
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/aggregate/';
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/figures/';

control2_basenames = {'04052022_Exp2', '04052022_Exp1', '04042022_Exp1', '04112022_Exp1', '04052022_Exp3', '02122022_Exp1', '02122022_Exp2', '02092022_Exp1', '01282022_Exp1'};
control2_labels = {'Frame Rate = 10 s', 'Frame Rate = 20 s', 'Frame Rate = 1 min', 'Frame Rate = 1 min', 'Frame Rate = 2 min', 'Frame Rate = 5 min', 'Frame Rate = 10 min', 'Frame Rate = 20 min', 'Frame Rate = 20 min'};

control3_basenames = {'10302021_Exp2', '10302021_Exp1', '10232021_Exp1', '10262021_Exp1', '10282021_Exp1'};
control3_labels = {'Frame Rate = 15 s', 'Frame Rate = 30 s', 'Frame Rate = 1 min', 'Frame Rate = 1 min', 'Frame Rate = 5 min'};

[control2_cellTrace, control2_bgTrace, control2_adjTrace, control2_times, control2_lcell, control2_frameRate] = loadData(dirpath, control2_basenames);
[control3_cellTrace, control3_bgTrace, control3_adjTrace, control3_times, control3_lcell, control3_frameRate] = loadData(dirpath, control3_basenames);

%% What is the final background fluor. value vs time
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:height(control2_bgTrace)
    x = control2_times{i}(end);
    y = control2_bgTrace{i}(:, end);
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, 'o', 'MarkerFaceColor', okabeIto{i}, 'Color', okabeIto{i})
end
for i=1:height(control3_bgTrace)
    x = control3_times{i}(end);
    y = control3_bgTrace{i}(:, end);
    ybar = mean(y, 1, 'omitnan');
    err = std(y, 0, 1, 'omitnan');
    errorbar(x, ybar, err, '^', 'MarkerFaceColor', okabeIto{i}, 'Color', okabeIto{i})
end
legend([control2_labels, control3_labels])
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

function [normTrace] = dataNormalize(adjTrace)
    
    %pre-allocate variables
    normTrace = cell(size(adjTrace));
    
    %define function
    normInitial = @(f)f./f(:, 1);
    
    for i=1:height(adjTrace)  
        %normalize to the initial value
        normTrace{i} = normInitial(adjTrace{i});    
    end
    
end

function [alpha, alpha_yhat] = alphaCalc(times, normTrace, alpha0)
    
    %pre-allocate variables
    alpha = cell(height(normTrace), 1);
    alpha_yhat = cell(height(normTrace), 1);
    
    for i=1:height(normTrace)
        
        normintensity = normTrace{i};
        tme = times{i};
        
        %convert negative values to NaN
        normintensity(normintensity<0) = NaN;
    
        %index traces that have at least 1 data point
        [nrow, ncol]=size(normintensity); 
        idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
        normintensity=normintensity(idx, :);
    
        %pre-allocate output variables
        [nrow, ~]=size(normintensity);
        alphaTemp=nan(nrow,1);    
        yhatTemp=nan(size(normintensity));
    
        %define exponential function
        modelfun = @(alpha, t)exp(-t./alpha);
    
        %fit data to the function and evaluate
        for j=1:nrow
            alphaTemp(j, 1)=nlinfit(tme, normintensity(j,:), modelfun, alpha0);
            yhatTemp(j, :)=modelfun(alphaTemp(j), tme);
        end
        
        alpha{i} = alphaTemp;
        alpha_yhat{i} = yhatTemp;
        
    end    
    
end

function [rho, rho_yhat, frames] = rhoCalc(times, frameRates, normTrace, rho0)
    
    %pre-allocate variables
    rho = cell(height(normTrace), 1);
    rho_yhat = cell(height(normTrace), 1);
    frames = cell(height(normTrace), 1);
    
    for i=1:height(normTrace)
        
        normintensity = normTrace{i};
        tme = times{i};
        frameRate = frameRates{i};
        
        %calculate frame number
        frame = [0, tme(2:end)./frameRate(2:end)];
        
        %convert negative values to NaN
        normintensity(normintensity<0) = NaN;
    
        %index traces that have at least 1 data point
        [nrow, ncol]=size(normintensity); 
        idx=find(sum(isnan(normintensity), 2) < ncol-1); %the sum of the NaNs should be 1 less than the number of time points
        normintensity=normintensity(idx, :);
    
        %pre-allocate output variables
        [nrow, ~]=size(normintensity);
        rhoTemp=nan(nrow,1);    
        yhatTemp=nan(size(normintensity));
    
        %define exponential function
        modelfun = @(rho, x)exp(-x*rho);
        
        %fit data to the function and evaluate
        for j=1:nrow
            rhoTemp(j, 1)=nlinfit(frame, normintensity(j,:), modelfun, rho0);
            yhatTemp(j, :)=modelfun(rhoTemp(j, 1), frame);
        end
            
        rho{i} = rhoTemp;
        rho_yhat{i} = yhatTemp;
        frames{i} = frame;        
    end  
    
end