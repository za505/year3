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

[control2_cellTrace, control2_bgTrace, control2_times] = loadData(dirpath, control2_basenames);
[control3_cellTrace, control3_bgTrace, control3_times] = loadData(dirpath, control3_basenames);

%since these controls are all taken in one continuous movie, the first 8-10
%minutes can be ignored when subtracting the background value and
%normalizing. How necessary is it to normlize to the pre-lysis vs
%post-lysis value? How much of a difference does it make (mathematically)?
cutoff = 9; %point in the movie at which to start the analyze

control2_sliceTime = dataSlice(control2_times, cutoff);
control2_sliceIntensity = dataSlice(control2_cellTrace, cutoff);
control2_sliceBackground = dataSlice(control2_bgTrace, cutoff);

control3_sliceTime = dataSlice(control3_times, cutoff);
control3_sliceIntensity = dataSlice(control3_cellTrace, cutoff);
control3_sliceBackground = dataSlice(control3_bgTrace, cutoff);

alpha0=1;
rho0=1;

[control2_alpha, control2_alphaYhat] = alphaCalc(control2_sliceTime, control2_sliceIntensity, control2_sliceBackground, alpha0)
[control3_alpha, control3_alphaYhat] = alphaCalc(control3_sliceTime, control3_sliceIntensity, control3_sliceBackground, alpha0)

[control2_rho, control2_rhoYhat, control2_frames] = rhoCalc(control2_sliceTime, control2_sliceIntensity, control2_sliceBackground, rho0)
[control3_rho, control3_rhoYhat, control3_frames] = rhoCalc(control3_sliceTime, control3_sliceIntensity, control3_sliceBackground, rho0)

%% what do the sliced traces look like?
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(control2_basenames)
    subplot(3, 3, i)
    plot(times{i}, cellTrace{i}, 'Color', okabeIto{i}), hold on
    plot(times{i}, bgTrace{i}, 'LineStyle', '--', 'Color', okabeIto{i})
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    ylim([0 Inf])
    title(labels{i})
end

cd(dirsave)
saveas(gcf, 'figure02.png')
saveas(gcf, 'figure02.fig')

%% Functions
function [cellTrace, bgTrace, times] = loadData(dirpath, basenames)

    cellTrace = cell(length(basenames), 1);
    bgTrace = cell(length(basenames), 1);
    times = cell(length(basenames), 1);

    for i=1:length(basenames)

        cd(dirpath);
        basename = basenames{i};
        datadir=dir([basename '*']);
        load(datadir(1).name, '-regexp', 'intensity$', 'time');

        if exist('icell_intensity', 'var')
            cellTrace{i} = icell_intensity;
            bgTrace{i} = bg_intensity;
            times{i} = time;
            clear icell_intensity bg_intensity time
        else
            cellTrace{i} = intensity;
            bgTrace{i} = bgintensity;
            times{i} = time;
            clear intensity bgintensity time
        end

    end

end

function [slicedData] = dataSlice(data, cutoff)

    cutoff = 9; %point in the movie at which to start the analyze
    dataSlice = @(x, i)x(:, i:end)-x(:, i);
    
    slicedData = cell(size(data));
    
    for i=1:height(data)
        slicedData{i} = dataSlice(data{i}, cutoff);
    end
    
end

function [normTrace] = dataNormalize(cellTrace, bgTrace)
    
    %pre-allocate variables
    adjTrace = cell(size(cellTrace));
    normTrace = cell(size(cellTrace));
    
    %define function
    bgMinus = @(f, f1)f-f1;
    normInitial = @(f)f./f(:, 1);
    
    for i=1:height(cellTrace)
        %subtract the background
        adjTrace{i} = bgMinus(cellTrace{i}, bgTrace{i});
    
        %normalize to the initial value
        normTrace{i} = normInitial(adjTrace{i});    
    end
    
end

function [alpha, alpha_yhat] = alphaCalc(times, cellTrace, bgTrace, alpha0)
    
    %pre-allocate variables
    alpha = cell(height(cellTrace), 1);
    alpha_yhat = cell(height(cellTrace), 1);
    
    %normalize traces
    normTrace = dataNormalize(cellTrace, bgTrace);
    
    for h=1:height(normTrace)
        
        normintensity = normTrace{h};
        tme = times{h};
        
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
        for i=1:nrow
            alphaTemp(i, 1)=nlinfit(tme, normintensity(i,:), modelfun, alpha0);
            yhatTemp(i, :)=modelfun(alphaTemp(i), tme);
        end
        
        alpha{h} = alphaTemp;
        alpha_yhat{h} = yhatTemp;
        
    end    
    
end

function [rho, rho_yhat, frames] = rhoCalc(times, cellTrace, bgTrace, rho0)

    %pre-allocate variables
    rho = cell(height(cellTrace), 1);
    rho_yhat = cell(height(cellTrace), 1);
    frames = cell(height(cellTrace), 1);
    
    %normalize traces
    normTrace = dataNormalize(cellTrace, bgTrace);
    
    for h=1:height(normTrace)
        
        normintensity = normTrace{h};
        tme = times{h};
        
        %calculate frame number
        dt = diff(tme, 1);
        frame = [0, tme(2:end)./dt];
    
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
        for i=1:nrow
            rhoTemp(i, 1)=nlinfit(frame, normintensity(i,:), modelfun, rho0);
            yhatTemp(i, :)=modelfun(rhoTemp(i), frame);
        end
        
        rho{h} = rhoTemp;
        rho_yhat{h} = yhatTemp;
        frames = frame;
        
    end    
end