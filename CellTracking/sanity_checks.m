%sanity_checks.m
%Author: Zarina Akbary
%Date: 04/15/2022
%Purpose: this is basically a sandbox for me to look through my data and
%keep track of how I do it. 

clear, close all

okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);
okabeIto = [okabeIto, okabeIto];

%% for the 20% intensity controls
dirpath = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/aggregate/';
dirsave = '/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/figures/';
basenames = {'02192022_Exp2', '02192022_Exp3', '02192022_Exp4', '02212022_Exp2', '02212022_Exp1'};
labels = {'Frame Rate = 1 s', 'Frame Rate = 2 s', 'Frame Rate = 3 s', 'Frame Rate = 1 min', 'Frame Rate = 20 min'};

cellTrace = cell(length(basenames), 1);
bgTrace = cell(length(basenames), 1);
times = cell(length(basenames), 1);

for i=1:length(basenames)
    
    cd(dirpath);
    basename = basenames{i}
    datadir=dir([basename '*']);
    datadir(1).name
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

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize', 15), hold on
for i=1:length(basenames)
    subplot(1, 5, i)
    plot(times{i}, cellTrace{i}, 'Color', okabeIto{i}), hold on
    plot(times{i}, bgTrace{i}, 'LineStyle', '--', 'Color', okabeIto{i})
    xlabel('Time (minutes)')
    ylabel('Fluorescence (A.U.)')
    ylim([0 Inf])
    title(labels{i})
end

cd(dirsave)
saveas(gcf, 'figure01.png')
saveas(gcf, 'figure01.fig')