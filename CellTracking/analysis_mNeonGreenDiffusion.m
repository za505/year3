%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

%% user input and load data
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09222021_analysis";
basename='09182021_Exp1_colony';
eidx=79;

data = struct('experiment', [], 'colony', [], 'length', [], 'intensity', [], 'norm_intensity', [], 'time', [], 'fitModel', {});
cd(datadir)

files=dir([basename '*.mat']);

for i=1:height(files)
    
    data(i).experiment=files(i).name(1:end-15);
    data(i).colony=files(i).name(15:21);
    load(files(i).name, 'dataTable', 'time')
    data(i).length=table2array(dataTable(dataTable.halfie==0, 1));
    data(i).intensity=table2array(dataTable(dataTable.halfie==0, 2));
    data(i).norm_intensity=table2array(dataTable(dataTable.halfie==0, 3));
    
    data(i).length=data(i).length(:, 1:eidx);
    data(i).intensity=data(i).intensity(:, 1:eidx);
    data(i).norm_intensity=data(i).norm_intensity(:, 1:eidx);
    
    if time(1)~=0 
        time = time-time(1);
    end
  
    data(i).time=time(1:eidx)
    
end

%% plot length and intensity
cd(dirsave)

figure(1), hold on
for i=1:height(data)
    for j=1:height(data(i).length)
        plot(data(i).time, data(i).length(j,:), '-r')
    end
end

figure(2), hold on
for i=1:height(data)
    for j=1:height(data(i).intensity)
        plot(data(i).time, data(i).intensity(j,:), '-r')
    end
end

figure(3), hold on
for i=1:height(data)
    for j=1:height(data(i).norm_intensity)
        plot(data(i).time, data(i).norm_intensity(j,:), '-r')
    end
end

%% fit data to exponential function

