%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

%% load sodium azide
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09272021_analysis";
basename='09222021_Exp1_colony';

%store directory
%data = struct('experiment', [], 'colony', [], 'length', [], 'intensity', [], 'norm_intensity', [], 'time', [], 'fitModel', {});
cd(datadir)
intensity_files=dir([basename '*dm.mat']);

%pre-allocate matrices
data_intensity=[];
data_norm_intensity=[];

%populate matrices
for i=1:height(intensity_files)
    
    load(intensity_files(i).name, 'dataTable', 'time', 'tidx')
    
    data_intensity=[data_intensity; table2array(dataTable(dataTable.halfie==0, 2))];
    data_norm_intensity=[data_norm_intensity; table2array(dataTable(dataTable.halfie==0, 3))];

    intensity_time=time;
end

%store directory
length_files=dir([basename '*BTphase.mat']);

%pre-allocate matrices
data_length=[];

%populate matrices
for i=1:height(length_files)
    
    load(length_files(i).name, 'lcell', 'time')
    
    data_length=[data_length; lcell];

    length_time=time;
end
%% plot length and intensity
cd(dirsave)

figure(1), hold on
for i=1:height(data_length)
    plot(length_time/60, data_length(i,:), '-g')
end

figure(2), hold on
for i=1:height(data_intensity)
    plot(intensity_time, data_intensity(i,:), '-g')
end

figure(3), hold on
for i=1:height(data_norm_intensity)
    plot(intensity_time, data_norm_intensity(i,:), '-g')
end

%% fit data to exponential function
segment1=data_norm_intensity(:, 1:11);
segment2=data_norm_intensity(:, 12:32);
segment3=data_norm_intensity(:, 33:end);

beta=0.99;
modelfun=@(coeff,x)(1-beta)*exp(-x./coeff)+beta;

coeff0=0.1;
time=intensity_time(1:11);
for i=1:height(segment1)
    mdl1=fitnlm(time, segment1(i,:), modelfun, coeff0);
    y_hat=modelfun(mdl1.Coefficients, time)
    figure
    plot(time, segment1(i,:),...
        time, y_hat)
    pause, close
end
