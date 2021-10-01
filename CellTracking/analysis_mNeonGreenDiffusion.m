%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

%% load sodium azide
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09302021_analysis";
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
xlabel('Time (min)')
ylabel('Length \mum')
title('Length vs Time')
saveas(gcf,'length.png')
saveas(gcf, 'length.fig')

figure(2), hold on
for i=1:height(data_intensity)
    plot(intensity_time, data_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'intensity.png')
saveas(gcf, 'intensity.fig')

figure(3), hold on
for i=1:height(data_norm_intensity)
    plot(intensity_time, data_norm_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'norm_intensity.png')
saveas(gcf, 'norm_intensity.fig')

%% fit data to exponential function
cd(dirsave)

segment1=data_norm_intensity(:, 1:11);
segment2=data_norm_intensity(:, 12:32);
segment3=data_norm_intensity(:, 33:end);

beta=0.5;
modelfun=@(coeff,x)(1-beta)*exp(-x./coeff)+beta;

%segment 1
mdl1={};
y_hat1=segment1;
coeff0=1;
time=intensity_time(1:11);
for i=1:height(segment1)
    mdl1{i}=nlinfit(time, segment1(i,:), modelfun, coeff0);
    y_hat1(i,:)=modelfun(mdl1{i}, time);   
end

figure, hold on
for i=1:height(segment1)
    plot(time, segment1(i,:),'-g',...
        time, y_hat1(i,:), '--k')
end
legend({'data', 'model'})
title(['1 min frame rate, tau = ' num2str(coeff0)])
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,['segment1_' num2str(coeff0) '.png'])
saveas(gcf, ['segment1_' num2str(coeff0) '.fig'])
    
%segment2
mdl2={};
y_hat2=segment2;
coeff0=0.1;
time=intensity_time(12:32)-intensity_time(12);
for i=1:height(segment2)
    mdl2{i}=nlinfit(time, segment2(i,:), modelfun, coeff0);
    y_hat2(i,:)=modelfun(mdl2{i}, time);   
end

figure, hold on
for i=1:height(segment2)
    plot(time, segment2(i,:),'-g',...
        time, y_hat2(i,:), '--k')
end
legend({'data', 'model'})
title(['30 second frame rate, tau = ' num2str(coeff0)])
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,['segment2_' num2str(coeff0) '.png'])
saveas(gcf,['segment2_' num2str(coeff0) '.fig'])

%segment3
mdl3={};
y_hat3=segment3;
coeff0=0.01;
time=intensity_time(33:end)-intensity_time(33);
for i=1:height(segment3)
    mdl3{i}=nlinfit(time, segment3(i,:), modelfun, coeff0);
    y_hat3(i,:)=modelfun(mdl3{i}, time);   
end

figure, hold on
for i=1:height(segment3)
    plot(time, segment3(i,:),'-g', ...
        time, y_hat3(i,:), '--k')
end
legend({'data', 'model'})
title(['15 second frame rate, tau = ' num2str(coeff0)])
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,['segment3_' num2str(coeff0) '.png'])
saveas(gcf,['segment3_' num2str(coeff0) '.fig'])