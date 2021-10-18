%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

%% 10022021_Exp1
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/10182021_analysis";
basename='10022021_Exp1_colony';

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
% plot length and intensity
cd(dirsave)

figure(1), hold on
for i=1:height(data_length)
    plot(length_time/60, data_length(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Length \mum')
title('Length vs Time')
saveas(gcf,'10022021_length.png')
saveas(gcf, '10022021_length.fig')

figure(2), hold on
for i=1:height(data_intensity)
    plot(intensity_time, data_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'10022021_intensity.png')
saveas(gcf, '10022021_intensity.fig')

figure(3), hold on
for i=1:height(data_norm_intensity)
    plot(intensity_time, data_norm_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'10022021_norm_intensity.png')
saveas(gcf, '10022021_norm_intensity.fig')

% fit data to exponential function
cd(dirsave)

segment1_avg=mean(data_intensity(:, 1));
segment1=data_intensity(:, 1:11)./segment1_avg;
segment1_initial=mean(segment1(:,1));

segment2_avg=mean(data_intensity(:, 12));
segment2=data_intensity(:, 12:32)./segment2_avg;
segment2_initial=mean(segment2(:,1));

segment3_avg=mean(data_intensity(:, 33));
segment3=data_intensity(:, 33:end)./segment3_avg;
segment3_initial=mean(segment3(:,1));

%segment 1
beta=0.4;
modelfun1=@(coeff,x)(segment1_initial-beta)*exp(-x./coeff)+beta;

mdl1=[];
y_hat1=segment1;
coeff0=0.1;
time=intensity_time(1:11);
for i=1:height(segment1)
    mdl1(i)=nlinfit(time, segment1(i,:), modelfun1, coeff0);
    y_hat1(i,:)=modelfun1(mdl1(i), time);   
end

figure, hold on
for i=1:height(segment1)
    plot(time, segment1(i,:),'-g',...
        time, y_hat1(i,:), '--k')
end
legend({'data', 'model'})
title('1 min frame rate')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10022021_segment1.png')
saveas(gcf, '10022021_segment1.fig')
    
%segment2
modelfun2=@(coeff,x)(segment2_initial-beta)*exp(-x./coeff)+beta;

mdl2=[];
y_hat2=segment2;
coeff0=0.1;
time=intensity_time(12:32)-intensity_time(12);
for i=1:height(segment2)
    mdl2(i)=nlinfit(time, segment2(i,:), modelfun2, coeff0);
    y_hat2(i,:)=modelfun2(mdl2(i), time);   
end

figure, hold on
for i=1:height(segment2)
    plot(time, segment2(i,:),'-g',...
        time, y_hat2(i,:), '--k')
end
legend({'data', 'model'})
title('30 second frame rate')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10022021_segment2.png')
saveas(gcf,'10022021_segment2.fig')

%segment3
modelfun3=@(coeff,x)(segment3_initial-beta)*exp(-x./coeff)+beta;

mdl3=[];
y_hat3=segment3;
coeff0=0.01;
time=intensity_time(33:end)-intensity_time(33);
for i=1:height(segment3)
    mdl3(i)=nlinfit(time, segment3(i,:), modelfun3, coeff0);
    y_hat3(i,:)=modelfun3(mdl3(i), time);   
end

figure, hold on
for i=1:height(segment3)
    plot(time, segment3(i,:),'-g', ...
        time, y_hat3(i,:), '--k')
end
legend({'data', 'model'})
title('15 second frame rate')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10022021_segment3.png')
saveas(gcf,'10022021_segment3.fig')

%plot tau as a histogram
figure, hold on
histogram(mdl1(mdl1<1000), 'BinWidth', 20)
histogram(mdl2, 'BinWidth', 20)
histogram(mdl3, 'BinWidth', 20)
legend({'1 min', '30 s', '15 s'})
saveas(gcf,['10022021_Exp1_tau.png'])
saveas(gcf,['10022021_Exp1_tau.fig'])
%% 09092021_Exp2
clear, close all

datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09302021_analysis";
basename='09092021_Exp2_colony';

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

%plot data
cd(dirsave)

figure(1), hold on
for i=1:height(data_length)
    plot(length_time/60, data_length(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Length \mum')
title('Length vs Time')
saveas(gcf,'length_16min.png')
saveas(gcf, 'length_16min.fig')

figure(2), hold on
for i=1:height(data_intensity)
    plot(intensity_time, data_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'intensity_16min.png')
saveas(gcf, 'intensity_16min.fig')

figure(3), hold on
for i=1:height(data_norm_intensity)
    plot(intensity_time, data_norm_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'norm_intensity_16min.png')
saveas(gcf, 'norm_intensity_16min.fig')

%% 09222021_Exp2
clear, close all

datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09302021_analysis";
basename='09222021_Exp2_colony';

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

%plot data
cd(dirsave)

figure(1), hold on
for i=1:height(data_length)
    plot(length_time/60, data_length(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Length \mum')
title('Length vs Time')
saveas(gcf,'length_20min.png')
saveas(gcf, 'length_20min.fig')

figure(2), hold on
for i=1:height(data_intensity)
    plot(intensity_time, data_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'intensity_20min.png')
saveas(gcf, 'intensity_20min.fig')

figure(3), hold on
for i=1:height(data_norm_intensity)
    plot(intensity_time, data_norm_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'norm_intensity_20min.png')
saveas(gcf, 'norm_intensity_20min.fig')

%% 09232021_Exp2
clear, close all

datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09302021_analysis";
basename='09232021_Exp2_colony';

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

%plot data
cd(dirsave)

figure(1), hold on
for i=1:height(data_length)
    plot(length_time/60, data_length(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Length \mum')
title('Length vs Time')
saveas(gcf,'length_20min_Mg.png')
saveas(gcf, 'length_20min_Mg.fig')

figure(2), hold on
for i=1:height(data_intensity)
    plot(intensity_time, data_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'intensity_20min_Mg.png')
saveas(gcf, 'intensity_20min_Mg.fig')

figure(3), hold on
for i=1:height(data_norm_intensity)
    plot(intensity_time, data_norm_intensity(i,:), '-g')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'norm_intensity_20min_Mg.png')
saveas(gcf, 'norm_intensity_20min_Mg.fig')

%% 09282021_Exp1
clear, close all

datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09302021_analysis";
basename='09282021_Exp1_colony';

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
% length_files=dir([basename '*BTphase.mat']);
% 
% %pre-allocate matrices
% data_length=[];
% 
% %populate matrices
% for i=1:height(length_files)
%     
%     load(length_files(i).name, 'lcell', 'time')
%     
%     data_length=[data_length; lcell];
% 
%     length_time=time;
% end

%plot data
cd(dirsave)

% figure(1), hold on
% for i=1:height(data_length)
%     plot(length_time/60, data_length(i,:), '-g')
% end
% xlabel('Time (min)')
% ylabel('Length \mum')
% title('Length vs Time')
% saveas(gcf,'length_20min_Mg.png')
% saveas(gcf, 'length_20min_Mg.fig')

figure(2), hold on
for i=1:height(data_intensity)
    plot(intensity_time, data_intensity(i,:), '-r')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'intensity_20%.png')
saveas(gcf, 'intensity_20%.fig')

figure(3), hold on
for i=1:height(data_norm_intensity)
    plot(intensity_time, data_norm_intensity(i,:), '-r')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'norm_intensity_20%.png')
saveas(gcf, 'norm_intensity_20%.fig')

%% 09282021_Exp2
clear, close all

datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09302021_analysis";
basename='09282021_Exp2_colony';

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
% length_files=dir([basename '*BTphase.mat']);
% 
% %pre-allocate matrices
% data_length=[];
% 
% %populate matrices
% for i=1:height(length_files)
%     
%     load(length_files(i).name, 'lcell', 'time')
%     
%     data_length=[data_length; lcell];
% 
%     length_time=time;
% end

%plot data
cd(dirsave)

% figure(1), hold on
% for i=1:height(data_length)
%     plot(length_time/60, data_length(i,:), '-g')
% end
% xlabel('Time (min)')
% ylabel('Length \mum')
% title('Length vs Time')
% saveas(gcf,'length_20min_Mg.png')
% saveas(gcf, 'length_20min_Mg.fig')

figure(2), hold on
for i=1:height(data_intensity)
    plot(intensity_time, data_intensity(i,:), '-r')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'intensity_5%.png')
saveas(gcf, 'intensity_5%.fig')

figure(3), hold on
for i=1:height(data_norm_intensity)
    plot(intensity_time, data_norm_intensity(i,:), '-r')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'norm_intensity_5%.png')
saveas(gcf, 'norm_intensity_5%.fig')