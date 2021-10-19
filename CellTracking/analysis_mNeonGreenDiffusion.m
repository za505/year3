%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

%% Mock Model--30 min, full
beta=0.7;
t1=[0:1:30];
tau1=100;

t2=[0:0.5:30];
tau2=10;

t3=[0:0.25:30];
tau3=1;

sim1 = (1-beta)*exp(-t1/tau1)+beta;
sim2 = (1-beta)*exp(-t2/tau2)+beta;
sim3 = (1-beta)*exp(-t3/tau3)+beta;

figure, hold on
plot(t1, sim1)
plot(t2, sim2)
plot(t3, sim3)

%% Mock Model--30 min stepwise
beta=0.4;
t1=[0:1:10];
tau1=250;

t2=[0:0.5:10];
tau2=100;

t3=[0:0.25:10];
tau3=10;

sim1 = (1-beta)*exp(-t1/tau1)+beta;
sim2 = (sim1(end)-beta)*exp(-t2/tau2)+beta;
sim3 = (sim2(end)-beta)*exp(-t3/tau3)+beta;

simT = [sim1, sim2, sim3];
time = [t1, t2+(t1(end)+0.5), t3+(t1(end)+t2(end)+0.25)];

figure, hold on
plot(time, simT)

%% Working backwards
beta=0.7;
coeff0=1;

modelfun1=@(coeff,x)(sim1(1)-beta)*exp(-x./coeff)+beta;
modelfun2=@(coeff,x)(sim2(1)-beta)*exp(-x./coeff)+beta;
modelfun3=@(coeff,x)(sim3(1)-beta)*exp(-x./coeff)+beta;

mdl1=nlinfit(t1, sim1, modelfun1, coeff0);
mdl2=nlinfit(t2, sim2, modelfun2, coeff0);
mdl3=nlinfit(t3, sim3, modelfun3, coeff0);

y_hat1=modelfun1(mdl1, t1); 
y_hat2=modelfun2(mdl2, t2); 
y_hat3=modelfun3(mdl3, t3); 

figure, hold on
plot(t1, sim1,'-g', ...
        t1, y_hat1, '--k')

figure, hold on
plot(t2, sim2,'-g', ...
        t2, y_hat2, '--k')

figure, hold on
plot(t3, sim3,'-g', ...
        t3, y_hat3, '--k')

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
% 
% figure(1), hold on
% for i=1:height(data_length)
%     plot(length_time/60, data_length(i,:), '-g')
% end
% xlabel('Time (min)')
% ylabel('Length \mum')
% title('Length vs Time')
% saveas(gcf,'10022021_length.png')
% saveas(gcf, '10022021_length.fig')
% 
% figure(2), hold on
% for i=1:height(data_intensity)
%     plot(intensity_time, data_intensity(i,:), '-g')
% end
% xlabel('Time (min)')
% ylabel('Intensity (A.U.)')
% title('Intensity vs Time')
% saveas(gcf,'10022021_intensity.png')
% saveas(gcf, '10022021_intensity.fig')
% 
% figure(3), hold on
% for i=1:height(data_norm_intensity)
%     plot(intensity_time, data_norm_intensity(i,:), '-g')
% end
% xlabel('Time (min)')
% ylabel('Normalized Intensity (A.U.)')
% title('Normalized Intensity vs Time')
% saveas(gcf,'10022021_norm_intensity.png')
% saveas(gcf, '10022021_norm_intensity.fig')

% fit data to exponential function
cd(dirsave)

segment1=data_norm_intensity(:, 1:10);
segment2=data_norm_intensity(:, 11:31);
segment3=data_norm_intensity(:, 32:end);

%segment 1
beta=0.3;

mdl1=[];
y_hat1=segment1;
coeff0=1;
time=intensity_time(1:10);
for i=1:height(segment1)
    modelfun1=@(coeff,x)(segment1(i,1)-beta)*exp(-x./coeff)+beta;
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
mdl2=[];
y_hat2=segment2;
time=intensity_time(11:31)-intensity_time(11);
for i=1:height(segment2)
    modelfun2=@(coeff,x)(segment2(i,1)-beta)*exp(-x./coeff)+beta;
    mdl2(i)=nlinfit(time, segment2(i,:), modelfun2, coeff0);
    y_hat2(i,:)=modelfun2(mdl2(i), time);   
end

figure, hold on
for i=1:height(segment2)
    plot(time, segment2(i,:),'-g',...
        time, y_hat2(i,:), '--k')
end
legend({'data', 'model'})
title('30 s frame rate')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10022021_segment2.png')
saveas(gcf, '10022021_segment2.fig')

% %segment3
mdl3=[];
y_hat3=segment3;
time=intensity_time(32:end)-intensity_time(32);
for i=1:height(segment3)
    modelfun3=@(coeff,x)(segment3(i,1)-beta)*exp(-x./coeff)+beta;
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
