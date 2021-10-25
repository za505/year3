%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

%% Compare sodium azide controls
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/10182021_analysis";
basename={'07242021_Exp1', '08172021_Exp1', '09222021_Exp1_colony', '10022021_Exp1_colony'};

%concerted = one frame rate, stepwise = multiple frame rates
cd(datadir)
concerted_files1=dir([basename{1} '*dm.mat']);
concerted_files2=dir([basename{2} '*dm.mat']);
stepwisefiles_noIPTG=dir([basename{3} '*dm.mat']);
stepwisefiles_IPTG=dir([basename{4} '*dm.mat']);

%% Stepwise Controls
%1 = noIPTG, 2 = IPTG

%pre-allocate matrices
stepwise_length1=[];
stepwise_intensity1=[];
stepwise_normintensity1=[];

stepwise_length2=[];
stepwise_intensity2=[];
stepwise_normintensity2=[];

cd(datadir)
%populate matrices
for i=1:height(stepwisefiles_noIPTG)
    
    load(stepwisefiles_noIPTG(i).name, 'dataTable', 'time', 'tidx')
    
    stepwise_length1=[stepwise_length1; table2array(dataTable(dataTable.halfie==0, 1))];
    stepwise_intensity1=[stepwise_intensity1; table2array(dataTable(dataTable.halfie==0, 2))];
    stepwise_normintensity1=[stepwise_normintensity1; table2array(dataTable(dataTable.halfie==0, 3))];

    stepwise_time1=time;
end

for i=1:height(stepwisefiles_IPTG)
    
    load(stepwisefiles_IPTG(i).name, 'dataTable', 'time', 'tidx')
    
    stepwise_length2=[stepwise_length2; table2array(dataTable(dataTable.halfie==0, 1))];
    stepwise_intensity2=[stepwise_intensity2; table2array(dataTable(dataTable.halfie==0, 2))];
    stepwise_normintensity2=[stepwise_normintensity2; table2array(dataTable(dataTable.halfie==0, 3))];

    stepwise_time2=time;
end

% plot length and intensity
cd(dirsave)

pretime=[0:2:15];
figure(1), hold on
for i=1:height(stepwise_length1)
    plot([pretime stepwise_time1+15], stepwise_length1(i,:), '-r')
end
for i=1:height(stepwise_length2)
    plot([pretime stepwise_time2+16], stepwise_length2(i,:), '-b')
end
xlabel('Time (min)')
ylabel('Length \mum')
title('Length vs Time')
saveas(gcf,'stepwise_length.png')
saveas(gcf, 'stepwise_length.fig')

figure(2), hold on
for i=1:height(stepwise_length1)
    plot(stepwise_time1, stepwise_intensity1(i,:), '-r')
end
for i=1:height(stepwise_length2)
    plot(stepwise_time2, stepwise_intensity2(i,:), '-b')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'stepwise_intensity.png')
saveas(gcf, 'stepwise_intensity.fig')

figure(3), hold on
for i=1:height(stepwise_length1)
    plot(stepwise_time1, stepwise_normintensity1(i,:), '-r')
end
for i=1:height(stepwise_length2)
    plot(stepwise_time2, stepwise_normintensity2(i,:), '-b')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'stepwise_normintensity.png')
saveas(gcf, 'stepwise_normintensity.fig')

%combine intensity plots and segment by frame rate
stepwise_normintensity = [stepwise_normintensity1; stepwise_normintensity2];

% fit data to exponential function
cd(dirsave)

segment1=stepwise_normintensity(:, 1:11);
segment2=stepwise_normintensity(:, 12:32);
segment3=stepwise_normintensity(:, 33:end);

%% fit to an exponential model to determine tau for each segment
beta=0.3;

%segment 1
mdl1=[];
y_hat1=segment1;
coeff0=1;
time=stepwise_time1(1:11);
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
saveas(gcf,'stepwise_segment1.png')
saveas(gcf, 'stepwise_segment1.fig')
    
%segment 2
mdl2=[];
y_hat2=segment2;
coeff0=1;
time=stepwise_time1(12:32)-stepwise_time1(12);
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
saveas(gcf,'stepwise_segment2.png')
saveas(gcf, 'stepwise_segment2.fig')

%segment 3
mdl3=[];
y_hat3=segment3;
coeff0=1;
time=stepwise_time1(33:end)-stepwise_time1(33);
for i=1:height(segment3)
    modelfun3=@(coeff,x)(segment3(i,1)-beta)*exp(-x./coeff)+beta;
    mdl3(i)=nlinfit(time, segment3(i,:), modelfun3, coeff0);
    y_hat3(i,:)=modelfun3(mdl3(i), time);   
end

figure, hold on
for i=1:height(segment3)
    plot(time, segment3(i,:),'-g',...
        time, y_hat3(i,:), '--k')
end
legend({'data', 'model'})
title('15 s frame rate')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'stepwise_segment3.png')
saveas(gcf, 'stepwise_segment3.fig')

%plot tau as a histogram
figure, hold on
histogram(mdl1(mdl1<1000), 'BinWidth', 20)
histogram(mdl2, 'BinWidth', 20)
histogram(mdl3, 'BinWidth', 20)
legend({'1 min', '30 s', '15 s'})
saveas(gcf,['stepwise_tau.png'])
saveas(gcf,['stepwise_tau.fig'])

%% Concentrated Controls
%pre-allocate matrices
concerted_length1=[];
concerted_intensity1=[];
concerted_normintensity1=[];

concerted_length2=[];
concerted_intensity2=[];
concerted_normintensity2=[];

%populate matrices
for i=1:height(concerted_files1)
    
    load(concerted_files1(i).name, 'dataTable', 'time', 'tidx')
    
    concerted_length1=[concerted_length1; table2array(dataTable(dataTable.halfie==0, 1))];
    concerted_intensity1=[concerted_intensity1; table2array(dataTable(dataTable.halfie==0, 2))];
    concerted_normintensity1=[concerted_normintensity1; table2array(dataTable(dataTable.halfie==0, 3))];

    concerted_time1=time;
end

for i=1:height(concerted_files2)
    
    load(concerted_files2(i).name, 'dataTable', 'time', 'tidx')
    
    concerted_length2=[concerted_length2; table2array(dataTable(dataTable.halfie==0, 1))];
    concerted_intensity2=[concerted_intensity2; table2array(dataTable(dataTable.halfie==0, 2))];
    concerted_normintensity2=[concerted_normintensity2; table2array(dataTable(dataTable.halfie==0, 3))];

    concerted_time2=time;
end

%what is the lowest fluorescent intensity the cells reach?
cd(dirsave)

pretime=[0:1:15];
figure(1), hold on
for i=1:height(concerted_length1)
    plot([pretime concerted_time1+16], concerted_length1(i,:), '-r')
end
for i=1:height(concerted_length2)
    plot([pretime concerted_time2+16], concerted_length2(i,:), '-b')
end
xlabel('Time (min)')
ylabel('Length \mum')
title('Length vs Time')
saveas(gcf,'concerted_length.png')
saveas(gcf, 'concerted_length.fig')

figure(2), hold on
for i=1:height(concerted_intensity1)
    plot(concerted_time1, concerted_intensity1(i,:), '-r')
end
for i=1:height(concerted_intensity2)
    plot(concerted_time2, concerted_intensity2(i,:), '-b')
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
saveas(gcf,'concerted_intensity.png')
saveas(gcf, 'concerted_intensity.fig')

figure(3), hold on
for i=1:height(concerted_normintensity1)
    plot(concerted_time1, concerted_normintensity1(i,:), '-r')
end
for i=1:height(concerted_normintensity2)
    plot(concerted_time2, concerted_normintensity2(i,:), '-b')
end
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
title('Normalized Intensity vs Time')
saveas(gcf,'concerted_normintensity.png')
saveas(gcf, 'concerted_normintensity.fig')

