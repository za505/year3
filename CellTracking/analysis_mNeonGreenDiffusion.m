%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all

%% user input
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09202021_analysis";
presdir="/Users/zarina/Downloads/NYU/Year3_2021_Fall/updates/09142021";

data = struct('experiment', [], 'colony', [], 'length', [], 'intensity', [], 'norm_intensity', [], 'time', [], 'fitModel', {});
cd(datadir)

files=dir('*.mat');

experiments=[];
colonies=[];

for i=1:height(files)
    
    data(i).experiment=files(i).name(1:end-15);
    data(i).colony=files(i).name(15:21);
    load(files(i).name, 'dataTable', 'time')
    data(i).length=table2array(dataTable(dataTable.halfie==0, 1));
    data(i).intensity=table2array(dataTable(dataTable.halfie==0, 2));
    data(i).norm_intensity=table2array(dataTable(dataTable.halfie==0, 3));
    
    if time(1)~=0 
        time = time-time(1);
    end
  
    data(i).time=time;
    
    line1=convertCharsToStrings(data(i).experiment);
    line2=convertCharsToStrings(data(i).colony);
    
    mat1=repmat(line1, height(data(i).length),1);
    mat2=repmat(line2, height(data(i).length),1);
    
    experiments=[experiments; mat1];
    colonies=[colonies;mat2];
    
end

%% Notes
%Mg2+, time = [0:2:180]
%LB, frame rate = 2, time = [0:2:180]
%EDTA, time = [0:2:180]
%LB + sodium azide, frame rate = 2, time = [0:2:180]
%LB + sodium azide, frame rate = 2, time = [0:2:180]
%LB, frame rate = 1, time = [0:91]
%PBS, frame rate = 1, incubation = 2 min, time = [0:90]
%PBS, frame rate = 1, incubation = 8 min, time = [0:1:45]
%PBS, frame rate = 1, incubation = 16 min, time = [0:1:45]
%PBS, frame rate = 1, incubation = 20 min, time = [0:1:45]

%annotation=["Mg2+", "LB, fr=2", "EDTA", "sodium azide", "sodium azide", "LB, fr=1", "PBS, 2 min", "PBS, 8 min", "PBS, 16 min"];
SA1=116:123;
SA2=124:151;
PBS1=239:340;
PBS2=341:417;
PBS3=418:462;
%% fit data to exp2
for i=26:46 %1:length(data)
    for j=1:height(data(i).norm_intensity)
        j
        if isnan(data(i).norm_intensity(j,:))==0
            [xData, yData]=prepareCurveData(data(i).time/60, data(i).norm_intensity(j,:));
            f = fit(xData, yData, '(1-B)*exp(-x/tau)+B', 'Lower', [0, 0], 'Upper', [2, Inf], 'StartPoint', [0.3,0.8]);
            
            data(i).fitModel{j,1}=f;
            
%             figure
%             plot(data(i).fitModel{j,1}, data(i).time/60, data(i).norm_intensity(j,:))
%             title([data(i).experiment ',' data(i).colony])
%             xlabel('Time (hours)')
%             ylabel('Normalized Intensity (A.U.)')
%             pause, close
        else
            [data(i).fitModel{j,1}] = [];
        end
    end
end

%% plot tau for each distribution, then plot tau over w (where w is incubation time in PBS)
pbs1=26:31;
pbs2=32:37;
pbs3=38:42;
pbs4=43:46;

cd(dirsave)

figure, hold on  
for i=pbs1
    for j=1:height(data(i).norm_intensity)
        if isnan(data(i).norm_intensity(j,:))==0
           plot(data(i).time, data(i).norm_intensity(j,:), 'r')
        end
    end
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Time vs Normalized Fluorescence, 2 min incubation')
saveas(gcf, 'pbs_2min.png')
saveas(gcf, 'pbs_2min.fig')

figure, hold on  
for i=pbs2
    for j=1:height(data(i).norm_intensity)
        if isnan(data(i).norm_intensity(j,:))==0
           plot(data(i).time, data(i).norm_intensity(j,:), 'b')
        end
    end
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Time vs Normalized Fluorescence, 8 min incubation')
saveas(gcf, 'pbs_8min.png')
saveas(gcf, 'pbs_8min.fig')

figure, hold on  
for i=pbs3
    for j=1:height(data(i).norm_intensity)
        if isnan(data(i).norm_intensity(j,:))==0
           plot(data(i).time, data(i).norm_intensity(j,:), 'g')
        end
    end
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Time vs Normalized Fluorescence, 16 min incubation')
saveas(gcf, 'pbs_16min.png')
saveas(gcf, 'pbs_16min.fig')

figure, hold on  
for i=pbs4
    for j=1:height(data(i).norm_intensity)
        if isnan(data(i).norm_intensity(j,:))==0
           plot(data(i).time, data(i).norm_intensity(j,:), 'k')
        end
    end
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Time vs Normalized Fluorescence, 20 min incubation')
saveas(gcf, 'pbs_20min.png')
saveas(gcf, 'pbs_20min.fig')

[b1, tau1] = coeffExt(data, pbs1)
figure(1), hold on
h1=histogram(tau1)
h1.FaceColor = 'r';
h1.EdgeColor = 'k';
h1.BinWidth = 0.1;
xlabel('tau')
ylabel('Counts')
title('Distribution of tau Coefficient for 2 min incubation')
saveas(gcf, 'tau1_2min.png')
saveas(gcf, 'tau1_2min.fig')

[b2, tau2] = coeffExt(data, pbs2)
figure(2), hold on
h2=histogram(tau2)
h2.FaceColor = 'b';
h2.EdgeColor = 'k';
h2.BinWidth = 0.1;
xlabel('tau')
ylabel('Counts')
title('Distribution of tau Coefficient for 8 min incubation')
saveas(gcf, 'tau1_8min.png')
saveas(gcf, 'tau1_8min.fig')

[b3, tau3] = coeffExt(data, pbs3)
figure(3), hold on
h3=histogram(tau3)
h3.FaceColor = 'c';
h3.EdgeColor = 'k';
h3.BinWidth = 0.1;
xlabel('tau')
ylabel('Counts')
title('Distribution of tau Coefficient for 16 min incubation')
saveas(gcf, 'tau1_16min.png')
saveas(gcf, 'tau1_16min.fig')

[b4, tau4] = coeffExt(data, pbs4)
figure(4), hold on
h4=histogram(tau4)
h4.FaceColor = 'm';
h4.EdgeColor = 'k';
h4.BinWidth = 0.1;
xlabel('tau')
ylabel('Counts')
title('Distribution of tau Coefficient for 20 min incubation')
saveas(gcf, 'tau1_20min.png')
saveas(gcf, 'tau1_20min.fig')

tau=[tau1; tau2; tau3; tau4];
w=[repelem(2/60, length(tau1))'; repelem(8/60, length(tau2))'; repelem(16/60, length(tau3))'; repelem(20/60, length(tau4))'];

figure(5)
scatter(w, tau)
xlabel('omega (incubation time in hours)')
ylabel('tau (h^-1)')
saveas(gcf, 'tau_omega.png')
saveas(gcf, 'tau_omega.fig')

tau_avg=[mean(tau1); mean(tau2); mean(tau3); mean(tau4)];
tau_std=[std(tau1); std(tau2); std(tau3); std(tau4)];
w_avg=[2; 8; 16; 20];

figure(6), hold on
% errorbar(w_avg, tau_avg,-tau_std,tau_std,'o')
errorbar(2, mean(tau1),-std(tau1), std(tau1),'o')
errorbar(8, mean(tau2),-std(tau2), std(tau2),'o')
errorbar(16, mean(tau3),-std(tau3), std(tau3),'o')
errorbar(20, mean(tau4),-std(tau4), std(tau4),'o')
legend({'2 min', '8 min', '16 min', '20 min'})
xlabel('omega (incubation time in minutes)')
ylabel('average tau (h^-1)')
xlim([0, 22])
saveas(gcf, 'tau_omega_avg.png')
saveas(gcf, 'tau_omega_avg.fig')

%% save data
cd(dirsave)
save('fitting.mat')

%% Sodium Azide Experiment
clear, close all
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/09202021_analysis";
cd('/Users/zarina/Downloads/NYU/Year3_2021_Fall/09182021_analysis/09182021_Exp1/09182021_Exp1_colony1/09182021_Exp1_mNeonGreen/09182021_Exp1_figures')
load('09182021_Exp1_dm.mat')

cd(dirsave)
figure(7), hold on
for i=1:height(norm_intensity)
    plot(time, norm_intensity(i,:), '-b')
end
xlabel('Time (minutes)')
ylabel('Intensity (A.U.)')
saveas(gcf, 'azide_normintensity.png')
saveas(gcf, 'azide_normintensity.fig')

figure(8), hold on
for i=1:height(icell_intensity)
    plot(time, icell_intensity(i,:), '-b')
end
xlabel('Time (minutes)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf, 'azide_intensity.png')
saveas(gcf, 'azide_intensity.fig')
 
dm=nan(ncells, T-1);
dx = diff(norm_intensity, 1, 2);
dy = diff(time/60);
for i=1:height(dx)
    dm(i, :) = dy./dx(i,:);
    idx = find(abs(dm(i,:))>10);
    dm(i, idx)=NaN;
end

dm_avg = mean(dm, 1, 'omitnan');
p1 = [mean(dm_avg(:, 4:6), 2, 'omitnan'), std(dm_avg(:, 4:6), 0, 2, 'omitnan')];
p2 = [mean(dm_avg(:, 7:11), 2, 'omitnan'), std(dm_avg(:, 7:11), 0, 2, 'omitnan')];
p3 = [mean(dm_avg(:, 12:22), 2, 'omitnan'), std(dm_avg(:, 12:22), 0, 2, 'omitnan')];
p4 = [mean(dm_avg(:, 23:41), 2, 'omitnan'), std(dm_avg(:, 23:41), 0, 2, 'omitnan')];

dm_comp = [p1;p2;p3;p4];

figure(9), hold on
for i=1:ncells
   scatter(time(1,1:end-1)/60, dm(i,:), 'b')
end

figure(10)
%scatter(time(1,1:end-1)/60, mean(dm, 1))
errorbar(time(1,1:end-1)/60, mean(dm, 1),-std(dm, 0,1),std(dm, 0,1),'o')
ylabel('Slope, Normalized Intensity (A.U.)/Time (h)')
xlabel('Time (h)')

%fr = [repelem(120, 6), repelem(60, 6), repelem(30, 11), repelem(15,21), repelem(7.5,41)];
fr = [120, 60, 30, 15, 7.5];
figure(11), hold on
%scatter(time(1,1:end-1)/60, mean(dm, 1))
errorbar(fr(1), p1(1), -p1(2), p1(2),'o')
errorbar(fr(2), p2(1), -p2(2), p2(2),'o')
errorbar(fr(3), p3(1), -p3(2), p3(2),'o')
errorbar(fr(4), p4(1), -p4(2), p4(2),'o')
ylabel('Slope, Normalized Intensity (A.U.)/Time (h)')
xlabel('Frame Rate (s)')
xlim([0, 122])
legend({'120 s', '60 s', '30 s', '15 s', '7.5 s'})
title('Slope as a Function of Frame Rate')
saveas(gcf, 'azide_fr.png')
saveas(gcf, 'azide_fr.fig')
%% functions
function [b, tau]=coeffExt(data, range)
    b=[];
    tau=[];
    for i=range
        for j=1:height(data(i).fitModel)
            if isempty(data(i).fitModel{j,1})==0
                coeffs = coeffvalues(data(i).fitModel{j,1})
                b=[b; coeffs(1)];
                tau=[tau; coeffs(2)];
            end
        end
    end
end

function [th] = halfLife(data, range)
    figure, hold on
    
    for i=range
        for j=1:height(data(i).norm_intensity)
            if isnan(data(i).norm_intensity(j,:))==0
               plot(data(i).time, data(i).norm_intensity(j,:))
            end
        end
    end
    
    pause, close
end