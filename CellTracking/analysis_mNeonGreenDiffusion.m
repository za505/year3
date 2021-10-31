%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

% load the data
dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/10282021_analysis";
cd(dirsave)

datadir1=dir(['10232021_Exp1' '*dm.mat']); %LB, frame rate=1 min
datadir2=dir(['10232021_Exp2' '*dm.mat']); %PBS, frame rate=1 min
datadir3=dir(['10262021_Exp1' '*dm.mat']); %LB, frame rate=1 min
datadir4=dir(['10262021_Exp2' '*dm.mat']); %PBS, frame rate=1 min

intensity1=[];
adjintensity1=[];
normintensity1=[];
positions1=[0];
for i=1:length(datadir1)
    load(datadir1(i).name, 'icell_intensity', 'adj_intensity', 'norm_intensity', 'time')
    intensity1=[intensity1; icell_intensity];
    adjintensity1=[adjintensity1; adj_intensity];
    normintensity1=[normintensity1; norm_intensity];
    positions1(i+1)=height(norm_intensity)+positions1(i);
end
time1=time;

intensity2=[];
adjintensity2=[];
normintensity2=[];
positions2=[0];
for i=1:length(datadir2)
    load(datadir2(i).name, 'icell_intensity', 'adj_intensity', 'norm_intensity', 'time')
    intensity2=[intensity2; icell_intensity];
    adjintensity2=[adjintensity2; adj_intensity];
    normintensity2=[normintensity2; norm_intensity];
    positions2(i+1)=height(norm_intensity)+positions2(i);
end
time2=time;

intensity3=[];
adjintensity3=[];
normintensity3=[];
positions3=[0];
for i=1:length(datadir3)
    load(datadir3(i).name, 'icell_intensity', 'adj_intensity', 'norm_intensity', 'time')
    intensity3=[intensity3; icell_intensity];
    adjintensity3=[adjintensity3; adj_intensity];
    normintensity3=[normintensity3; norm_intensity];
    positions3(i+1)=height(norm_intensity)+positions3(i);
end
time3=time;

intensity4=[];
adjintensity4=[];
normintensity4=[];
positions4=[0];
for i=1:length(datadir4)
    load(datadir4(i).name, 'icell_intensity', 'adj_intensity', 'norm_intensity', 'time')
    intensity4=[intensity4; icell_intensity];
    adjintensity4=[adjintensity4; adj_intensity];
    normintensity4=[normintensity4; norm_intensity];
    positions4(i+1)=height(norm_intensity)+positions4(i);
end
time4=time;

% plot the data
figure(1), hold on
for i=1:height(intensity1)
    plot(time1(1:98), intensity1(i,1:98), '-r')
end
for i=1:height(intensity3)
    plot(time3, intensity3(i,:), '-b')
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
title('Diffusion in LB')
text(75, 5000, 'replicate 1');
text(68, 5050, '_____', 'Color', 'red');
text(75, 4800, 'replicate 2');
text(68, 4850, '_____', 'Color', 'blue');
saveas(gcf,'LB_intensity.png')
saveas(gcf, 'LB_intensity.fig')

figure(2), hold on
for i=1:height(intensity2)
    plot(time2, intensity2(i,:), '-r')
end
for i=1:height(intensity4)
    plot(time4(1:51), intensity4(i,1:51), '-b')
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
title('Diffusion After Incubation in PBS')
text(35, 5000, 'replicate 1');
text(32, 5050, '_____', 'Color', 'red');
text(35, 4600, 'replicate 2');
text(32, 4650, '_____', 'Color', 'blue');
saveas(gcf,'PBS_intensity.png')
saveas(gcf, 'PBS_intensity.fig')

figure(3), hold on
for i=1:height(adjintensity1)
    plot(time1(1:98), adjintensity1(i,1:98), '-r')
end
for i=1:height(adjintensity3)
    plot(time3, adjintensity3(i,:), '-b')
end
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
title('Diffusion in LB')
text(75, 4000, 'replicate 1');
text(68, 4050, '_____', 'Color', 'red');
text(75, 3800, 'replicate 2');
text(68, 3850, '_____', 'Color', 'blue');
saveas(gcf,'LB_adjintensity.png')
saveas(gcf, 'LB_adjintensity.fig')

figure(4), hold on
for i=1:height(adjintensity2)
    plot(time2, adjintensity2(i,:), '-r')
end
for i=1:height(adjintensity4)
    plot(time4(1:51), adjintensity4(i,1:51), '-b')
end
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence (A.U.)')
title('Diffusion After Incubation in PBS')
text(35, 4000, 'replicate 1');
text(32, 4050, '_____', 'Color', 'red');
text(35, 3800, 'replicate 2');
text(32, 3850, '_____', 'Color', 'blue');
saveas(gcf,'PBS_adjintensity.png')
saveas(gcf, 'PBS_adjintensity.fig')

figure(5), hold on
for i=1:height(normintensity1)
    plot(time1(1:98), normintensity1(i,1:98), '-r')
end
for i=1:height(normintensity3)
    plot(time3, normintensity3(i,:), '-b')
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion in LB')
text(75, .8, 'replicate 1');
text(68, .82, '_____', 'Color', 'red');
text(75, .75, 'replicate 2');
text(68, .77, '_____', 'Color', 'blue');
saveas(gcf,'LB_normintensity.png')
saveas(gcf, 'LB_normintensity.fig')

figure(6), hold on
for i=1:height(normintensity2)
    plot(time2, normintensity2(i,:), '-r')
end
for i=1:height(normintensity4)
    plot(time4(1:51), normintensity4(i,1:51), '-b')
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion After Incubation in PBS')
text(35, .8, 'replicate 1');
text(32, .82, '_____', 'Color', 'red');
text(35, .75, 'replicate 2');
text(32, .77, '_____', 'Color', 'blue');
saveas(gcf,'PBS_normintensity.png')
saveas(gcf, 'PBS_normintensity.fig')

%% determine whether position affects fluorescence traces

%set the color for the traces
colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
%set the position of the text in the legend
ty=[.8, .75, .7, .65, .6, .55];

figure(7), hold on
for i=1:height(normintensity1)
    for h=2:length(positions1)
        if i>positions1(h-1)&i<=positions1(h)
            plot(time1, normintensity1(i,:), 'Color', colorcode{h-1})
        end
    end
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion in LB, replicate 1')
for i=1:length(positions1)-1
    text(time1(end)-10, ty(i), ['position ' num2str(i)]);
    text(time1(end)-13, ty(i)+.02, '_____', 'Color', colorcode{i}, 'FontWeight', 'bold');
end

figure(8), hold on
for i=1:height(normintensity2)
    for h=2:length(positions2)
        if i>positions2(h-1)&i<=positions2(h)
            plot(time2, normintensity2(i,:), 'Color', colorcode{h-1})
        end
    end
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion After PBS Incubation, replicate 1')
for i=1:length(positions2)-1
    text(time2(end)-10, ty(i), ['position ' num2str(i)]);
    text(time2(end)-13, ty(i)+.02, '_____', 'Color', colorcode{i}, 'FontWeight', 'bold');
end

figure(9), hold on
for i=1:height(normintensity3)
    for h=2:length(positions3)
        if i>positions3(h-1)&i<=positions3(h)
            plot(time3, normintensity3(i,:), 'Color', colorcode{h-1})
        end
    end
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion in LB, replicate 2')
for i=1:length(positions3)-1
    text(time3(end)-10, ty(i), ['position ' num2str(i)]);
    text(time3(end)-13, ty(i)+.02, '_____', 'Color', colorcode{i}, 'FontWeight', 'bold');
end

figure(10), hold on
for i=1:height(normintensity4)
    for h=2:length(positions4)
        if i>positions4(h-1)&i<=positions4(h)
            plot(time4, normintensity4(i,:), 'Color', colorcode{h-1})
        end
    end
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion After PBS Incubation, replicate 2')
for i=1:length(positions4)-1
    text(time4(end)-10, ty(i), ['position ' num2str(i)]);
    text(time4(end)-13, ty(i)+.02, '_____', 'Color', colorcode{i}, 'FontWeight', 'bold');
end

%% fit the data to an exponential model
modelfun=@(tau,x)exp(-x./tau);
tau0=0.1;

%10232021_Exp1
fit1=[];
tau1=[];
y_hat1=normintensity1;
for i=1:height(normintensity1)
    if ~isnan(normintensity1(i,:))
        fit1=[fit1, i];
        tau1(i)=nlinfit(time1, normintensity1(i,:), modelfun, tau0);
        y_hat1(i,:)=modelfun(tau1(i), time1);   
    end
end

figure, hold on
for i=1:height(normintensity1)
    plot(time1, normintensity1(i,:),'-g',...
    time1, y_hat1(i,:), '--k')
end
legend({'data', 'model'})
title('Diffusion in LB, replicate 1')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10232021_Exp1_fit.png')
saveas(gcf, '10232021_Exp1_fit.fig')

%10232021_Exp2
fit2=[];
tau2=[];
y_hat2=normintensity2;
for i=1:height(normintensity2)
    if ~isnan(normintensity2(i,:))
        fit2=[fit2, i];
        tau2(i)=nlinfit(time2, normintensity2(i,:), modelfun, tau0);
        y_hat2(i,:)=modelfun(tau2(i), time2);   
    end
end

figure, hold on
for i=1:height(normintensity2)
    plot(time2, normintensity2(i,:),'-g',...
        time2, y_hat2(i,:), '--k')
end
legend({'data', 'model'})
title('Diffusion After Incubation in PBS, replicate 1')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10232021_Exp2_fit.png')
saveas(gcf, '10232021_Exp2_fit.fig')

%10262021_Exp1
fit3=[];
tau3=[];
y_hat3=normintensity3;
for i=1:height(normintensity3)
    if ~isnan(normintensity3(i,:))
        fit3=[fit3, i];
        tau3(i)=nlinfit(time3, normintensity3(i,:), modelfun, tau0);
        y_hat3(i,:)=modelfun(tau3(i), time3);   
    end
end

figure, hold on
for i=1:height(normintensity3)
    plot(time3, normintensity3(i,:),'-g',...
        time3, y_hat3(i,:), '--k')
end
legend({'data', 'model'})
title('Diffusion in LB, replicate 2')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10262021_Exp1_fit.png')
saveas(gcf, '10262021_Exp1_fit.fig')

%10262021_Exp2
fit4=[];
tau4=[];
y_hat4=normintensity4;
for i=1:height(normintensity4)
    if ~isnan(normintensity4(i,:))
        fit4=[fit4, i];
        tau4(i)=nlinfit(time4, normintensity4(i,:), modelfun, tau0);
        y_hat4(i,:)=modelfun(tau4(i), time4);  
    end
end

figure, hold on
for i=1:height(normintensity4)
    plot(time4, normintensity4(i,:),'-g',...
        time4, y_hat4(i,:), '--k')
end
legend({'data', 'model'})
title('Diffusion After Incubation in PBS, replicate 2')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10262021_Exp2_fit.png')
saveas(gcf, '10262021_Exp2_fit.fig')

% compare the time constants
figure, hold on
histogram(tau1, 'BinWidth', 1)
histogram(tau3, 'BinWidth', 1)
histogram(tau2, 'BinWidth', 1)
histogram(tau4, 'BinWidth', 1)
legend({'LB, rep1', 'LB, rep2', 'PBS, rep1', 'PBS, rep2'})
title('Time Constants')
xlabel('\tau (min^{-1})')
ylabel('Count')
saveas(gcf,'timescale_hist.png')
saveas(gcf, 'timescale_hist.fig')

% find the mean and standard deviation of the time constants
tau_means = [mean(tau1), mean(tau3), mean(tau2), mean(tau4)];
tau_std = [std(tau1), std(tau3), std(tau2), std(tau4)];

figure
scatter(categorical({'LB, rep1', 'LB, rep2', 'PBS, rep1', 'PBS, rep2'}), tau_means)

