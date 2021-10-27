%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

dirsave="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/10242021_analysis";
%% Make mock plot for 10/10/2021 Exp6 expected results for Group Meeting
cd(dirsave)

% x = [0:23];
% y = [repelem(0, 5), repelem(-5, 5), repelem(0,4), [-5,-3.5, -2.5, -1, 0], [0:-0.5:-2]];
% 
% figure(1)
% plot(x,y)
% xlabel('Time (minutes)')
% ylabel('C_{in}-C_{out} (A.U.)')
% title('Predicted C_{in}-C_{out} Profile')
% ylim([-10,5])
% xline(5, '--k', 'dextran')
% xline(10, '--k', 'detergent')
% xline(14, '--k', 'dextran')
% xline(19, '--k', 'dextran + Mg_{2+}')
% saveas(gcf,'10102021_Exp6_mock.png')
% saveas(gcf, '10102021_Exp6_mock.fig')

%% Perform analysis of 10/10/2021 Exp6 data
% datadir=dir('*BTfluo.mat');
% 
% diffC_total=[];
% for i=1:length(datadir)
%     load(datadir(i).name, 'diffC', 'time')
%     diffC_total=[diffC_total; diffC];
% end
% 
% figure(2), hold on
% for i=1:height(diffC_total)
%     plot(time./60, diffC_total(i, :))
% end
% xlabel('Time (minutes)')
% ylabel('C_{in}-C_{out} (A.U.)')
% title('C_{in}-C_{out} vs Time')
% xline(6, '--k', 'dextran')
% xline(10, '--k', 'detergent')
% xline(14, '--k', 'dextran')
% xline(19, '--k', 'dextran + Mg_{2+}')
% saveas(gcf,'10102021_Exp6_diffC.png')
% saveas(gcf, '10102021_Exp6_diffC.fig')

%% Perform analysis of 08/26/2021 Exp1 and 10/23/2021 Exp1 
datadir1=dir(['08262021' '*dm.mat']);
datadir2=dir(['10232021_Exp1' '*dm2.mat']);
datadir3=dir(['09092021_Exp2' '*dm.mat']);
datadir4=dir(['10232021_Exp2' '*dm2.mat']);

intensity1=[];
normintensity1=[];
for i=1:length(datadir1)
    load(datadir1(i).name, 'icell_green', 'norm_green', 'time')
    intensity1=[intensity1; icell_green];
    normintensity1=[normintensity1; norm_green];
end
time1=time;

intensity2=[];
normintensity2=[];
for i=1:length(datadir2)
    load(datadir2(i).name, 'icell_intensity', 'norm_intensity', 'time')
    intensity2=[intensity2; icell_intensity];
    normintensity2=[normintensity2; norm_intensity];
end
time2=time;

intensity3=[];
normintensity3=[];
for i=1:length(datadir3)
    load(datadir3(i).name, 'icell_intensity', 'norm_intensity', 'time')
    intensity3=[intensity3; icell_intensity];
    normintensity3=[normintensity3; norm_intensity];
end
time3=time;

intensity4=[];
normintensity4=[];
for i=1:length(datadir4)
    load(datadir4(i).name, 'icell_intensity', 'norm_intensity', 'time')
    intensity4=[intensity4; icell_intensity];
    normintensity4=[normintensity4; norm_intensity];
end
time4=time;

figure(3), hold on
for i=1:height(intensity1)
    plot(time1, intensity1(i,:))
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
title('Diffusion in LB Without IPTG')
saveas(gcf,'08262021_Exp1_intensity.png')
saveas(gcf, '08262021_Exp1_intensity.fig')

figure(4), hold on
for i=1:height(intensity2)
    plot(time2, intensity2(i,:))
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
title('Diffusion in LB With IPTG')
saveas(gcf,'10232021_Exp1_intensity.png')
saveas(gcf, '10232021_Exp1_intensity.fig')

figure(5), hold on
for i=1:height(intensity3)
    plot(time3, intensity3(i,:))
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
title('Diffusion in LB After 16-minute PBS Incubation')
saveas(gcf,'09092021_Exp2_intensity.png')
saveas(gcf, '09092021_Exp2_intensity.fig')

figure(6), hold on
for i=1:height(intensity4)
    plot(time4, intensity4(i,:))
end
xlabel('Time (minutes)')
ylabel('Fluorescence (A.U.)')
title('Diffusion in LB After 1-hour PBS Incubation')
saveas(gcf,'10232021_Exp2_intensity.png')
saveas(gcf, '10232021_Exp2_intensity.fig')

figure(7), hold on
for i=1:height(normintensity1)
    plot(time1, normintensity1(i,:))
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion in LB Without IPTG')
saveas(gcf,'08262021_Exp1_norm.png')
saveas(gcf, '08262021_Exp1_norm.fig')

figure(8), hold on
for i=1:height(normintensity2)
    plot(time2, normintensity2(i,:))
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion in LB With IPTG')
saveas(gcf,'10232021_Exp1_norm.png')
saveas(gcf, '10232021_Exp1_norm.fig')

figure(9), hold on
for i=1:height(normintensity3)
    plot(time3, normintensity3(i,:))
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion in LB After 16-minute PBS Incubation')
saveas(gcf,'09092021_Exp2_norm.png')
saveas(gcf, '09092021_Exp2_norm.fig')

figure(10), hold on
for i=1:height(normintensity4)
    plot(time4, normintensity4(i,:))
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion in LB After 1-hour PBS Incubation')
saveas(gcf,'10232021_Exp2_norm.png')
saveas(gcf, '10232021_Exp2_norm.fig')

%% Modeling
beta=0.25;
modelfun=@(coeff,x)(1-beta)*exp(-x./coeff)+beta;
coeff0=1;

%08262021_Exp1
coeff1=[];
y_hat1=normintensity1;
for i=1:height(normintensity1)
    coeff1(i)=nlinfit(time1, normintensity1(i,:), modelfun, coeff0);
    y_hat1(i,:)=modelfun(coeff1(i), time1);   
end

figure, hold on
for i=1:height(normintensity1)
    plot(time1, normintensity1(i,:),'-g',...
        time1, y_hat1(i,:), '--k')
end
legend({'data', 'model'})
title('Diffusion in LB Without IPTG')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'08262021_Exp1_fit.png')
saveas(gcf, '08262021_Exp1_fit.fig')

%10232021_Exp1
coeff2=[];
y_hat2=normintensity2;
for i=1:height(normintensity2)
    if isnan(normintensity2(i,:))==0
    coeff2(i)=nlinfit(time2, normintensity2(i,:), modelfun, coeff0);
    y_hat2(i,:)=modelfun(coeff2(i), time2);   
    else
    end
end

figure, hold on
for i=1:height(normintensity2)
    plot(time2, normintensity2(i,:),'-g',...
        time2, y_hat2(i,:), '--k')
end
legend({'data', 'model'})
title('Diffusion in LB With IPTG')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10232021_Exp1_fit.png')
saveas(gcf, '10232021_Exp1_fit.fig')

%09092021_Exp2
coeff3=[];
y_hat3=normintensity3;
for i=1:height(normintensity3)
    coeff3(i)=nlinfit(time3, normintensity3(i,:), modelfun, coeff0);
    y_hat3(i,:)=modelfun(coeff3(i), time3);   
end

figure, hold on
for i=1:height(normintensity3)
    plot(time3, normintensity3(i,:),'-g',...
        time3, y_hat3(i,:), '--k')
end
legend({'data', 'model'})
title('Diffusion in LB After 16-minute PBS Incubation')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'09092021_Exp2_fit.png')
saveas(gcf, '09092021_Exp2_fit.fig')

%10232021_Exp2
coeff4=[];
y_hat4=normintensity4;
for i=1:height(normintensity4)
    coeff4(i)=nlinfit(time4, normintensity4(i,:), modelfun, coeff0);
    y_hat4(i,:)=modelfun(coeff4(i), time4);   
end

figure, hold on
for i=1:height(normintensity4)
    plot(time4, normintensity4(i,:),'-g',...
        time4, y_hat4(i,:), '--k')
end
legend({'data', 'model'})
title('Diffusion in LB After 1-hour PBS Incubation')
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
saveas(gcf,'10232021_Exp2_fit.png')
saveas(gcf, '10232021_Exp2_fit.fig')

figure, hold on
histogram(coeff1, 'BinWidth', 5)
histogram(coeff2, 'BinWidth', 5)
histogram(coeff3(coeff3<500), 'BinWidth', 5)
histogram(coeff4, 'BinWidth', 5)
legend({'LB without induction', 'LB with induction', '16-minute PBS incubation', '1-hour PBS incubation'})
title('Time Constants')
xlabel('\tau (min^{-1})')
ylabel('Count')
saveas(gcf,'timescale_hist.png')
saveas(gcf, 'timescale_hist.fig')

figure, hold on
histogram(coeff1, 'BinWidth', 5)
histogram(coeff2, 'BinWidth', 5)
legend({'LB without induction', 'LB with induction'})
title('Time Constants')
xlabel('\tau (min^{-1})')
ylabel('Count')
saveas(gcf,'timescale_hist_LB.png')
saveas(gcf, 'timescale_hist_LB.fig')

%% Check
cellLength3=[];
for i=1:length(datadir3)
    load(datadir3(i).name, 'lcell', 'T')
    cellLength3=[cellLength3; lcell];
end
timeL3=[0:T-1];

figure(1), hold on
for i=1:height(cellLength3)
    plot(timeL3, cellLength3(i,:))
end

figure(2), hold on
for i=1:height(normintensity3)
    plot(time3, normintensity3(i,:))
end
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence (A.U.)')
title('Diffusion in LB After 16-minute PBS Incubation')