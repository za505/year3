%Author: Zarina Akbary
%Date: 03/03/2022
%Purpose: To model mNeonGreen diffusion out of membrane-lysed cells

clear, close all

%directory guide
%untreated, 1 minute frame rate, 100% intensity = '10232021_Exp1' & '10262021_Exp1' 
%untreated, 1 minute frame rate, 20% intensity = '02212022_Exp2' 
%untreated, 5 minute frame rate, 100% intensity = '02122022_Exp1' 
%untreated, 10 minute frame rate, 100% intensity = '02122022_Exp2'
%untreated, 20 minute frame rate, 100% intensity = '02092022_Exp1'
%untreated, 20 minute frame rate, 20% intensity = '02212022_Exp1'
 
%2 minute PBS incubation, 1 minute frame rate, 100% intensity = '11192021_Exp2' & '12082021_Exp3' 
%20 minute PBS incubation, 1 minute frame rate, 100% intensity = '11192021_Exp1' & '11302021_Exp1'
%60 minute PBS incubation, 1 minute frame rate, 100% intensity = '10232021_Exp2' & '10262021_Exp2'
%120 minute PBS incubation, 1 minute frame rate, 100% intensity = '01142022_Exp1'

% LBMga=dir(['01172022_Exp1' '*dm.mat']); %LB, 20 mM Mg2+
% LBEa=dir(['01172022_Exp2' '*dm.mat']); %LB, 10 mM EDTA
% LBtuna=dir(['01242022_Exp1' '*dm.mat']); %LB, 0.5 ug/mL tunicamycin
% LBvana=dir(['01262022_Exp1' '*dm.mat']); %LB, 1 ug/mL vancomycin
% LBspnta=dir(['02122022_Exp3' '*dm.mat']); %LB, spent media (10 min frame rate)
% LBspntb=dir(['02192022_Exp1' '*dm.mat']); %LB, spent media (10 min frame rate)
% 
%    
% LB1s=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
% LB2s=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
% LB3s=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s
% 
% LB1t=dir(['02192022_Exp2' '*dm.mat']); %LB, frame rate = 1.76 s
% LB2t=dir(['02192022_Exp3' '*dm.mat']); %LB, frame rate = 2.3 s
% LB3t=dir(['02192022_Exp4' '*dm.mat']); %LB, frame rate = 3 s

%for photobleach correction (100% intensity)
% alpha=32.2114;
% intercept=0.1614;

%for photobleach correction (20% intensity)
% alpha=151.7544;
% intercept=-1.0566;

%% User Input 

%color codes must be in RGB [0-1] format to be used in ciplot
colorcode={[204 0 0], [204 102 0], [204 204 0], [102 204 0], [0 204 204], [0 0 204], [102 0 204], [204 0 204], [204 0 102], [255 102 102], [255 178 102], [102 255 102], [102 255 255], [102 178 255], [178 102 255], [255 102 255],[255 102 178]};
colorcode2={[255 51 51], [255 153 51], [255 255 51], [153 255 51], [51 255 255], [51 51 255], [153 51 255], [255 51 255], [255 51 153], [255 204 204], [255 229 204], [204 255 204], [204 255 255], [204 229 255], [229 204 255], [255 204 255],[255 204 229]};

colorcode=cellfun(@(x)(x./255), colorcode, 'UniformOutput', false);
colorcode2=cellfun(@(x)(x./255), colorcode2, 'UniformOutput', false);

m=10;
c1=[0.847 0.427 0.933];
c2=[0.101 0.101 0.737];
color_p=[linspace(c1(1), c2(1), m)', linspace(c1(2), c2(2), m)', linspace(c1(3), c2(3), m)'];

transparency = 0.3; %this is the alpha argument for the ciplot function

%location and names of of the *.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis';
basenames={'10232021_Exp1', '10262021_Exp1', '02212022_Exp2', '02122022_Exp1','02122022_Exp2', '02092022_Exp1', '02212022_Exp1', '11192021_Exp2', '12082021_Exp3', '11192021_Exp1', '11302021_Exp1', '10232021_Exp2', '10262021_Exp2', '01142022_Exp1', '01172022_Exp1', '01172022_Exp2', '01242022_Exp1', '01262022_Exp1', '02122022_Exp3', '02192022_Exp1','11202021_Exp1', '12082021_Exp1', '12082021_Exp2', '02192022_Exp2', '02192022_Exp3', '02192022_Exp4'};

%% normalize data
imstarts=[6, 5, 6, 6, 6, 6, 6, 10, 10, 4, 4, 3, 5, 4, 6, 6, 3, 5, 6, 6, NaN, NaN, NaN, NaN, NaN, NaN]; %based on time vector
%idxes={[2:5], [2:6], [2:3], [2:3], [2:3], [2:3], [2:3], [2:6], [2:6], [2:6], [2:4], [2:4], [2:4], [], [2:7], [2:4], [2:4], [2:3], [2:3], [], [], [], [], [], [], []}; %based on tme vector
idxes={[2:4], [2:6], [2:4], [2:3], [2:3], [2:3], [2:3], [2:4], [2:5], [2:4], [2:4], [2:4], [2:4], [], [2:7], [2:4], [2:4], [], [2:3], [2:3], [], [], [], [], [], []}; %based on tme vector

for b=1:length(basenames)

    basename=basenames{b};
    intp=idxes{b};
    
    cd([dirsave '/rawFiles'])
    
    datadir=dir([basename '*']);
    imstart=imstarts(b);
    
    if ismember(b, [21:26]) %controls
        
        intensity=[];
        bgintensity=[];
        adjintensity=[];
        normintensity=[];
        lcell=[];
        
        for j=1:height(datadir)
            cd([dirsave '/rawFiles']);
            load(datadir(j).name);
            
            
            [intensity1, bgintensity1, adjintensity1, normintensity1, lCell, time, tme, imstart]=controlsNormalize(datadir);
            intensity=[intensity;intensity1];
            bgintensity=[bgintensity; bgintensity1];
            adjintensity=[adjintensity; adjintensity1];
            normintensity=[normintensity; normintensity1];
            lcell=[lcell; lCell];
            
        end
        
        cd([dirsave '/normalizedFiles'])
        save([basename '_norm.mat'], 'intensity', 'bgintensity', 'adjintensity', 'normintensity', 'lcell', 'imstart', 'time', 'tme');
    
    else
        
        intensity=[];
        adjintensity=[];
        normintensity=[];
        lcell=[];
        
        for j=1:height(datadir)
            cd([dirsave '/rawFiles'])
            load(datadir(j).name);
            
            
            [intensity1, adjintensity1, normintensity1, lCell, time, tme, imstart]=dataNormalize(datadir, imstart, intp);
            intensity=[intensity;intensity1];
            adjintensity=[adjintensity; adjintensity1];
            normintensity=[normintensity; normintensity1];
            lcell=[lcell; lCell];
            
        end
        
        cd([dirsave '/normalizedFiles'])
        save([basename '_norm.mat'], 'intensity', 'adjintensity', 'normintensity', 'lcell', 'imstart', 'time', 'tme', 'intp');
    end
    
end


%% correct data
for i=1:length(basenames)
    
    basename=basenames{i};
    cd([dirsave '/normalizedFiles'])
    datadir=dir([basename '*']);
    
    if ismember(i, [3, 7, 24:26]) % 20% intensity
        
        %for photobleach correction (20% intensity)
%         alpha=152.6455;
%         intercept=-1.0976;

        alpha=151.7544;
        intercept=-1.0566;
        load(datadir.name)
        [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept);

        cd([dirsave '/correctedFiles'])
        save([basename '_corrected.mat'], 'time', 'tme', 'alpha', 'intercept', 'Cnew', 'dCB', 'dCT', 'dCP', 'Cbl_exp', 'unb_frac');
    
    else
        
        %for photobleach correction (100% intensity)
        alpha=32.2114;
        intercept=0.1614;
        
        load(datadir.name);
        [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept);

        cd([dirsave '/correctedFiles'])
        save([basename '_corrected.mat'], 'time', 'tme', 'alpha', 'intercept', 'Cnew', 'dCB', 'dCT', 'dCP', 'Cbl_exp', 'unb_frac');
    end
    
end

%% compare the controls traces
%each row is a different frame rate (1, 2, and 3 second), each column is a
%different variable (time, tme, intensity, bgintensity, adjintensity,
%normintensity, lcell, Cnew, imstart, tau, and yhat)

controls_100 = cell(3, 11);
controls_20 = cell(3, 11);

idx1=1;
idx2=1;

for i=1:length(basenames)
    basename=basenames{i};
    
    if ismember(i, [21:23])
        cd([dirsave '/normalizedFiles'])
        load([basename '_norm.mat'])
        controls_100{idx1, 1}=time;
        controls_100{idx1, 2}=tme;
        controls_100{idx1, 3}=intensity;
        controls_100{idx1, 4}=bgintensity;
        controls_100{idx1, 5}=adjintensity;
        controls_100{idx1, 6}=normintensity;
        controls_100{idx1, 7}=lcell;
        controls_100{idx1, 9}=imstart;
        
        [controls_100{idx1, 10}, controls_100{idx1, 11}]=tauCalc(tme, normintensity, 1);
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        controls_100{idx1, 8}=Cnew;
        
        idx1=idx1+1;
        
    elseif ismember(i, [24:26])
        
        cd([dirsave '/normalizedFiles'])
        load([basename '_norm.mat'])
        controls_20{idx2, 1}=time;
        controls_20{idx2, 2}=tme;
        controls_20{idx2, 3}=intensity;
        controls_20{idx2, 4}=bgintensity;
        controls_20{idx2, 5}=adjintensity;
        controls_20{idx2, 6}=normintensity;
        controls_20{idx2, 7}=lcell;
        controls_20{idx2, 9}=imstart;
        
        [controls_20{idx2, 10}, controls_20{idx2, 11}]=tauCalc(tme, normintensity, 1);
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        controls_20{idx2, 8}=Cnew;
        
        idx2=idx2+1;
    else
        continue
    end
    
end

%% compare the untreated traces
%each row is a different frame rate (1 min, 1 min, 5 min, 10 min, and 20 min for 100% intensity and 1 min and 20 min for 20% intensity), each column is a
%different variable (time, tme, intensity, adjintensity, normintensity, lcell, and Cnew)

untreated_100 = cell(5, 7);
untreated_20 = cell(2, 7);

idx1=1;
idx2=1;

for i=1:length(basenames)
    basename=basenames{i};
    
     if ismember(i, [1,2,4,5,6])
        cd([dirsave '/normalizedFiles'])
        load([basename '_norm.mat'])
        untreated_100{idx1, 1}=time;
        untreated_100{idx1, 2}=tme;
        untreated_100{idx1, 3}=intensity;
        untreated_100{idx1, 4}=adjintensity;
        untreated_100{idx1, 5}=normintensity;
        untreated_100{idx1, 6}=lcell;
        untreated_100{idx1, 8}=imstart;
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        untreated_100{idx1, 7}=Cnew;
        
        idx1=idx1+1;
        
    elseif ismember(i, [3,7])
        cd([dirsave '/normalizedFiles'])
        load([basename '_norm.mat'])
        untreated_20{idx2, 1}=time;
        untreated_20{idx2, 2}=tme;
        untreated_20{idx2, 3}=intensity;
        untreated_20{idx2, 4}=adjintensity;
        untreated_20{idx2, 5}=normintensity;
        untreated_20{idx2, 6}=lcell;
        untreated_20{idx2, 8}=imstart;
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        untreated_20{idx2, 7}=Cnew;
        
        idx2=idx2+1;
    else
        continue
    end
    
end

%% compare the PBS traces
%each row is a different frame rate (1 min, 1 min, 5 min, 10 min, and 20 min for 100% intensity and 1 min and 20 min for 20% intensity), each column is a
%different variable (time, tme, intensity, adjintensity, normintensity, lcell, and Cnew)

PBS_100 = cell(7, 7);

idx1=1;

for i=1:length(basenames)
    basename=basenames{i};
    
     if ismember(i, [8:14])
        cd([dirsave '/normalizedFiles'])
        load([basename '_norm.mat'])
        PBS_100{idx1, 1}=time;
        PBS_100{idx1, 2}=tme;
        PBS_100{idx1, 3}=intensity;
        PBS_100{idx1, 4}=adjintensity;
        PBS_100{idx1, 5}=normintensity;
        PBS_100{idx1, 6}=lcell;
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        PBS_100{idx1, 7}=Cnew;
        
        idx1=idx1+1;
        
     end
    
end

%% compare the LB treated traces
%each row is a different treatment in LB (100% intensity)

treated_100 = cell(6, 7);

idx1=1;

for i=1:length(basenames)
    basename=basenames{i};
    
     if ismember(i, [15:20])
        cd([dirsave '/normalizedFiles'])
        load([basename '_norm.mat'])
        treated_100{idx1, 1}=time;
        treated_100{idx1, 2}=tme;
        treated_100{idx1, 3}=intensity;
        treated_100{idx1, 4}=adjintensity;
        treated_100{idx1, 5}=normintensity;
        treated_100{idx1, 6}=lcell;
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        treated_100{idx1, 7}=Cnew;
        
        idx1=idx1+1;
        
     end
    
end

%% plot the fit of the exponential on the control traces
cd([dirsave '/03072022_groupMeeting']);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',16)
subplot(2,3,1)
plot(controls_100{1,2}, controls_100{1, 6}, 'Color', colorcode{1}), hold on
plot(controls_100{1,2}, controls_100{1, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('100% Intensity, 1 s Frame Rate')
% saveas(gcf, 'fit100_1s.png')
% saveas(gcf, 'fit100_1s.fig')

subplot(2,3,2)
plot(controls_100{2,2}, controls_100{2, 6}, 'Color', colorcode{3}), hold on
plot(controls_100{2,2}, controls_100{2, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('100% Intensity, 2 s Frame Rate')
% saveas(gcf, 'fit100_2s.png')
% saveas(gcf, 'fit100_2s.fig')

subplot(2,3,3)
plot(controls_100{3,2}, controls_100{3, 6}, 'Color', colorcode{5}), hold on
plot(controls_100{3,2}, controls_100{3, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('100% Intensity, 3 s Frame Rate')
% saveas(gcf, 'fit100_3s.png')
% saveas(gcf, 'fit100_3s.fig')

subplot(2,3,4)
plot(controls_20{1,2}, controls_20{1, 6}, 'Color', colorcode{2}), hold on
plot(controls_20{1,2}, controls_20{1, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('20% Intensity, 1 s Frame Rate')
% saveas(gcf, 'fit20_1s.png')
% saveas(gcf, 'fit20_1s.fig')

subplot(2,3,5)
plot(controls_20{2,2}, controls_20{2, 6}, 'Color', colorcode{4}), hold on
plot(controls_20{2,2}, controls_20{2, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('20% Intensity, 2 s Frame Rate')
% saveas(gcf, 'fit20_2s.png')
% saveas(gcf, 'fit20_2s.fig')

subplot(2,3,6)
plot(controls_20{3,2}, controls_20{3, 6}, 'Color', colorcode{6}), hold on
plot(controls_20{3,2}, controls_20{3, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('20% Intensity, 3 s Frame Rate')
% saveas(gcf, 'fit20_3s.png')
% saveas(gcf, 'fit20_3s.fig')

saveas(gcf, 'fitControls.png')
saveas(gcf, 'fitControls.fig')
%% calculate the slope of tau as a function of frame rate
tau={controls_100{1, 10}; controls_100{2, 10}; controls_100{3, 10}; controls_20{1, 10}; controls_20{2, 10}; controls_20{3, 10}};
dt=[1.2/60, 2/60, 3/60, 1.76/60, 2.3/60, 3/60]; 

linearCoef1 = polyfit([repelem(dt(1), length(tau{1})), repelem(dt(2), length(tau{2})), repelem(dt(3), length(tau{3}))],[tau{1}', tau{2}', tau{3}'],1);
linearFit1= polyval(linearCoef1,[0 1.2/60 2/60 3/60]);

linearCoef2 = polyfit([repelem(dt(4), length(tau{4})), repelem(dt(5), length(tau{5})), repelem(dt(6), length(tau{6}))],[tau{4}', tau{5}', tau{6}'],1);
linearFit2= polyval(linearCoef2,[0 1.76/60 2.3/60 3/60]);

%% plot the slope as a function of frame rate
tau_means = cellfun(@(x)mean(x, 1, 'omitnan'), tau);
tau_std = cellfun(@(x)std(x, 0, 1, 'omitnan'), tau);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
scatter(repelem(dt(1), length(tau{1})), tau{1}, 'MarkerFaceColor', '#1267A1', 'MarkerEdgeColor', '#1267A1')
scatter(repelem(dt(2), length(tau{2})), tau{2}, 'MarkerFaceColor', '#1267A1', 'MarkerEdgeColor', '#1267A1')
scatter(repelem(dt(3), length(tau{3})), tau{3}, 'MarkerFaceColor', '#1267A1', 'MarkerEdgeColor', '#1267A1')

scatter(repelem(dt(4), length(tau{4})), tau{4}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(dt(5), length(tau{5})), tau{5}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(dt(6), length(tau{6})), tau{6}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')

scatter([dt(1), dt(2), dt(3), dt(4), dt(5), dt(6)], tau_means, 'MarkerFaceColor', 'black')
errorbar(dt(1), tau_means(1), tau_std(1), 'Color', 'black')
errorbar(dt(2), tau_means(2), tau_std(2), 'Color', 'black')
errorbar(dt(3), tau_means(3), tau_std(3), 'Color', 'black')
errorbar(dt(4), tau_means(4), tau_std(4), 'Color', 'black')
errorbar(dt(5), tau_means(5), tau_std(5), 'Color', 'black')
errorbar(dt(6), tau_means(6), tau_std(6), 'Color', 'black')

plot([0, dt(1), dt(2), dt(3)], linearFit1, '--k', 'LineWidth', 1)
plot([0, dt(4), dt(5), dt(6)], linearFit2, '--k', 'LineWidth', 1)
xlim([0, 0.06])
ylabel('\tau (min^{-1})')
%xticks([0, dt(1),  dt(4), dt(2),  dt(5), dt(3)])
%xticklabels({'0', '1.2 s', '1.76 s', '2 s', '2.3 s', '3 s'})
xlabel('Time (minutes)')

caption1 = sprintf('%f * frame rate + %f', linearCoef1(1), linearCoef1(2));
caption2 = sprintf('%f * frame rate + %f', linearCoef2(1), linearCoef2(2));
text(dt(2), tau_means(2)+1.5, ['\tau = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
text(dt(1), tau_means(5)+1.5, ['\tau = ' caption2], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

title('Tau vs Frame Rate')
saveas(gcf, 'tau_vs_frameRate.png')
saveas(gcf, 'tau_vs_frameRate.fig')

%% calculate the asymptote value as a function of frame rate
%untreated 100% intensity with frame rates at: 1 s, 2 s, 3 s, 1 min, 5 min,
%10 min and 20 min (although exclude the 20 minute if need be since it
%doesn't quite asymptote)

asymptote=cell(7, 2); %column 1 = dt, column 2 = final raw fluor values

dt2=[1.2/60, 2/60, 3/60, 1, 5, 10, 20];

asymptote{1,2}=controls_100{1,3}(:, end);
asymptote{2,2}=controls_100{2,3}(:, end);
asymptote{3,2}=controls_100{3,3}(:, end);
asymptote{4,2}=untreated_100{1,3}(:, end);
asymptote{5,2}=untreated_100{3,3}(:, end);
asymptote{6,2}=untreated_100{4,3}(:, end);
asymptote{7,2}=untreated_100{5,3}(:, end);

for i=1:length(dt2)
    asymptote{i,1}=repelem(dt2(i), height(asymptote{i,2}))';
end

% asymptote_x=cell2mat(asymptote(:, 1))';
% asymptote_y=cell2mat(asymptote(:, 2))';
% linearCoef3 = polyfit(asymptote_x, asymptote_y, 1);
% linearFit3 = polyval(linearCoef3,[0 dt2]);

%% generate plots to illustrate the asymptote value as a function of frame rate
asymptote_means = cellfun(@(x)mean(x, 1, 'omitnan'), asymptote(:, 2));
asymptote_std = cellfun(@(x)std(x, 0, 1, 'omitnan'), asymptote(:, 2));

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
scatter(asymptote{1, 1}, asymptote{1, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(asymptote{2, 1}, asymptote{2, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(asymptote{3, 1}, asymptote{3, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(asymptote{4, 1}, asymptote{4, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(asymptote{5, 1}, asymptote{5, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(asymptote{6, 1}, asymptote{6, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(asymptote{7, 1}, asymptote{7, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')

scatter([dt2(1), dt2(2), dt2(3), dt2(4), dt2(5), dt2(6), dt2(7)], asymptote_means, 'MarkerFaceColor', 'black')
errorbar(dt2(1), asymptote_means(1), asymptote_std(1), 'Color', 'black')
errorbar(dt2(2), asymptote_means(2), asymptote_std(2), 'Color', 'black')
errorbar(dt2(3), asymptote_means(3), asymptote_std(3), 'Color', 'black')
errorbar(dt2(4), asymptote_means(4), asymptote_std(4), 'Color', 'black')
errorbar(dt2(5), asymptote_means(5), asymptote_std(5), 'Color', 'black')
errorbar(dt2(6), asymptote_means(6), asymptote_std(6), 'Color', 'black')
errorbar(dt2(7), asymptote_means(7), asymptote_std(7), 'Color', 'black')

plot([0 dt2], linearFit3, '--k', 'LineWidth', 1)
ylabel('Asymptote Value (A.U.)')
xlabel('Time (minutes)')
ylim([0 1800])
xlim([-1 dt2(end)+1])

% caption3 = sprintf('%f * frame rate + %f', linearCoef3(1), linearCoef3(2));
% text(dt2(2), asymptote_means(2)+1.5, ['\asymptote = ' caption3], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

title('Asymptote Value vs Frame Rate')
saveas(gcf, 'asymptote_vs_frameRate.png')
saveas(gcf, 'asymptote_vs_frameRate.fig')
%% generate plots for the controls
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
p1=plot(controls_100{1,2}, controls_100{1,3}(:, controls_100{1,9}:end), 'Color', colorcode{1}, 'LineWidth', 0.75)
p2=plot(controls_100{2,2}, controls_100{2,3}(:, controls_100{2,9}:end), 'Color', colorcode{3}, 'LineWidth', 0.75)
p3=plot(controls_100{3,2}, controls_100{3,3}(:, controls_100{3,9}:end), 'Color', colorcode{5}, 'LineWidth', 0.75)
p4=plot(controls_100{1,2}, controls_100{1,4}(:, controls_100{1,9}:end), '--', 'Color', colorcode{1}, 'LineWidth', 1.5)
p5=plot(controls_100{2,2}, controls_100{2,4}(:, controls_100{2,9}:end), '--', 'Color', colorcode{3}, 'LineWidth', 1.5)
p6=plot(controls_100{3,2}, controls_100{3,4}(:, controls_100{3,9}:end), '--', 'Color', colorcode{5}, 'LineWidth', 1.5)
ylim([0 Inf])
ylabel('Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'rawControls_100.png')
saveas(gcf, 'rawControls_100.fig')
%set(gca,'DefaultTextFontSize', 42)
% hleg=legend([p4, p5, p6], {'1 s', '2 s', '3 s'})
% title(hleg,'Frame Rate')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
p1=plot(controls_20{1,2}, controls_20{1,3}(:, controls_20{1,9}:end), 'Color', colorcode{2}, 'LineWidth', 0.75)
p2=plot(controls_20{2,2}, controls_20{2,3}(:, controls_20{2,9}:end), 'Color', colorcode{4}, 'LineWidth', 0.75)
p3=plot(controls_20{3,2}, controls_20{3,3}(:, controls_20{3,9}:end), 'Color', colorcode{6}, 'LineWidth', 0.75)
p4=plot(controls_20{1,2}, controls_20{1,4}(:, controls_20{1,9}:end), '--', 'Color', colorcode{2}, 'LineWidth', 1.5)
p5=plot(controls_20{2,2}, controls_20{2,4}(:, controls_20{2,9}:end), '--', 'Color', colorcode{4}, 'LineWidth', 1.5)
p6=plot(controls_20{3,2}, controls_20{3,4}(:, controls_20{3,9}:end), '--', 'Color', colorcode{6}, 'LineWidth', 1.5)
ylim([0 Inf])
ylabel('Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'rawControls_20.png')
saveas(gcf, 'rawControls_20.fig')
%set(gca,'DefaultTextFontSize', 42)
% hleg=legend([p4, p5, p6], {'1 s', '2 s', '3 s'})
% title(hleg,'Frame Rate')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
meanPlot(controls_100{1,2}, controls_100{1, 6}, colorcode{1}, colorcode2{1})
meanPlot(controls_100{2,2}, controls_100{2, 6}, colorcode{3}, colorcode2{3})
meanPlot(controls_100{3,2}, controls_100{3, 6}, colorcode{5}, colorcode2{5})
meanPlot(controls_20{1,2}, controls_20{1, 6}, colorcode{2}, colorcode2{2})
meanPlot(controls_20{2,2}, controls_20{2, 6}, colorcode{4}, colorcode2{4})
meanPlot(controls_20{3,2}, controls_20{3, 6}, colorcode{6}, colorcode2{6})
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'normControls.png')
saveas(gcf, 'normControls.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
meanPlot(controls_100{1,2}, controls_100{1, 8}, colorcode{1}, colorcode2{1}, transparency)
meanPlot(controls_100{2,2}, controls_100{2, 8}, colorcode{3}, colorcode2{3}, transparency)
meanPlot(controls_100{3,2}, controls_100{3, 8}, colorcode{5}, colorcode2{5}, transparency)
meanPlot(controls_20{1,2}, controls_20{1, 8}, colorcode{2}, colorcode2{2}, transparency)
meanPlot(controls_20{2,2}, controls_20{2, 8}, colorcode{4}, colorcode2{4}, transparency)
meanPlot(controls_20{3,2}, controls_20{3, 8}, colorcode{6}, colorcode2{6}, transparency)
ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'correctedControls.png')
saveas(gcf, 'correctedControls.fig')

%% generate plots for PBS
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{1,1}, untreated_100{1, 3}, colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{2,1}, untreated_100{2, 3}, colorcode{2}, colorcode2{2}, transparency)
meanPlot(PBS_100{1,1}, PBS_100{1, 3}, colorcode{3}, colorcode2{3}, transparency)
meanPlot(PBS_100{2,1}, PBS_100{2, 3}, colorcode{4}, colorcode2{4}, transparency)
meanPlot(PBS_100{3,1}, PBS_100{3, 3}, colorcode{5}, colorcode2{5}, transparency)
meanPlot(PBS_100{4,1}, PBS_100{4, 3}, colorcode{6}, colorcode2{6}, transparency)
meanPlot(PBS_100{5,1}, PBS_100{5, 3}, colorcode{7}, colorcode2{7}, transparency)
meanPlot(PBS_100{6,1}, PBS_100{6, 3}, colorcode{8}, colorcode2{8}, transparency)
meanPlot(PBS_100{7,1}, PBS_100{7, 3}, colorcode{9}, colorcode2{9}, transparency)
ylim([0 Inf])
ylabel('Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'rawPBS.png')
saveas(gcf, 'rawPBS.fig')

%% Generate plots for untreated LB
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{1,1}, untreated_100{1, 3}, colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{2,1}, untreated_100{2, 3}, colorcode{2}, colorcode2{2}, transparency)
meanPlot(untreated_100{3,1}, untreated_100{3, 3}, colorcode{3}, colorcode2{3}, transparency)
meanPlot(untreated_100{4,1}, untreated_100{4, 3}, colorcode{4}, colorcode2{4}, transparency)
meanPlot(untreated_100{5,1}, untreated_100{5, 3}, colorcode{5}, colorcode2{5}, transparency)
ylim([0 Inf])
ylabel('Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'rawUntreated.png')
saveas(gcf, 'rawUntreated.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{1,1}, untreated_100{1, 4}, colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{2,1}, untreated_100{2, 4}, colorcode{2}, colorcode2{2}, transparency)
meanPlot(untreated_100{3,1}, untreated_100{3, 4}, colorcode{3}, colorcode2{3}, transparency)
meanPlot(untreated_100{4,1}, untreated_100{4, 4}, colorcode{4}, colorcode2{4}, transparency)
meanPlot(untreated_100{5,1}, (untreated_100{5, 3}-mean(untreated_100{4, 3}(:,end), 1, 'omitnan')), colorcode{5}, colorcode2{5}, transparency)
ylim([0 Inf])
ylabel('Adjusted Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'adjUntreated.png')
saveas(gcf, 'adjUntreated.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{1,2}, untreated_100{1, 5}, colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{2,2}, untreated_100{2, 5}, colorcode{2}, colorcode2{2}, transparency)
meanPlot(untreated_100{3,2}, untreated_100{3, 5}, colorcode{3}, colorcode2{3}, transparency)
meanPlot(untreated_100{4,2}, untreated_100{4, 5}, colorcode{4}, colorcode2{4}, transparency)
meanPlot(untreated_100{5,2}, untreated_100{5, 5}, colorcode{5}, colorcode2{5}, transparency)
ylim([0 Inf])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'normUntreated.png')
saveas(gcf, 'normUntreated.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
meanPlot(untreated_100{2,2}, [untreated_100{1, 7}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{3,2}, untreated_100{3, 7}, colorcode{3}, colorcode2{3}, transparency)
meanPlot(untreated_100{4,2}, untreated_100{4, 7}, colorcode{4}, colorcode2{4}, transparency)
meanPlot(untreated_100{5,2}, untreated_100{5, 7}, colorcode{5}, colorcode2{5}, transparency)
ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'correctedUntreated.png')
saveas(gcf, 'correctedUntreated.fig')

%% when do the untreated cells hit the asymptote?
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
dervPlot(untreated_100{1,1}, untreated_100{1, 3}, colorcode{1})
dervPlot(untreated_100{2,1}, untreated_100{2, 3}, colorcode{2})
dervPlot(untreated_100{3,1}, untreated_100{3, 3}, colorcode{3})
dervPlot(untreated_100{4,1}, untreated_100{4, 3}, colorcode{4})
dervPlot(untreated_100{5,1}, untreated_100{5, 3}, colorcode{5})
%ylim([0 Inf])
ylabel('Relative Change in Fluorescence (A.U.)')
xlabel('Time (minutes)')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
dervPlot(untreated_100{1,2}, untreated_100{1, 5}, colorcode{1})
dervPlot(untreated_100{2,2}, untreated_100{2, 5}, colorcode{2})
dervPlot(untreated_100{3,2}, untreated_100{3, 5}, colorcode{3})
dervPlot(untreated_100{4,2}, untreated_100{4, 5}, colorcode{4})
dervPlot(untreated_100{5,2}, untreated_100{5, 5}, colorcode{5})
%ylim([0 Inf])
ylabel('Relative Change in Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
dervPlot(untreated_100{1,2}, untreated_100{1, 7}, colorcode{1})
dervPlot(untreated_100{2,2}, untreated_100{2, 7}, colorcode{2})
dervPlot(untreated_100{3,2}, untreated_100{3, 7}, colorcode{3})
dervPlot(untreated_100{4,2}, untreated_100{4, 7}, colorcode{4})
dervPlot(untreated_100{5,2}, untreated_100{5, 7}, colorcode{5})
%ylim([0 Inf])
ylabel('Relative Change in Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')

%% compare treated to untreated and 2 minute PBS
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
meanPlot(untreated_100{1,2}, untreated_100{1, 7}, colorcode{1}, colorcode2{1}, transparency)
meanPlot(untreated_100{2,2}, untreated_100{2, 7}, colorcode{2}, colorcode2{2}, transparency)
meanPlot(PBS_100{1,2}, PBS_100{1, 7}, colorcode{3}, colorcode2{3}, transparency)
meanPlot(PBS_100{2,2}, PBS_100{2, 7}, colorcode{4}, colorcode2{4}, transparency)
meanPlot(treated_100{1,2}, treated_100{1, 7}, colorcode{5}, colorcode2{5}, transparency)
meanPlot(treated_100{2,2}, treated_100{2, 7}, colorcode{6}, colorcode2{6}, transparency)
meanPlot(treated_100{3,2}, treated_100{3, 7}, colorcode{7}, colorcode2{7}, transparency)
meanPlot(treated_100{4,2}, treated_100{4, 7}, colorcode{8}, colorcode2{8}, transparency)
meanPlot(treated_100{5,2}, treated_100{5, 7}, colorcode{9}, colorcode2{9}, transparency)
ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')
saveas(gcf, 'treated_vs_untreated.png')
saveas(gcf, 'treated_vs_untreated.fig')

%% Functions
%to aggregate data and normalize
function [intensity, adjintensity, normintensity, lCell, time, tme, imstart]=dataNormalize(datadir, imstart, idx)
        
    %pre-allocate variables
    intensity=[];
    lCell=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, 'icell_intensity', 'time', 'lcell')
        
        [nrow, ~]=size(icell_intensity);
        
        for n=1:nrow
            if ~isnan(icell_intensity(n, imstart))
                intensity=[intensity; icell_intensity(n, :)];
                lCell=[lCell; lcell(n, :)];
            end
        end

        if i==1
            tme=time; %pre-set new time vector
        end

    end 
       
    %subtract final fluor. value from the trace
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<0)=NaN;      
        
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:end)./adjintensity(:, imstart);
    [nrow, ncol]=size(normintensity);
    
    %interpolate the fluor values during detergent perfusion 
    for n=1:height(normintensity)        
            if ~isempty(idx)
                x=setdiff(1:ncol, idx); %x=time
                v=normintensity(n, x); %v=intensity
                vq=interp1(x, v, idx); %vq=interpolation at query time pts
                normintensity(n,idx)=vq;
            end
    end
    
    %adjust the time vector
    tme=tme(imstart:end)-tme(imstart);     
        
end

%to aggregate controls data and normalize
function [intensity, bgintensity, adjintensity, normintensity, lCell, time, tme, imstart]=controlsNormalize(datadir)

    %pre-allocate variables
    intensity=[];
    lCell=[];
    bgintensity=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, 'icell_intensity', 'bg_intensity', 'time', 'lcell')
        
        intensity=[intensity; icell_intensity];
        lCell=[lCell; lcell];
        bgintensity=[bgintensity; bg_intensity];

        if i==1
            tme=time; %pre-set new time vector
        end

    end
        

    %find the initial post-lysis frame by identifying the first
    %time point where dt matches the final dt
    dt=round(diff(time));
    if dt(1)==dt(end)
%         dl=diff(lCell, 1);
%         lvg=mean(dl, 1, 'omitnan');
%         imstart=find(lvg<0, 1, 'first')+2;
        imstart=1;
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        intensity=intensity(idx,:);
        lCell=lCell(idx,:);    

        %subtract final fluor. value from the trace and set limit of
        %detection
        adjintensity=intensity-bgintensity; 
        adjintensity(adjintensity<200)=NaN; 
    
    else
        if dt(end)<1 %discrepancies in dt might be found in the meta data, this is the ad hoc fix
            imstart=find(dt==dt(end), 1, 'first')+2;
        else 
            imstart=find(dt==dt(end), 1, 'first');
        end
    
  
    %remove cells without a value for imstart
    idx=find(~isnan(intensity(:, imstart)));
    intensity=intensity(idx,:);
    lCell=lCell(idx,:);    
        
    %subtract final fluor. value from the trace and set limit of
    %detection
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<200)=NaN; 
    
    end
        
    %find when all values are below the limit of detection
    [~, ncol]=size(adjintensity);
    nsum=sum(isnan(adjintensity));
    if max(nsum)==ncol
        imend=min(find(nsum==ncol));
    else 
        imend=ncol;
    end
    
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:imend)./adjintensity(:, imstart);

    %adjust the time vector
    tme=tme(imstart:imend)-tme(imstart);     
        
end

%do not fit to y=(1-beta)*exp(-t./tau)+beta (beta=normalized limit of detection)
%fit to y=Ae^(-t/tau)
function [tau, yhat]=tauCalc(tme, normintensity, tau0)
    
    [nrow, ~]=size(normintensity);
    tau=nan(nrow,1);
    yhat=nan(size(normintensity));
    modelfun = @(tau, x)exp(-x./tau);
    for i=1:nrow
        tau(i, 1)=nlinfit(tme, normintensity(i,:), modelfun, tau0);
         yhat(i, :)=modelfun(tau(i), tme);
    end
    
end

%to correct for photobleaching
function [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept)
        
        %pre-allocate variables
        %assume that the initial 'measured' fluorescence values and corrected
        %fluor. values will be equal. I prefer to pre-allocate with nan in case 
        %some values are missing in the raw data
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:,1);
             
        dCB=nan(height(normintensity), length(tme)-1); %change in fluor. due to photobleaching
        dCT=nan(height(normintensity), length(tme)-1); %total change in fluor.
        dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability

        unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
        unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

        Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
        Cbl_exp(:, 1)=0;
        
        %calculate dt (the dt between frames may vary in a single run).
        %Note: this is also the frame rate
        dt=round(diff(tme), 2);
        
        %this formula comes from the slope and intercept calculated for the 1.2, 2,
        %and 3 second tau vs frame rate controls 
       dC=@(C, alpha, dt, b)(C/(alpha*dt+b))*dt;
       %dC=@(C, alpha, intercept, dt, b)(C/(alpha*dt+b));
        
        %the correction
        for n=1:height(normintensity)
           
            for i=1:length(tme)-1
                
                dCB(n,i) = dC(normintensity(n,i), alpha, dt(i), intercept); %this is the amount of photobleaching that occured in our measured value
   
                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value

                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)*unb_frac(n,i); %Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
                
            end  
            
        end      

end

%to find the time point at which the corrected trace asymptotes and calculate
%the average
function [avg, imend]=traceAvg(Cnew)
    cmean=mean(Cnew, 1, 'omitnan');
    imend=find(cmean-cmean(end)<0.01, 1, 'first');
    if isempty(imend)
        lidx=find(cmean>=0, 1, 'last');
        imend=find(cmean-cmean(lidx)<0.01, 1, 'first');
    end
    avg=cmean(:, 1:imend);
end

%to calculate strain
function [strain]=strainCalc(lcell, imstart)

    lcell=lcell(:, 1:imstart+3);
%     dl=diff(lcell, 1, 2);
% 
%     strain=dl./lcell(:, 1:end-1);

    [nrow, ncol]=size(lcell);
    strain=nan(nrow, ncol);

    for i=1:ncol-1
        dl=lcell(:,i)-lcell(:, i+1);
        strain(:, i)=dl./lcell(:, i+1);
    end

end

%to calculation growth rate
function [growthRate]=gRateCalc(time, lcell, imstart)

    lcell=lcell(:, 1:imstart+3);
    time=time(1:imstart+3);

    dl=diff(lcell, 1, 2);
    dt=diff(time, 1, 2);

    sl=(lcell(:, 1:end-1)+(dl+lcell(:, 1:end-1)))/2;

    growthRate=(dl./dt)./sl;

end

%to predict a normintensity trace
function [Pnew]=tracePredict(time, tau, alpha, intercept)
    
    diffusion=@(time, tau)exp(-time/tau);
    trueTrace=diffusion(time, tau);
    Pnew=trueTrace;
    
    dt=diff(time, 1, 2);
    dC=@(C, alpha, dt, b)(C/(alpha*dt+b))*dt;
    unb_frac=1;
    
    for i=1:length(time)-1
        dCP=trueTrace(1, i+1)-trueTrace(1,i);
        dCP=dCP*unb_frac;
     
        dCB = dC(trueTrace(1,i), alpha, dt(i), intercept);
        Pnew(1, i+1) = Pnew(1, i) - dCB;
        
        Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
        unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(n,i+1));%Calculate the new fraction of unbleached fluorophores
    
    end

end

%to plot mean and std 
function meanPlot(tme, normintensity, colorcode1, colorcode2, transparency)

nmean=mean(normintensity, 1, 'omitnan');
nstd=std(normintensity, 0, 1, 'omitnan');

if nargin==5
    ciplot(nmean-nstd, nmean+nstd, tme, colorcode2, transparency)
    plot(tme, nmean, 'Color', colorcode1, 'LineWidth', 1.5)
elseif nargin==4
    errorbar(tme, nmean, -nstd, nstd, 'Color', colorcode1, 'LineWidth', 1.5)
end

end

%to plot the derivative of traces
function dervPlot(tme, normintensity, colorcode)

    %nmean=mean(normintensity, 1, 'omitnan');
    %nstd=std(normintensity, 0, 1, 'omitnan');
    %dmean=diff(nmean, 1);
    
    dtrace=diff(normintensity, 1, 2);
    dt=diff(tme);
    %rdtrace=(dtrace./dt)./normintensity(:, 1:end-1); %relative change
    rdtrace=dtrace./normintensity(:, 1:end-1);
    rdmean=mean(rdtrace./dt, 1, 'omitnan'); 

    plot(tme(1:end-1), rdmean, 'Color', colorcode, 'LineWidth', 1.5)

end

