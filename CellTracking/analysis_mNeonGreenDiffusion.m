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

%LB + 20 mM Mg2+, 1 minute frame rate, 100% intensity = '01172022_Exp1' 
%LB + 10 mM EDTA, 1 minute frame rate, 100% intensity = '01172022_Exp2' 
%LB + 0.5 ug/mL tunicamycin, 1 minute frame rate, 100% intensity = '01242022_Exp1' 
%LB + 1 ug/mL vancomycin, 1 minute frame rate, 100% intensity = '01262022_Exp1'
%spent LB, 10 minute frame rate, 100% intensity = '02122022_Exp3' & '02192022_Exp1'
%spent LB, 1 minute frame rate, 100% intensity = '03012022_Exp1'

%untreated, 1.2 s frame rate, 100% intensity = '11202021_Exp1' 
%untreated, 2 s frame rate, 100% intensity = '12082021_Exp1'
%untreated, 3 s frame rate, 100% intensity = '12082021_Exp2'

%untreated, 1.76 s frame rate, 20% intensity = '02192022_Exp2' 
%untreated, 2.3 s frame rate, 20% intensity = '02192022_Exp3'
%untreated, 3 s frame rate, 20% intensity = '02192022_Exp4'

%for photobleach correction (100% intensity)
% alpha=32.2114;
% intercept=0.1614;

%for photobleach correction (20% intensity)
% alpha=151.7544;
% intercept=-1.0566;

%% User Input 
%color codes must be in RGB [0-1] format to be used in ciplot
color_gray = {[204 204 204], [150, 150, 150], [82, 82, 82]};
color_gray = cellfun(@(x)x./255, color_gray, 'UniformOutput', false);

colorcode={[204 0 0], [204 102 0], [204 204 0], [102 204 0], [0 204 204], [0 0 204], [102 0 204], [204 0 204], [204 0 102], [255 102 102], [255 178 102], [102 255 102], [102 255 255], [102 178 255], [178 102 255], [255 102 255],[255 102 178]};
colorcode2={[255 51 51], [255 153 51], [255 255 51], [153 255 51], [51 255 255], [51 51 255], [153 51 255], [255 51 255], [255 51 153], [255 204 204], [255 229 204], [204 255 204], [204 255 255], [204 229 255], [229 204 255], [255 204 255],[255 204 229]};

colorcode=cellfun(@(x)(x./255), colorcode, 'UniformOutput', false);
colorcode2=cellfun(@(x)(x./255), colorcode2, 'UniformOutput', false);

okabeIto = {[230, 159, 0], [86, 180, 233], [0, 158, 115], [240, 228, 66], [0, 114, 178], [213, 94, 0], [204, 121, 167], [0, 0, 0]};
okabeIto = cellfun(@(x)(x./255), okabeIto, 'UniformOutput', false);

% m=10;
% c1=[0.847 0.427 0.933];
% c2=[0.101 0.101 0.737];
% color_p=[linspace(c1(1), c2(1), m)', linspace(c1(2), c2(2), m)', linspace(c1(3), c2(3), m)'];

color_red={[246, 177, 195], [222, 38, 76], [162, 13, 30]};
color_red=cellfun(@(x)x./255, color_red, 'UniformOutput', false);

color_blue={[173, 213, 247], [78, 122, 199], [0, 48, 86]};
color_blue=cellfun(@(x)x./255, color_blue, 'UniformOutput', false);

color_purple={[254 196 254], [255 102 209],[228 91 255], [179 49 255], [109 0 175]};
color_purple=cellfun(@(x)x./255, color_purple, 'UniformOutput', false);

color_green={[0.3804 1.0000 0.9020], [0.4196 0.6824 0.8392], [0.8392 0.8549 0], [0.4549 0.7686 0.4627], [0.1569  0.7804  0.6667], [0 0.4275 0.1725], [0.5725 0.4706 1.0000]};
%color_green=cellfun(@(x)x./255, color_green, 'UniformOutput', false);

transparency = 0.3; %this is the alpha argument for the ciplot function

%location and names of of the *.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis';
basenames={'10232021_Exp1', '10262021_Exp1', '02212022_Exp2', '02122022_Exp1','02122022_Exp2', '02092022_Exp1', '02212022_Exp1', '11192021_Exp2', '12082021_Exp3', '11192021_Exp1', '11302021_Exp1', '10232021_Exp2', '10262021_Exp2', '01142022_Exp1', '01172022_Exp1', '01172022_Exp2', '01242022_Exp1', '01262022_Exp1', '02122022_Exp3', '02192022_Exp1', '03012022_Exp1', '11202021_Exp1', '12082021_Exp1', '12082021_Exp2', '02192022_Exp2', '02192022_Exp3', '02192022_Exp4'};

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
    
    cidx1 = [length(basenames)-5:length(basenames)-3];
    cidx2 = [length(basenames)-2:length(basenames)];
    
    if ismember(i, cidx1)
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
        
    elseif ismember(i, cidx2)
        
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
%different variable (time, tme, intensity, adjintensity, normintensity, lcell, Cnew, and strain)

untreated_100 = cell(5, 8);
untreated_20 = cell(2, 8);

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
        untreated_100{idx1, 9}=bgintensity;
        untreated_100{idx1, 10}=gRateCalc(time, lcell, imstart);
        
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
        untreated_20{idx2, 9}=bgintensity;
        
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
%different variable (time, tme, intensity, adjintensity, normintensity, lcell, Cnew, and strain)

PBS_100 = cell(7, 8);

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
        PBS_100{idx1, 8}=gRateCalc(time, lcell, imstart);
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        PBS_100{idx1, 7}=Cnew;
        
        idx1=idx1+1;
        
     end
    
end

%% compare the LB treated traces
%each row is a different treatment in LB (100% intensity)

treated_100 = cell(6, 8);

idx1=1;

for i=1:length(basenames)
    basename=basenames{i};
    
     if ismember(i, [15:21])
        cd([dirsave '/normalizedFiles'])
        load([basename '_norm.mat'])
        treated_100{idx1, 1}=time;
        treated_100{idx1, 2}=tme;
        treated_100{idx1, 3}=intensity;
        treated_100{idx1, 4}=adjintensity;
        treated_100{idx1, 5}=normintensity;
        treated_100{idx1, 6}=lcell;
        treated_100{idx1, 8}=gRateCalc(time, lcell, imstart);
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        treated_100{idx1, 7}=Cnew;
        
        idx1=idx1+1;
        
     end
    
end

%% generate plots for the controls
cd([dirsave '/03262022_analysis']);

labels={'1 s', '2 s', '3 s', '1 s background', '2 s background', '3 s background'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',20)
subplot(1, 2, 1), hold on
p1=meanPlot(controls_100{1,2}, controls_100{1,3}(:, controls_100{1,9}:end), color_red{1}, color_red{1}) %plot starting from the initial post-lysis value
p2=meanPlot(controls_100{2,2}, controls_100{2,3}(:, controls_100{2,9}:end), color_red{2}, color_red{2})
p3=meanPlot(controls_100{3,2}, controls_100{3,3}(:, controls_100{3,9}:end), color_red{3}, color_red{3})
p4=plot(controls_100{1,2}, controls_100{1,4}(:, controls_100{1,9}:end), '--', 'Color', color_red{1}, 'LineWidth', 1.5)
p5=plot(controls_100{2,2}, controls_100{2,4}(:, controls_100{2,9}:end), '--', 'Color', color_red{2}, 'LineWidth', 1.5)
p6=plot(controls_100{3,2}, controls_100{3,4}(:, controls_100{3,9}:end), '--', 'Color', color_red{3}, 'LineWidth', 1.5)
ylim([0 Inf])
ylabel('Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('100% Power')
hleg=legend([p1(1), p2(1), p3(1), p4, p5, p6], labels)
%hleg.NumColumns=2;
title(hleg,'Settings')

% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
subplot(1, 2, 2), hold on
p1=meanPlot(controls_20{1,2}, controls_20{1,3}(:, controls_20{1,9}:end), color_blue{1}, color_blue{1})
p2=meanPlot(controls_20{2,2}, controls_20{2,3}(:, controls_20{2,9}:end), color_blue{2}, color_blue{2})
p3=meanPlot(controls_20{3,2}, controls_20{3,3}(:, controls_20{3,9}:end), color_blue{3}, color_blue{3})
p4=plot(controls_20{1,2}, controls_20{1,4}(:, controls_20{1,9}:end), '--', 'Color', color_blue{1}, 'LineWidth', 1.5)
p5=plot(controls_20{2,2}, controls_20{2,4}(:, controls_20{2,9}:end), '--', 'Color', color_blue{2}, 'LineWidth', 1.5)
p6=plot(controls_20{3,2}, controls_20{3,4}(:, controls_20{3,9}:end), '--', 'Color', color_blue{3}, 'LineWidth', 1.5)
ylim([0 Inf])
ylabel('Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('20% Power')
hleg=legend([p1(1), p2(1), p3(1), p4, p5, p6], labels)
title(hleg,'Settings')
%hleg.NumColumns=2;

saveas(gcf, 'controlsRaw.png')
saveas(gcf, 'controlsRaw.fig')

labels={'100% 1 s', '100% 2 s', '100% 3 s', '20% 1 s', '20% 2 s', '20% 3 s'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
p1=meanPlot(controls_100{1,2}, controls_100{1, 6}, color_red{1}, color_red{1})
p2=meanPlot(controls_100{2,2}, controls_100{2, 6}, color_red{2}, color_red{2})
p3=meanPlot(controls_100{3,2}, controls_100{3, 6}, color_red{3}, color_red{3})
p4=meanPlot(controls_20{1,2}, controls_20{1, 6}, color_blue{1}, color_blue{1})
p5=meanPlot(controls_20{2,2}, controls_20{2, 6}, color_blue{2}, color_blue{2})
p6=meanPlot(controls_20{3,2}, controls_20{3, 6}, color_blue{3}, color_blue{3})
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
hleg=legend([p1, p2, p3, p4, p5, p6], labels)
%hleg.NumColumns=2;
title(hleg,'Settings')
saveas(gcf, 'controlsNorm.png')
saveas(gcf, 'controlsNorm.fig')


pt=8;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
idx=find(controls_100{1,2}<pt);
p1=meanPlot(controls_100{1,2}(:, idx), controls_100{1, 8}(:, idx), color_red{1}, color_red{1}, transparency)
idx=find(controls_100{2,2}<pt);
p2=meanPlot(controls_100{2,2}(:, idx), controls_100{2, 8}(:, idx), color_red{2}, color_red{2}, transparency)
idx=find(controls_100{3,2}<pt);
p3=meanPlot(controls_100{3,2}(:, idx), controls_100{3, 8}(:, idx), color_red{3}, color_red{3}, transparency)
idx=find(controls_20{1,2}<pt);
p4=meanPlot(controls_20{1,2}(:, idx), controls_20{1, 8}(:, idx), color_blue{1}, color_blue{1}, transparency)
idx=find(controls_20{2,2}<pt);
p5=meanPlot(controls_20{2,2}(:, idx), controls_20{2, 8}(:, idx), color_blue{2}, color_blue{2}, transparency)
idx=find(controls_20{3,2}<pt);
p6=meanPlot(controls_20{3,2}(:, idx), controls_20{3, 8}(:, idx), color_blue{3}, color_blue{3}, transparency)
%ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')
hleg=legend([p1, p2, p3, p4, p5, p6], labels, 'Location', 'east')
title(hleg,'Settings')
%hleg.NumColumns=2;
saveas(gcf, 'controlsCorrected_zoom.png')
saveas(gcf, 'controlsCorrected_zoom.fig')

pt=8;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
idx=find(controls_100{1,2}<pt);
p1=meanPlot(controls_100{1,2}(:, idx), controls_100{1, 8}(:, idx), color_red{1}, color_red{1}, transparency)
idx=find(controls_100{2,2}<pt);
p2=meanPlot(controls_100{2,2}(:, idx), controls_100{2, 8}(:, idx), color_red{2}, color_red{2}, transparency)
idx=find(controls_100{3,2}<pt);
p3=meanPlot(controls_100{3,2}(:, idx), controls_100{3, 8}(:, idx), color_red{3}, color_red{3}, transparency)
idx=find(controls_20{1,2}<pt);
p4=meanPlot(controls_20{1,2}(:, idx), controls_20{1, 8}(:, idx), color_blue{1}, color_blue{1}, transparency)
idx=find(controls_20{2,2}<pt);
p5=meanPlot(controls_20{2,2}(:, idx), controls_20{2, 8}(:, idx), color_blue{2}, color_blue{2}, transparency)
idx=find(controls_20{3,2}<pt);
p6=meanPlot(controls_20{3,2}(:, idx), controls_20{3, 8}(:, idx), color_blue{3}, color_blue{3}, transparency)
ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')
hleg=legend([p1, p2, p3, p4, p5, p6], labels, 'Location', 'east')
title(hleg,'Settings')
%hleg.NumColumns=2;
saveas(gcf, 'controlsCorrected_pan.png')
saveas(gcf, 'controlsCorrected_pan.fig')

%% generate plots for untreated LB
cd([dirsave '/03262022_analysis']);

labels={'1 min', '1 min', '5 min', '10 min', '20 min'};

% for i=1:height(untreated_100)
%     figure, hold on
%     plot(untreated_100{i,1}, untreated_100{i, 3}, '-k')
%     plot(untreated_100{i,1}, untreated_100{i, 9}, '--k')
%     pause, close all
% end

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
p1=meanPlot(untreated_100{1,1}, untreated_100{1, 3}, color_purple{1}, color_purple{1}, transparency)
p2=meanPlot(untreated_100{2,1}, untreated_100{2, 3}, color_purple{3}, color_purple{3}, transparency)
p3=meanPlot(untreated_100{3,1}, untreated_100{3, 3}, color_gray{1}, color_gray{1}, transparency)
p4=meanPlot(untreated_100{4,1}, untreated_100{4, 3}, color_gray{2}, color_gray{2}, transparency)
p5=meanPlot(untreated_100{5,1}, untreated_100{5, 3}, color_gray{3}, color_gray{3}, transparency)
ylim([0 Inf])
ylabel('Fluorescence (A.U.)')
xlabel('Time (minutes)')
hleg=legend([p1(1), p2(1), p3(1), p4(1), p5(1)], labels)
title(hleg,'Frame Rate')
saveas(gcf, 'untreatedRaw.png')
saveas(gcf, 'untreatedRaw.fig')
% 
% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
% %meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
% meanPlot(untreated_100{1,1}, untreated_100{1, 4}, colorcode{1}, colorcode2{1}, transparency)
% meanPlot(untreated_100{2,1}, untreated_100{2, 4}, colorcode{2}, colorcode2{2}, transparency)
% meanPlot(untreated_100{3,1}, untreated_100{3, 4}, colorcode{3}, colorcode2{3}, transparency)
% meanPlot(untreated_100{4,1}, untreated_100{4, 4}, colorcode{4}, colorcode2{4}, transparency)
% meanPlot(untreated_100{5,1}, untreated_100{5, 4}, colorcode{5}, colorcode2{5}, transparency)
% %meanPlot(untreated_100{5,1}, (untreated_100{5, 3}-mean(untreated_100{4, 3}(:,end), 1, 'omitnan')), colorcode{5}, colorcode2{5}, transparency)
% ylim([0 Inf])
% ylabel('Adjusted Fluorescence (A.U.)')
% xlabel('Time (minutes)')
% hleg=legend(labels)
% title(hleg,'Frame Rate')
% saveas(gcf, 'untreatedAdj.png')
% saveas(gcf, 'untreatedAdj.fig')
% 
% % 
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24)
subplot(1,2,1)
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 5}], colorcode{1}, colorcode2{1}, transparency)
p1=meanPlot(untreated_100{1,2}, untreated_100{1, 5}, color_purple{1}, color_purple{1}, transparency), hold on
p2=meanPlot(untreated_100{2,2}, untreated_100{2, 5}, color_purple{3}, color_purple{3}, transparency), hold on
p3=meanPlot(untreated_100{3,2}, untreated_100{3, 5}, color_gray{1}, color_gray{1}, transparency), hold on
p4=meanPlot(untreated_100{4,2}, untreated_100{4, 5}, color_gray{2}, color_gray{2}, transparency), hold on
p5=meanPlot(untreated_100{5,2}, untreated_100{5, 5}, color_gray{3}, color_gray{3}, transparency), hold on
ylim([0 Inf])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
hleg=legend([p1(1), p2(1), p3(1), p4(1), p5(1)], labels)
title(hleg,'Frame Rate')
% saveas(gcf, 'normUntreated.png')
% saveas(gcf, 'normUntreated.fig')
% 
% figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
% %meanPlot(untreated_100{2,2}, [untreated_100{1, 7}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
subplot(1,2,2)
p1=meanPlot(untreated_100{1,2}, untreated_100{1, 7}, color_purple{1}, color_purple{1}, transparency), hold on
p2=meanPlot(untreated_100{2,2}, untreated_100{2, 7}, color_purple{3}, color_purple{3}, transparency), hold on
p3=meanPlot(untreated_100{3,2}, untreated_100{3, 7}, color_gray{1}, color_gray{1}, transparency), hold on
p4=meanPlot(untreated_100{4,2}, untreated_100{4, 7}, color_gray{2}, color_gray{2}, transparency), hold on
p5=meanPlot(untreated_100{5,2}, untreated_100{5, 7}, color_gray{3}, color_gray{3}, transparency)
ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')
hleg=legend([p1(1), p2(1), p3(1), p4(1), p5(1)], labels)
title(hleg,'Frame Rate')
% % saveas(gcf, 'correctedUntreated.png')
% % saveas(gcf, 'correctedUntreated.fig')
% 
saveas(gcf, 'untreatedNorm&Corrected.png')
saveas(gcf, 'untreatedNorm&Corrected.fig')

% figure, hold on
% meanPlot(untreated_100{2,2}, [untreated_100{1, 7}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
% meanPlot(controls_100{1,2}, controls_100{1, 8}, colorcode{1}, colorcode2{1}, transparency)
% meanPlot(controls_100{2,2}, controls_100{2, 8}, colorcode{3}, colorcode2{3}, transparency)
% meanPlot(controls_100{3,2}, controls_100{3, 8}, colorcode{5}, colorcode2{5}, transparency)

%% plot the fit of the exponential on the control traces
cd([dirsave '/03262022_analysis']);

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',16)
subplot(2,3,1)
plot(controls_100{1,2}, controls_100{1, 6}, 'Color', color_red{1}), hold on
plot(controls_100{1,2}, controls_100{1, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('100% Intensity, 1 s Frame Rate')
% saveas(gcf, 'fit100_1s.png')
% saveas(gcf, 'fit100_1s.fig')

subplot(2,3,2)
plot(controls_100{2,2}, controls_100{2, 6}, 'Color', color_red{2}), hold on
plot(controls_100{2,2}, controls_100{2, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('100% Intensity, 2 s Frame Rate')
% saveas(gcf, 'fit100_2s.png')
% saveas(gcf, 'fit100_2s.fig')

subplot(2,3,3)
plot(controls_100{3,2}, controls_100{3, 6}, 'Color', color_red{3}), hold on
plot(controls_100{3,2}, controls_100{3, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('100% Intensity, 3 s Frame Rate')
% saveas(gcf, 'fit100_3s.png')
% saveas(gcf, 'fit100_3s.fig')

subplot(2,3,4)
plot(controls_20{1,2}, controls_20{1, 6}, 'Color', color_blue{1}), hold on
plot(controls_20{1,2}, controls_20{1, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('20% Intensity, 1 s Frame Rate')
% saveas(gcf, 'fit20_1s.png')
% saveas(gcf, 'fit20_1s.fig')

subplot(2,3,5)
plot(controls_20{2,2}, controls_20{2, 6}, 'Color', color_blue{2}), hold on
plot(controls_20{2,2}, controls_20{2, 11}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')
title('20% Intensity, 2 s Frame Rate')
% saveas(gcf, 'fit20_2s.png')
% saveas(gcf, 'fit20_2s.fig')

subplot(2,3,6)
plot(controls_20{3,2}, controls_20{3, 6}, 'Color', color_blue{3}), hold on
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
cd([dirsave '/03262022_analysis']);

tau_means = cellfun(@(x)mean(x, 1, 'omitnan'), tau);
tau_std = cellfun(@(x)std(x, 0, 1, 'omitnan'), tau);

labels={'\tau 100% power', '\tau 20% power', 'mean \tau', 'line of best fit'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
p1=scatter(repelem(dt(1), length(tau{1})), tau{1}, 'MarkerFaceColor', '#E74C3C', 'MarkerEdgeColor', '#E74C3C')
scatter(repelem(dt(2), length(tau{2})), tau{2}, 'MarkerFaceColor', '#E74C3C', 'MarkerEdgeColor', '#E74C3C')
scatter(repelem(dt(3), length(tau{3})), tau{3}, 'MarkerFaceColor', '#E74C3C', 'MarkerEdgeColor', '#E74C3C')

p2=scatter(repelem(dt(4), length(tau{4})), tau{4}, 'MarkerFaceColor', '#3498DB', 'MarkerEdgeColor', '#3498DB')
scatter(repelem(dt(5), length(tau{5})), tau{5}, 'MarkerFaceColor', '#3498DB', 'MarkerEdgeColor', '#3498DB')
scatter(repelem(dt(6), length(tau{6})), tau{6}, 'MarkerFaceColor', '#3498DB', 'MarkerEdgeColor', '#3498DB')

p3=scatter([dt(1), dt(2), dt(3), dt(4), dt(5), dt(6)], tau_means, 'MarkerFaceColor', 'black')
errorbar(dt(1), tau_means(1), tau_std(1), 'Color', 'black')
errorbar(dt(2), tau_means(2), tau_std(2), 'Color', 'black')
errorbar(dt(3), tau_means(3), tau_std(3), 'Color', 'black')
errorbar(dt(4), tau_means(4), tau_std(4), 'Color', 'black')
errorbar(dt(5), tau_means(5), tau_std(5), 'Color', 'black')
errorbar(dt(6), tau_means(6), tau_std(6), 'Color', 'black')

p4=plot([0, dt(1), dt(2), dt(3)], linearFit1, '--k', 'LineWidth', 1)
plot([0, dt(4), dt(5), dt(6)], linearFit2, '--k', 'LineWidth', 1)
xlim([0, 0.06])
ylabel('\tau (min^{-1})')
hleg=legend([p1(1), p2(1), p3, p4], labels, 'Location', 'southeast')
%xticks([0, dt(1),  dt(4), dt(2),  dt(5), dt(3)])
%xticklabels({'0', '1.2 s', '1.76 s', '2 s', '2.3 s', '3 s'})
xlabel('Frame Rate (minutes)')

caption1 = sprintf('%f * frame rate + %f', linearCoef1(1), linearCoef1(2));
caption2 = sprintf('%f * frame rate + %f', linearCoef2(1), linearCoef2(2));
text(dt(2), tau_means(2)+1.5, ['\tau = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
text(dt(1), tau_means(5)+1.5, ['\tau = ' caption2], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

title('Tau vs Frame Rate')
saveas(gcf, 'tau_vs_frameRate.png')
saveas(gcf, 'tau_vs_frameRate.fig')

%% calculate the autofluorescence value as a function of frame rate
%untreated 100% intensity with frame rates at: 1 s, 2 s, 3 s, 1 min, 5 min,
%10 min and 20 min (although exclude the 20 minute if need be since it
%doesn't quite autofluorescence)

autofluorescence=cell(7, 2); %column 1 = dt, column 2 = final raw fluor values

dt2=[1.2/60, 2/60, 3/60, 1, 5, 10, 20];

autofluorescence{1,2}=controls_100{1,4}(:, end);
autofluorescence{2,2}=controls_100{2,4}(:, end);
autofluorescence{3,2}=controls_100{3,4}(:, end);
autofluorescence{4,2}=untreated_100{1,9}(:, end);
autofluorescence{5,2}=untreated_100{3,9}(:, end);
autofluorescence{6,2}=untreated_100{4,9}(:, end);
autofluorescence{7,2}=untreated_100{5,9}(:, end);

for i=1:length(dt2)
    autofluorescence{i,1}=repelem(dt2(i), height(autofluorescence{i,2}))';
end

% autofluorescence_x=dt2([4, 5, 6, 7]);
% autofluorescence_y=cellfun(@(x)mean(x, 1, 'omitnan'), autofluorescence([4, 5, 6, 7], 2))';
% % linearCoef3 = polyfit(autofluorescence_x, autofluorescence_y, 1);
% % linearFit3 = polyval(linearCoef3,[0 dt2]);
% 
% modelfun=@(b,x)(b(1)*x)./(b(2)+x); %why does order matter here?
% beta0=[1300, 0.01];
% beta = nlinfit(autofluorescence_x, autofluorescence_y, modelfun, beta0);
% autofluorescence_yhat=modelfun(beta, [0 autofluorescence_x]);

%% generate plots to illustrate the autofluorescence value as a function of frame rate
cd([dirsave '/03262022_analysis']);

labels={'autofluorescence value for single cell', 'mean autofluorescence value'};

autofluorescence_means = cellfun(@(x)mean(x, 1, 'omitnan'), autofluorescence(:, 2));
autofluorescence_std = cellfun(@(x)std(x, 0, 1, 'omitnan'), autofluorescence(:, 2));

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
p1=scatter(autofluorescence{1, 1}, autofluorescence{1, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(autofluorescence{2, 1}, autofluorescence{2, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(autofluorescence{3, 1}, autofluorescence{3, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(autofluorescence{4, 1}, autofluorescence{4, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(autofluorescence{5, 1}, autofluorescence{5, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(autofluorescence{6, 1}, autofluorescence{6, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(autofluorescence{7, 1}, autofluorescence{7, 2}, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')

p2=scatter([dt2(1), dt2(2), dt2(3), dt2(4), dt2(5), dt2(6), dt2(7)], autofluorescence_means, 'MarkerFaceColor', 'black')
errorbar(dt2(1), autofluorescence_means(1), autofluorescence_std(1), 'Color', 'black')
errorbar(dt2(2), autofluorescence_means(2), autofluorescence_std(2), 'Color', 'black')
errorbar(dt2(3), autofluorescence_means(3), autofluorescence_std(3), 'Color', 'black')
errorbar(dt2(4), autofluorescence_means(4), autofluorescence_std(4), 'Color', 'black')
errorbar(dt2(5), autofluorescence_means(5), autofluorescence_std(5), 'Color', 'black')
errorbar(dt2(6), autofluorescence_means(6), autofluorescence_std(6), 'Color', 'black')
errorbar(dt2(7), autofluorescence_means(7), autofluorescence_std(7), 'Color', 'black')

%plot([0 autofluorescence_x], autofluorescence_yhat, '--k', 'LineWidth', 1)
ylabel('Autofluorescence Value (A.U.)')
xlabel('Frame Rate (minutes)')
ylim([0 1800])
xlim([-1 dt2(end)+1])

hleg=legend([p1(1), p2], labels, 'Location', 'east')

% caption3 = sprintf('%f * frame rate + %f', linearCoef3(1), linearCoef3(2));
% text(dt2(2), autofluorescence_means(2)+1.5, ['\autofluorescence = ' caption3], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

title('Autofluorescence Value vs Frame Rate')
% saveas(gcf, 'autofluorescence_vs_frameRate.png')
% saveas(gcf, 'autofluorescence_vs_frameRate.fig')

%% generate plots for PBS
cd([dirsave '/03262022_analysis']);

%labels={'untreated', 'untreated', '2 min PBS', '2 min PBS', '20 min PBS', '20 min PBS', '60 min PBS', '60 min PBS', '120 min PBS'}; 
labels={'untreated', '2 min PBS', '20 min PBS', '60 min PBS', '120 min PBS'}; 
imend=48;

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
p1=meanPlot(untreated_100{1,1}(1:imend), [untreated_100{1, 3}(:, 1:imend);untreated_100{2, 3}(:, 1:imend)], colorcode{1}, colorcode2{1}, transparency)
%p2=meanPlot(untreated_100{2,1}, untreated_100{2, 3}, color_purple{2})
p3=meanPlot(PBS_100{1,1}(1:imend), [PBS_100{1, 3}(:, 1:imend); PBS_100{2, 3}(:, 1:imend)], colorcode{2}, colorcode2{2}, transparency)
%p4=meanPlot(PBS_100{2,1}, PBS_100{2, 3}, color_green{2})
p5=meanPlot(PBS_100{3,1}(1:imend), [PBS_100{3, 3}(:, 1:imend); PBS_100{4, 3}(:, 1:imend)], colorcode{4}, colorcode2{4}, transparency)
%p6=meanPlot(PBS_100{4,1}, PBS_100{4, 3}, color_green{4})
p7=meanPlot(PBS_100{5,1}(1:imend), [PBS_100{5, 3}(:, 1:imend); PBS_100{6, 3}(:, 1:imend)], colorcode{5}, colorcode2{5}, transparency)
%p8=meanPlot(PBS_100{6,1}, PBS_100{6, 3}, color_green{6})
p9=meanPlot(PBS_100{7,1}, PBS_100{7, 3}, colorcode{6}, colorcode2{6}, transparency)
ylim([0 Inf])
ylabel('Fluorescence (A.U.)')
xlabel('Time (minutes)')

%hleg=legend([p1(1), p2(1), p3(1), p4(1), p5(1), p6(1), p7(1), p8(1), p9(1)], labels)
hleg=legend([p1(1), p3(1), p5(1),p7(1), p9(1)], labels)
title(hleg,'PBS Incubation')

saveas(gcf, 'rawPBS.png')
saveas(gcf, 'rawPBS.fig')

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1})
subplot(1,2,1)
p1=meanPlot(untreated_100{1,2}(1:imend), [untreated_100{1, 5}(:, 1:imend); untreated_100{2, 5}(:, 1:imend)], colorcode{1}, colorcode2{1}, transparency), hold on
%p2=meanPlot(untreated_100{2,2}, untreated_100{2, 5}, color_purple{2})
p3=meanPlot(PBS_100{1,2}(1:imend), [PBS_100{1, 5}(:, 1:imend); PBS_100{2, 5}(:, 1:imend)], colorcode{2}, colorcode2{2}, transparency)
%p4=meanPlot(PBS_100{2,2}, PBS_100{2, 5}, color_green{2})
p5=meanPlot(PBS_100{3,2}(1:imend), [PBS_100{3, 5}(:, 1:imend); PBS_100{4, 5}(:, 1:imend)], colorcode{4}, colorcode2{4}, transparency)
%p6=meanPlot(PBS_100{4,2}, PBS_100{4, 5}, color_green{4})
p7=meanPlot(PBS_100{5,2}(1:imend), [PBS_100{5, 5}(:, 1:imend); PBS_100{6, 5}(:, 1:imend)], colorcode{5}, colorcode2{5}, transparency)
%p8=meanPlot(PBS_100{6,2}, PBS_100{6, 5}, color_green{6})
p9=meanPlot(PBS_100{7,2}, PBS_100{7, 5}, colorcode{6}, colorcode2{6}, transparency)
ylim([0 Inf])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')

%hleg=legend([p1(1), p2(1), p3(1), p4(1), p5(1), p6(1), p7(1), p8(1), p9(1)], labels)
hleg=legend([p1(1), p3(1), p5(1),p7(1), p9(1)], labels)
title(hleg,'PBS Incubation')
% 
% saveas(gcf, 'PBSnorm.png')
% saveas(gcf, 'PBSnorm.fig')

%figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
%meanPlot(untreated_100{2,2}, [untreated_100{1, 5}(:, 1:94); untreated_100{2, 7}], colorcode{1}, colorcode2{1}, transparency)
subplot(1,2,2)
p1=meanPlot(untreated_100{1,2}(1:imend), [untreated_100{1, 7}(:, 1:imend); untreated_100{2, 7}(:, 1:imend)], colorcode{1}, colorcode2{1}, transparency), hold on
%p2=meanPlot(untreated_100{2,2}(1:50), untreated_100{2, 7}(:, 1:50), color_purple{2})
p3=meanPlot(PBS_100{1,2}(1:imend), [PBS_100{1, 7}(:, 1:imend); PBS_100{2, 7}(:, 1:imend)], colorcode{2}, colorcode2{2}, transparency)
%p4=meanPlot(PBS_100{2,2}(1:50), PBS_100{2, 7}(:, 1:50), color_green{2})
p5=meanPlot(PBS_100{3,2}(1:imend), [PBS_100{3, 7}(:, 1:imend); PBS_100{4, 7}(:, 1:imend)], colorcode{4}, colorcode2{4}, transparency)
%p6=meanPlot(PBS_100{4,2}(1:50), PBS_100{4, 7}(:, 1:50), color_green{4})
p7=meanPlot(PBS_100{5,2}(1:imend), [PBS_100{5, 7}(:, 1:imend); PBS_100{6, 7}(:, 1:imend)], colorcode{5}, colorcode2{5}, transparency)
%p8=meanPlot(PBS_100{6,2}, PBS_100{6, 7}, color_green{6})
p9=meanPlot(PBS_100{7,2}, PBS_100{7, 7}, colorcode{6}, colorcode2{6}, transparency)
ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')

% hleg=legend([p1(1), p2(1), p3(1), p4(1), p5(1), p6(1), p7(1), p8(1), p9(1)], labels)
% title(hleg,'PBS Incubation')

% saveas(gcf, 'PBScorrected.png')
% saveas(gcf, 'PBScorrected.fig')
% 
saveas(gcf, 'PBSnorm&corrected.png')
saveas(gcf, 'PBSnorm&corrected.fig')

%% when do the untreated cells hit the autofluorescence?
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
cd([dirsave '/03262022_analysis']);

labels={'untreated', 'untreated', '2 min PBS', '2 min PBS', 'Mg^{2+}', 'EDTA', 'tunicamycin', 'vancomycin'}; 

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24)
subplot(1,2,1)
p1=meanPlot(untreated_100{1,2}, untreated_100{1, 5}, colorcode{1}, colorcode2{1}, transparency), hold on
p2=meanPlot(untreated_100{2,2}, untreated_100{2, 5}, colorcode{2}, colorcode2{2}, transparency)
p3=meanPlot(PBS_100{1,2}, PBS_100{1, 5}, colorcode{3}, colorcode2{3}, transparency)
p4=meanPlot(PBS_100{2,2}, PBS_100{2, 5}, colorcode{4}, colorcode2{4}, transparency)
p5=meanPlot(treated_100{1,2}, treated_100{1, 5}, colorcode{5}, colorcode2{5}, transparency)
p6=meanPlot(treated_100{2,2}, treated_100{2, 5}, colorcode{6}, colorcode2{6}, transparency)
p7=meanPlot(treated_100{3,2}, treated_100{3, 5}, colorcode{7}, colorcode2{7}, transparency)
p8=meanPlot(treated_100{4,2}, treated_100{4, 5}, colorcode{8}, colorcode2{8}, transparency)
ylim([0 Inf])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')

hleg=legend([p1(1), p2(1), p3(1), p4(1), p5(1), p6(1), p7(1)], labels)
title(hleg,'Treatment')

%figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
subplot(1,2,2)
p1=meanPlot(untreated_100{1,2}, untreated_100{1, 7}, colorcode{1}, colorcode2{1}, transparency), hold on
p2=meanPlot(untreated_100{2,2}, untreated_100{2, 7}, colorcode{2}, colorcode2{2}, transparency)
p3=meanPlot(PBS_100{1,2}, PBS_100{1, 7}, colorcode{3}, colorcode2{3}, transparency)
p4=meanPlot(PBS_100{2,2}, PBS_100{2, 7}, colorcode{4}, colorcode2{4}, transparency)
p5=meanPlot(treated_100{1,2}, treated_100{1, 7}, colorcode{5}, colorcode2{5}, transparency)
p6=meanPlot(treated_100{2,2}, treated_100{2, 7}, colorcode{6}, colorcode2{6}, transparency)
p7=meanPlot(treated_100{3,2}, treated_100{3, 7}, colorcode{7}, colorcode2{7}, transparency)
p8=meanPlot(treated_100{3,2}, treated_100{3, 7}, colorcode{7}, colorcode2{7}, transparency)
ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')


saveas(gcf, 'treated_vs_untreated.png')
saveas(gcf, 'treated_vs_untreated.fig')

%% compare exponential vs spent (1 minute frame rate)
cd([dirsave '/03262022_analysis']);

labels={'rich media', 'spent media', 'PBS'}; 

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24)
subplot(1,2,1), hold on
p1=meanPlot(untreated_100{2,2}, untreated_100{2, 5}, okabeIto{6}, okabeIto{6}, transparency)
p2=meanPlot(treated_100{7,2}, treated_100{7, 5}, okabeIto{1}, okabeIto{1}, transparency)
p3=meanPlot(PBS_100{3,2}, PBS_100{3, 5}, okabeIto{2}, okabeIto{2}, transparency)
ylim([0 Inf])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Time (minutes)')

hleg=legend([p1(1), p2(1), p3(1)], labels)
title(hleg,'Media')

subplot(1,2,2), hold on
p1=meanPlot(untreated_100{1,2}, untreated_100{1, 7}, okabeIto{6}, okabeIto{6}, transparency)
p2=meanPlot(treated_100{7,2}, treated_100{7, 7},okabeIto{1}, okabeIto{1}, transparency)
p3=meanPlot(PBS_100{3,2}, PBS_100{3, 7}, okabeIto{2}, okabeIto{2}, transparency)
ylim([0 Inf])
ylabel('Corrected Fluorescence (A.U.)')
xlabel('Time (minutes)')

hleg=legend([p1(1), p2(1), p3(1)], labels)
title(hleg,'Media')

saveas(gcf, 'rich_vs_spent.png')
saveas(gcf, 'rich_vs_spent.fig')

%% compare differences in growth rate

figure, hold on
plot(untreated_100{1,1}(1:width(untreated_100{1,10})), mean(untreated_100{1, 10}, 1, 'omitnan')) %, color_purple{1})
plot(untreated_100{2,1}(1:width(untreated_100{2,10})), mean(untreated_100{2, 10}, 1, 'omitnan')) %, color_purple{2})
plot(PBS_100{1,1}(1:width(PBS_100{1,8})), mean(PBS_100{1, 8}, 1, 'omitnan')) %, color_green{1})
plot(PBS_100{3,1}(1:width(PBS_100{3,8})), mean(PBS_100{3, 8}, 1, 'omitnan')) %, color_green{2})
plot(PBS_100{5,1}(1:width(PBS_100{5,8})), mean(PBS_100{5, 8}, 1, 'omitnan')) %, color_green{3})
plot(PBS_100{7,1}(1:width(PBS_100{7,8})), mean(PBS_100{7, 8}, 1, 'omitnan')) %, color_green{4})
% scatter(repelem(7, height(treated_100{1,8})), treated_100{1, 8}, [], colorcode{1})
% scatter(repelem(8, height(treated_100{2,8})), treated_100{2, 8}, [], colorcode{2})
% scatter(repelem(9, height(treated_100{3,8})), treated_100{3, 8}, [], colorcode{3})
% scatter(repelem(10, height(treated_100{4,8})), treated_100{4, 8}, [], colorcode{4})
% scatter(repelem(11, height(treated_100{5,8})), treated_100{5, 8}, [], colorcode{5})
plot(treated_100{7,1}(1:width(treated_100{7,8})), mean(treated_100{7, 8}, 1, 'omitnan')) %, colorcode{6})

figure, hold on
plot(untreated_100{1,1}(1:width(untreated_100{1,10})), untreated_100{1, 10}, '-r') 
plot(untreated_100{2,1}(1:width(untreated_100{2,10})), untreated_100{2, 10}, '-b') 
plot(PBS_100{7,1}(1:width(PBS_100{7,8})), PBS_100{7, 8}, '-g') 
plot(treated_100{7,1}(1:width(treated_100{7,8})), treated_100{7, 8}, '-k') 

%% determine the difference between fluorescence traces as a fuction of frame rate
% figure, hold on
% p1=meanPlot(untreated_100{1,1}, untreated_100{1, 3}, okabeIto{6}, okabeIto{6}, transparency)
% p2=meanPlot(controls_100{1,1}, controls_100{1, 3}, okabeIto{1}, okabeIto{1}, transparency)
% ylabel('Fluorescence (A.U.)')
% xlabel('Time (minutes)')
% 
% hleg=legend([p1(1), p2(1)], {'1 min', '1 s'})
% title(hleg,'Frame Rate')
% 
% dt1 = [1, diff(untreated_100{1,1}, 1, 2)];
% dt2 = [1, diff(controls_100{1,1}, 1, 2)];
% 
% figure, hold on
% p1=meanPlot(untreated_100{1,1}./dt1, untreated_100{1, 3}, okabeIto{6}, okabeIto{6}, transparency)
% p2=meanPlot(controls_100{1,1}./dt2, controls_100{1, 3}, okabeIto{1}, okabeIto{1}, transparency)
% ylabel('Fluorescence (A.U.)')
% xlabel('Frames')
% 
% hleg=legend([p1(1), p2(1)], {'1 min', '1 s'})
% title(hleg,'Frame Rate')

dt1 = [1, diff(untreated_100{1,2}, 1, 2)];
dt2 = [0.02, diff(controls_100{1,2}, 1, 2)];
dt3 = [0.0333, diff(controls_100{2,2}, 1, 2)];
dt4 = [0.05, diff(controls_100{3,2}, 1, 2)];
dt5 = [1, diff(untreated_100{4,2}, 1, 2)];

figure, hold on
p2=meanPlot(controls_100{1,2}./dt2, controls_100{1, 6}, okabeIto{2}, okabeIto{2}, transparency)
p3=meanPlot(controls_100{2,2}./dt3, controls_100{2, 6}, okabeIto{3}, okabeIto{3}, transparency)
p4=meanPlot(controls_100{3,2}./dt4, controls_100{3, 6}, okabeIto{4}, okabeIto{4}, transparency)
p1=meanPlot(untreated_100{1,2}./dt1, untreated_100{1, 5}, okabeIto{1}, okabeIto{1}, transparency)
p5=meanPlot(untreated_100{4,2}./dt5, untreated_100{4, 5}, okabeIto{5}, okabeIto{5}, transparency)
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Frames')

hleg=legend([p2(1), p3(1), p4(1), p1(1), p5(1)], {'1 s', '2 s', '3 s', '1 min', '10 min'})
title(hleg,'Frame Rate')

figure, hold on
p2=meanPlot(controls_100{1,2}./dt2, controls_100{1, 6}, okabeIto{2}, okabeIto{2}, transparency)
p3=meanPlot(controls_100{2,2}./dt3, controls_100{2, 6}, okabeIto{3}, okabeIto{3}, transparency)
p4=meanPlot(controls_100{3,2}./dt4, controls_100{3, 6}, okabeIto{4}, okabeIto{4}, transparency)
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Frames')
ylim([0 Inf])

hleg=legend([p2(1), p3(1), p4(1)], {'1 s', '2 s', '3 s'})
title(hleg,'Frame Rate')

%% fit the normalized fluor. vs frame rate traces to an exponential
dt1 = [1, diff(untreated_100{1,2}, 1, 2)];
dt2 = [0.02, diff(controls_100{1,2}, 1, 2)];
dt3 = [0.0333, diff(controls_100{2,2}, 1, 2)];
dt4 = [0.05, diff(controls_100{3,2}, 1, 2)];
dt5 = [1, diff(untreated_100{4,2}, 1, 2)];

untreated_100{1,11} = untreated_100{1,2}./dt1;
controls_100{1,12} = controls_100{1,2}./dt2;
controls_100{2,12} = controls_100{2,2}./dt3;
controls_100{3,12} = controls_100{3,2}./dt4;
untreated_100{4,11} = untreated_100{4,2}./dt5;

[untreated_100{1,12}, untreated_100{1,13}]=rhoCalc(untreated_100{1,11}, untreated_100{1, 5}, 1);
[controls_100{1,13}, controls_100{1,14}]=rhoCalc(controls_100{1,12}, controls_100{1, 6}, 1);
[controls_100{2,13}, controls_100{2,14}]=rhoCalc(controls_100{2,12}, controls_100{2, 6}, 1);
[controls_100{3,13}, controls_100{3,14}]=rhoCalc(controls_100{3,12}, controls_100{3, 6}, 1);
[untreated_100{4,12}, untreated_100{4,13}]=rhoCalc(untreated_100{4, 11}, untreated_100{4, 5}, 1);

%% plot the fit of the exponential function
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',16)
subplot(1,4,1)
plot(controls_100{1,12}, controls_100{1, 6}, 'Color', okabeIto{2}), hold on
plot(controls_100{1,12}, controls_100{1, 14}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Frames')
title('100% Intensity, 1 s Frame Rate')
% saveas(gcf, 'fit100_1s.png')
% saveas(gcf, 'fit100_1s.fig')

subplot(1,4,2)
plot(controls_100{2,12}, controls_100{2, 6}, 'Color', okabeIto{3}), hold on
plot(controls_100{2,12}, controls_100{2, 14}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Frames')
title('100% Intensity, 2 s Frame Rate')
% saveas(gcf, 'fit100_2s.png')
% saveas(gcf, 'fit100_2s.fig')

subplot(1,4,3)
plot(controls_100{3,12}, controls_100{3, 6}, 'Color', okabeIto{4}), hold on
plot(controls_100{3,12}, controls_100{3, 14}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Frames')
title('100% Intensity, 3 s Frame Rate')
% saveas(gcf, 'fit100_3s.png')
% saveas(gcf, 'fit100_3s.fig')

%figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',16)
subplot(1,4,4)
plot(untreated_100{1,11}, untreated_100{1, 5}, 'Color', okabeIto{1}), hold on
plot(untreated_100{1,11}, untreated_100{1, 13}, '--k')
ylim([0 1.1])
ylabel('Normalized Fluorescence (A.U.)')
xlabel('Frames')
title('100% Intensity, 1 min Frame Rate')

% subplot(1,2,2)
% plot(untreated_100{4,11}, untreated_100{4, 5}, 'Color', color_blue{2}), hold on
% plot(untreated_100{4,11}, untreated_100{4, 13}, '--k')
% ylim([0 1.1])
% ylabel('Normalized Fluorescence (A.U.)')
% xlabel('Frame Rate')
% title('100% Intensity, 10 min Frame Rate')

%% plot rho as a function of frame rate
rho={controls_100{1, 13}; controls_100{2, 13}; controls_100{3, 13}; untreated_100{1,12}};
dt=[1.2/60, 2/60, 3/60, 1]; 

% linearCoef3 = polyfit([repelem(dt(1), length(rho{1})), repelem(dt(2), length(rho{2})), repelem(dt(3), length(rho{3}))],[rho{1}', rho{2}', rho{3}'],1);
% linearFit3= polyval(linearCoef3,[0 1.2/60 2/60 3/60]);

cd([dirsave '/03262022_analysis']);

rho_means = cellfun(@(x)mean(x, 1, 'omitnan'), rho);
rho_std = cellfun(@(x)std(x, 0, 1, 'omitnan'), rho);

%labels={'\tau 100% power', '\tau 20% power', 'mean \tau', 'line of best fit'};

figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'DefaultAxesFontSize',24), hold on
p1=scatter(repelem(dt(1), length(rho{1})), rho{1}, 'MarkerFaceColor', '#56B4E9', 'MarkerEdgeColor', '#56B4E9')
scatter(repelem(dt(2), length(rho{2})), rho{2}, 'MarkerFaceColor', '#009E73', 'MarkerEdgeColor', '#009E73')
scatter(repelem(dt(3), length(rho{3})), rho{3}, 'MarkerFaceColor', '#F0E442', 'MarkerEdgeColor', '#F0E442')
scatter(repelem(dt(4), length(rho{4})), rho{4}, 'MarkerFaceColor', '#E69F00', 'MarkerEdgeColor', '#E69F00')

p3=scatter(dt, rho_means, 'MarkerFaceColor', 'black')
errorbar(dt(1), rho_means(1), rho_std(1), 'Color', 'black')
errorbar(dt(2), rho_means(2), rho_std(2), 'Color', 'black')
errorbar(dt(3), rho_means(3), rho_std(3), 'Color', 'black')
errorbar(dt(4), rho_means(4), rho_std(4), 'Color', 'black')

% p4=plot([0, dt(1), dt(2), dt(3)], linearFit1, '--k', 'LineWidth', 1)
% plot([0, dt(4), dt(5), dt(6)], linearFit2, '--k', 'LineWidth', 1)
xlim([0, 1.1])
ylim([0 40])
ylabel('\rho')
%hleg=legend([p1(1), p2(1), p3, p4], labels, 'Location', 'southeast')
%xticks([0, dt(1),  dt(4), dt(2),  dt(5), dt(3)])
%xticklabels({'0', '1.2 s', '1.76 s', '2 s', '2.3 s', '3 s'})
xlabel('Frame Rate (minutes)')

% caption1 = sprintf('%f * frame rate + %f', linearCoef1(1), linearCoef1(2));
% caption2 = sprintf('%f * frame rate + %f', linearCoef2(1), linearCoef2(2));
% text(dt(2), tau_means(2)+1.5, ['\tau = ' caption1], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');
% text(dt(1), tau_means(5)+1.5, ['\tau = ' caption2], 'FontSize', 18, 'Color', 'k', 'FontWeight', 'bold');

title('Rho vs Frame Rate')
%% Functions
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

function [rho, yhat]=rhoCalc(tme, normintensity, rho0)
    
    [nrow, ncol]=size(normintensity);
    
    idx=find(sum(isnan(normintensity), 2) < ncol-1);
    normintensity=normintensity(idx, :);
    
    [nrow, ncol]=size(normintensity);
    rho=nan(nrow,1);
    yhat=nan(size(normintensity));
    modelfun = @(rho, x)exp(-x*rho);
    for i=1:nrow
        rho(i, 1)=nlinfit(tme, normintensity(i,:), modelfun, rho0);
        yhat(i, :)=modelfun(rho(i), tme);
    end
    
end


%to find the time point at which the corrected trace autofluorescences and calculate
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
%https://www.doitpoms.ac.uk/glossary/entry.php?term=strain
function [strain]=strainCalc(lcell, imstart)

    %negative strain=compression
    %positive strain=tensile
    %divide by unloaded dimension

    dl=lcell(:,imstart+1)-lcell(:, imstart);
    strain=dl./lcell(:, imstart+1);


end

%to calculation growth rate
function [growthRate]=gRateCalc(time, lcell, imstart)

    lcell=lcell(:, 1:imstart+3);
    time=time(1:imstart+3);

%     dl=diff(lcell, 1, 2);
%     dt=diff(time, 1, 2);
% 
%     sl=(lcell(:, 1:end-1)+(dl+lcell(:, 1:end-1)))/2;
% 
%     growthRate=(dl./dt)./sl;

    %Calculate the growth rate
    deltat=time(2:end)-time(1:end-1);
    v=(lcell(:,2:end)-lcell(:,1:end-1))./((lcell(:,1:end-1)+lcell(:,2:end))/2);
    for i=1:height(v)
        v(i,:)=v(i,:)./deltat;
    end

    growthRate=v;
    
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
function [p]=meanPlot(tme, normintensity, colorcode1, colorcode2, transparency)

[nrow, ~]=size(normintensity);

idxNaN=isnan(normintensity);
idx=find(sum(idxNaN, 1)<nrow/2);
idx=max(idx);

normintensity=normintensity(:, 1:idx);
tme=tme(:, 1:idx);

nmean=mean(normintensity, 1, 'omitnan');
nstd=std(normintensity, 0, 1, 'omitnan');

if nargin==5
    ciplot(nmean-nstd, nmean+nstd, tme, colorcode2, transparency)
    p=plot(tme, nmean, 'Color', colorcode1, 'LineWidth', 1.5)
elseif nargin==4
    p=errorbar(tme, nmean, -nstd, nstd, 'Color', colorcode1, 'LineWidth', 1.5)
elseif nargin==3
    p=plot(tme, nmean, 'Color', colorcode1, 'LineWidth', 1.5)
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

