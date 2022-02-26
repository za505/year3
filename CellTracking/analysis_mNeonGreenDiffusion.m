%Author: Zarina Akbary
%Date: 02/23/2022
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
% 
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
% alpha=152.6455;
% intercept=-1.0976;

%% User Input 

%color codes must be in RGB [0-1] format to be used in ciplot
colorcode={[204 0 0], [204 102 0], [204 204 0], [102 204 0], [0 204 204], [0 0 204], [102 0 204], [204 0 204], [204 0 102], [255 102 102], [255 178 102], [102 255 102], [102 255 255], [102 178 255], [178 102 255], [255 102 255],[255 102 178]};
colorcode2={[255 51 51], [255 153 51], [255 255 51], [153 255 51], [51 255 255], [51 51 255], [153 51 255], [255 51 255], [255 51 153], [255 204 204], [255 229 204], [204 255 204], [204 255 255], [204 229 255], [229 204 255], [255 204 255],[255 204 229]};

colorcode=cellfun(@(x)(x./255), colorcode, 'UniformOutput', false);
colorcode2=cellfun(@(x)(x./255), colorcode2, 'UniformOutput', false);

transparency = 0.3; %this is the alpha argument for the ciplot function

%location and names of of the *.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis';
basenames={'10232021_Exp1', '10262021_Exp1', '02212022_Exp2', '02122022_Exp1','02122022_Exp2', '02092022_Exp1', '02212022_Exp1', '11192021_Exp2', '12082021_Exp3', '11192021_Exp1', '11302021_Exp1', '10232021_Exp2', '10262021_Exp2', '01142022_Exp1', '01172022_Exp1', '01172022_Exp2', '01242022_Exp1', '01262022_Exp1', '02122022_Exp3', '02192022_Exp1','11202021_Exp1', '12082021_Exp1', '12082021_Exp2', '02192022_Exp2', '02192022_Exp3', '02192022_Exp4'};

%% normalize data
imstarts=[6, 5, 6, 6, 6, 6, 6, 10, 10, 4, 4, 3, 5, 4, 6, 6, 3, 5, 6, 6, NaN, NaN, NaN, NaN, NaN, NaN];

for i=1:length(basenames)
    basename=basenames{i};
    cd([dirsave '/rawFiles'])
    
    datadir=dir([basename '*']);
    imstart=imstarts(i);
    
    if ismember(i, [21:26]) %controls
        
        intensity=[];
        bgintensity=[];
        adjintensity=[];
        normintensity=[];
        lcell=[];
        
        for j=1:height(datadir)
            cd([dirsave '/rawFiles']);
            load(datadir(j).name);
            
            
            [intensity1, bgintensity1, adjintensity1, normintensity1, lCell, time, tme, imstart]=controlNormalize(datadir);
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
            
            
            [intensity1, adjintensity1, normintensity1, lCell, time, tme, imstart]=dataNormalize(datadir, imstart);
            intensity=[intensity;intensity1];
            adjintensity=[adjintensity; adjintensity1];
            normintensity=[normintensity; normintensity1];
            lcell=[lcell; lCell];
            
        end
        
        cd([dirsave '/normalizedFiles'])
        save([basename '_norm.mat'], 'intensity', 'adjintensity', 'normintensity', 'lcell', 'imstart', 'time', 'tme');
    end
    
end


%% correct data
for i=1:length(basenames)
    
    basename=basenames{i};
    cd([dirsave '/normalizedFiles'])
    datadir=dir([basename '*']);
    
    if ismember(i, [3, 7, 24:26]) % 20% intensity
        
        %for photobleach correction (20% intensity)
        alpha=152.6455;
        intercept=-1.0976;

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

%% compare the control traces
%each row is a different frame rate (1, 2, and 3 second), each column is a
%different variable (time, tme, intensity, bgintensity, adjintensity,
%normintensity, lcell, and Cnew)

controls_100 = cell(3, 8);
controls_20 = cell(3, 8);

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
        
        cd([dirsave '/correctedFiles'])
        load([basename '_corrected.mat'])
        controls_20{idx2, 8}=Cnew;
        
        idx2=idx2+1;
    else
        continue
    end
    
end

%% plot to compare the controls
figure, hold on
%ciplot(mean(controls_100{1,8}, 1, 'omitnan')-std(controls_100{1,8}, 0, 1, 'omitnan'), mean(controls_100{1,8}, 1, 'omitnan')+std(controls_100{1,8}, 0, 1, 'omitnan'), controls_100{1,2}, colorcode2{1}, transparency)
plot(controls_100{1,2}, mean(controls_100{1,8}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(controls_100{2,8}, 1, 'omitnan')-std(controls_100{2,8}, 0, 1, 'omitnan'), mean(controls_100{2,8}, 1, 'omitnan')+std(controls_100{2,8}, 0, 1, 'omitnan'), controls_100{2,2}, colorcode2{2}, transparency)
plot(controls_100{2,2}, mean(controls_100{2,8}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)

%ciplot(mean(controls_100{3,8}, 1, 'omitnan')-std(controls_100{3,8}, 0, 1, 'omitnan'), mean(controls_100{3,8}, 1, 'omitnan')+std(controls_100{3,8}, 0, 1, 'omitnan'), controls_100{3,2}, colorcode2{4}, transparency)
plot(controls_100{3,2}, mean(controls_100{3,8}, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

%ciplot(mean(controls_20{1,8}, 1, 'omitnan')-std(controls_20{1,8}, 0, 1, 'omitnan'), mean(controls_20{1,8}, 1, 'omitnan')+std(controls_20{1,8}, 0, 1, 'omitnan'), controls_20{1,2}, colorcode2{5}, transparency)
plot(controls_20{1,2}, mean(controls_20{1,8}, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(controls_20{2,8}, 1, 'omitnan')-std(controls_20{2,8}, 0, 1, 'omitnan'), mean(controls_20{2,8}, 1, 'omitnan')+std(controls_20{2,8}, 0, 1, 'omitnan'), controls_20{2,2}, colorcode2{6}, transparency)
plot(controls_20{2,2}, mean(controls_20{2,8}, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

%ciplot(mean(controls_20{3,8}, 1, 'omitnan')-std(controls_20{3,8}, 0, 1, 'omitnan'), mean(controls_20{3,8}, 1, 'omitnan')+std(controls_20{3,8}, 0, 1, 'omitnan'), controls_20{3,2}, colorcode2{7}, transparency)
plot(controls_20{3,2}, mean(controls_20{3,8}, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

%legend({'1 s, 100% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', '', '1 s, 20% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', ''})
legend({'1 s, 100% intensity, average', '2 s, 100% intensity, average', '3 s, 100% intensity, average', '1 s, 20% intensity, average', '2 s, 20% intensity, average', '3 s, 20% intensity, average'}, 'location', 'east')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
title('Control Averaged Corrected Traces')
% %saveas(gcf, 'frameRate_correctedAvg.png');

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

%% Compare the 20 minute frame rate data
figure, hold on
%ciplot(mean(untreated_100{1,3}, 1, 'omitnan')-std(untreated_100{1,3}, 0, 1, 'omitnan'), mean(untreated_100{1,3}, 1, 'omitnan')+std(untreated_100{1,3}, 0, 1, 'omitnan'), untreated_100{1,1}, colorcode2{1}, transparency)
plot(untreated_100{5,1}, mean(untreated_100{5,3}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
%ciplot(mean(untreated_100{2,3}, 1, 'omitnan')-std(untreated_100{2,3}, 0, 1, 'omitnan'), mean(untreated_100{2,3}, 1, 'omitnan')+std(untreated_100{2,3}, 0, 1, 'omitnan'), untreated_100{2,1}, colorcode2{2}, transparency)
plot(untreated_20{2,1}, mean(untreated_20{2,3}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)
title('Raw Intensity Plots')

figure, hold on
%ciplot(mean(untreated_100{1,3}, 1, 'omitnan')-std(untreated_100{1,3}, 0, 1, 'omitnan'), mean(untreated_100{1,3}, 1, 'omitnan')+std(untreated_100{1,3}, 0, 1, 'omitnan'), untreated_100{1,1}, colorcode2{1}, transparency)
plot(untreated_100{5,1}, mean(untreated_100{5,4}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
%ciplot(mean(untreated_100{2,3}, 1, 'omitnan')-std(untreated_100{2,3}, 0, 1, 'omitnan'), mean(untreated_100{2,3}, 1, 'omitnan')+std(untreated_100{2,3}, 0, 1, 'omitnan'), untreated_100{2,1}, colorcode2{2}, transparency)
plot(untreated_20{2,1}, mean(untreated_20{2,4}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)
title('Adjusted Intensity Plots')

figure, hold on
%ciplot(mean(untreated_100{1,3}, 1, 'omitnan')-std(untreated_100{1,3}, 0, 1, 'omitnan'), mean(untreated_100{1,3}, 1, 'omitnan')+std(untreated_100{1,3}, 0, 1, 'omitnan'), untreated_100{1,1}, colorcode2{1}, transparency)
plot(untreated_100{5,2}, mean(untreated_100{5,5}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
%ciplot(mean(untreated_100{2,3}, 1, 'omitnan')-std(untreated_100{2,3}, 0, 1, 'omitnan'), mean(untreated_100{2,3}, 1, 'omitnan')+std(untreated_100{2,3}, 0, 1, 'omitnan'), untreated_100{2,1}, colorcode2{2}, transparency)
plot(untreated_20{2,2}, mean(untreated_20{2,5}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)
title('Normalized Intensity Plots')


%normalize the adjusted trace
imstart=untreated_100{5,8};
sc100=mean(untreated_100{5,4}, 1, 'omitnan');
sc100=sc100./sc100(:, imstart);
sc100=sc100(imstart:end);
time100=untreated_100{5,1}(imstart:end);

imstart=untreated_20{2,8};
sc20=mean(untreated_20{2,4}, 1, 'omitnan');
sc20=sc20./sc20(:, imstart);
sc20=sc20(imstart:end);
time20=untreated_20{2,1}(imstart:end);

figure, hold on
plot(time100, sc100, 'Color', colorcode{1}, 'LineWidth', 1)
plot(time20, sc20, 'Color', colorcode{2}, 'LineWidth', 1)
title('Mean Normalized Intensity Plots')

pause
%% plot to compare the untreated experiments
cd(dirsave)

%compare raw traces
figure, hold on
%ciplot(mean(untreated_100{1,3}, 1, 'omitnan')-std(untreated_100{1,3}, 0, 1, 'omitnan'), mean(untreated_100{1,3}, 1, 'omitnan')+std(untreated_100{1,3}, 0, 1, 'omitnan'), untreated_100{1,1}, colorcode2{1}, transparency)
plot(untreated_100{1,1}, mean(untreated_100{1,3}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(untreated_100{2,3}, 1, 'omitnan')-std(untreated_100{2,3}, 0, 1, 'omitnan'), mean(untreated_100{2,3}, 1, 'omitnan')+std(untreated_100{2,3}, 0, 1, 'omitnan'), untreated_100{2,1}, colorcode2{2}, transparency)
plot(untreated_100{2,1}, mean(untreated_100{2,3}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)

%ciplot(mean(untreated_100{3,3}, 1, 'omitnan')-std(untreated_100{3,3}, 0, 1, 'omitnan'), mean(untreated_100{3,3}, 1, 'omitnan')+std(untreated_100{3,3}, 0, 1, 'omitnan'), untreated_100{3,1}, colorcode2{4}, transparency)
plot(untreated_100{3,1}, mean(untreated_100{3,3}, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

%ciplot(mean(untreated_100{4,3}, 1, 'omitnan')-std(untreated_100{4,3}, 0, 1, 'omitnan'), mean(untreated_100{4,3}, 1, 'omitnan')+std(untreated_100{4,3}, 0, 1, 'omitnan'), untreated_100{4,1}, colorcode2{5}, transparency)
plot(untreated_100{4,1}, mean(untreated_100{4,3}, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(untreated_100{5,3}, 1, 'omitnan')-std(untreated_100{5,3}, 0, 1, 'omitnan'), mean(untreated_100{5,3}, 1, 'omitnan')+std(untreated_100{5,3}, 0, 1, 'omitnan'), untreated_100{5,1}, colorcode2{6}, transparency)
plot(untreated_100{5,1}, mean(untreated_100{5,3}, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

%ciplot(mean(untreated_20{1,3}, 1, 'omitnan')-std(untreated_20{1,3}, 0, 1, 'omitnan'), mean(untreated_20{1,3}, 1, 'omitnan')+std(untreated_20{1,3}, 0, 1, 'omitnan'), untreated_20{1,1}, colorcode2{7}, transparency)
plot(untreated_20{1,1}, mean(untreated_20{1,3}, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

%ciplot(mean(untreated_20{2,3}, 1, 'omitnan')-std(untreated_20{2,3}, 0, 1, 'omitnan'), mean(untreated_20{2,3}, 1, 'omitnan')+std(untreated_20{2,3}, 0, 1, 'omitnan'), untreated_20{2,1}, colorcode2{8}, transparency)
plot(untreated_20{2,1}, mean(untreated_20{2,3}, 1, 'omitnan'), 'Color', colorcode{8}, 'LineWidth', 1)

%legend({'1 s, 100% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', '', '1 s, 20% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', ''})
legend({'1 min, 100% intensity, average', '1 min, 100% intensity, average', '5 min, 100% intensity, average', '10 min, 100% intensity, average', '20 min, 100% intensity, average', '1 min, 20% intensity, average', '20 min, 20% intensity, average'}, 'location', 'northeast')
xlabel('Time (minutes)')
ylabel('Fluorescence')
title('Untreated Averaged Traces')
saveas(gcf, 'untreated_rawAvg.png');
saveas(gcf, 'untreated_rawAvg.fig');

%compare adjusted traces
figure, hold on
%ciplot(mean(untreated_100{1,4}, 1, 'omitnan')-std(untreated_100{1,4}, 0, 1, 'omitnan'), mean(untreated_100{1,4}, 1, 'omitnan')+std(untreated_100{1,4}, 0, 1, 'omitnan'), untreated_100{1,1}, colorcode2{1}, transparency)
plot(untreated_100{1,1}, mean(untreated_100{1,4}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(untreated_100{2,4}, 1, 'omitnan')-std(untreated_100{2,4}, 0, 1, 'omitnan'), mean(untreated_100{2,4}, 1, 'omitnan')+std(untreated_100{2,4}, 0, 1, 'omitnan'), untreated_100{2,1}, colorcode2{2}, transparency)
plot(untreated_100{2,1}, mean(untreated_100{2,4}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)

%ciplot(mean(untreated_100{3,4}, 1, 'omitnan')-std(untreated_100{3,4}, 0, 1, 'omitnan'), mean(untreated_100{3,4}, 1, 'omitnan')+std(untreated_100{3,4}, 0, 1, 'omitnan'), untreated_100{3,1}, colorcode2{4}, transparency)
plot(untreated_100{3,1}, mean(untreated_100{3,4}, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

%ciplot(mean(untreated_100{4,4}, 1, 'omitnan')-std(untreated_100{4,4}, 0, 1, 'omitnan'), mean(untreated_100{4,4}, 1, 'omitnan')+std(untreated_100{4,4}, 0, 1, 'omitnan'), untreated_100{4,1}, colorcode2{5}, transparency)
plot(untreated_100{4,1}, mean(untreated_100{4,4}, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(untreated_100{5,4}, 1, 'omitnan')-std(untreated_100{5,4}, 0, 1, 'omitnan'), mean(untreated_100{5,4}, 1, 'omitnan')+std(untreated_100{5,4}, 0, 1, 'omitnan'), untreated_100{5,1}, colorcode2{6}, transparency)
plot(untreated_100{5,1}, mean(untreated_100{5,4}, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

%ciplot(mean(untreated_20{1,4}, 1, 'omitnan')-std(untreated_20{1,4}, 0, 1, 'omitnan'), mean(untreated_20{1,4}, 1, 'omitnan')+std(untreated_20{1,4}, 0, 1, 'omitnan'), untreated_20{1,1}, colorcode2{7}, transparency)
plot(untreated_20{1,1}, mean(untreated_20{1,4}, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

%ciplot(mean(untreated_20{2,4}, 1, 'omitnan')-std(untreated_20{2,4}, 0, 1, 'omitnan'), mean(untreated_20{2,4}, 1, 'omitnan')+std(untreated_20{2,4}, 0, 1, 'omitnan'), untreated_20{2,1}, colorcode2{8}, transparency)
plot(untreated_20{2,1}, mean(untreated_20{2,4}, 1, 'omitnan'), 'Color', colorcode{8}, 'LineWidth', 1)

%legend({'1 s, 100% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', '', '1 s, 20% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', ''})
legend({'1 min, 100% intensity, average', '1 min, 100% intensity, average', '5 min, 100% intensity, average', '10 min, 100% intensity, average', '20 min, 100% intensity, average', '1 min, 20% intensity, average', '20 min, 20% intensity, average'}, 'location', 'northeast')
xlabel('Time (minutes)')
ylabel('Adjusted Fluorescence')
title('Untreated Averaged Adjusted Traces')
saveas(gcf, 'untreated_adjAvg.png');
saveas(gcf, 'untreated_adjAvg.fig');

%compare normalized traces
figure, hold on
%ciplot(mean(untreated_100{1,5}, 1, 'omitnan')-std(untreated_100{1,5}, 0, 1, 'omitnan'), mean(untreated_100{1,5}, 1, 'omitnan')+std(untreated_100{1,5}, 0, 1, 'omitnan'), untreated_100{1,2}, colorcode2{1}, transparency)
plot(untreated_100{1,2}, mean(untreated_100{1,5}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(untreated_100{2,5}, 1, 'omitnan')-std(untreated_100{2,5}, 0, 1, 'omitnan'), mean(untreated_100{2,5}, 1, 'omitnan')+std(untreated_100{2,5}, 0, 1, 'omitnan'), untreated_100{2,2}, colorcode2{2}, transparency)
plot(untreated_100{2,2}, mean(untreated_100{2,5}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)

%ciplot(mean(untreated_100{3,5}, 1, 'omitnan')-std(untreated_100{3,5}, 0, 1, 'omitnan'), mean(untreated_100{3,5}, 1, 'omitnan')+std(untreated_100{3,5}, 0, 1, 'omitnan'), untreated_100{3,2}, colorcode2{4}, transparency)
plot(untreated_100{3,2}, mean(untreated_100{3,5}, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

%ciplot(mean(untreated_100{4,5}, 1, 'omitnan')-std(untreated_100{4,5}, 0, 1, 'omitnan'), mean(untreated_100{4,5}, 1, 'omitnan')+std(untreated_100{4,5}, 0, 1, 'omitnan'), untreated_100{4,2}, colorcode2{5}, transparency)
plot(untreated_100{4,2}, mean(untreated_100{4,5}, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(untreated_100{5,5}, 1, 'omitnan')-std(untreated_100{5,5}, 0, 1, 'omitnan'), mean(untreated_100{5,5}, 1, 'omitnan')+std(untreated_100{5,5}, 0, 1, 'omitnan'), untreated_100{5,2}, colorcode2{6}, transparency)
plot(untreated_100{5,2}, mean(untreated_100{5,5}, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

%ciplot(mean(untreated_20{1,5}, 1, 'omitnan')-std(untreated_20{1,5}, 0, 1, 'omitnan'), mean(untreated_20{1,5}, 1, 'omitnan')+std(untreated_20{1,5}, 0, 1, 'omitnan'), untreated_20{1,2}, colorcode2{7}, transparency)
plot(untreated_20{1,2}, mean(untreated_20{1,5}, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

%ciplot(mean(untreated_20{2,5}, 1, 'omitnan')-std(untreated_20{2,5}, 0, 1, 'omitnan'), mean(untreated_20{2,5}, 1, 'omitnan')+std(untreated_20{2,5}, 0, 1, 'omitnan'), untreated_20{2,2}, colorcode2{8}, transparency)
plot(untreated_20{2,2}, mean(untreated_20{2,5}, 1, 'omitnan'), 'Color', colorcode{8}, 'LineWidth', 1)

%legend({'1 s, 100% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', '', '1 s, 20% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', ''})
legend({'1 min, 100% intensity, average', '1 min, 100% intensity, average', '5 min, 100% intensity, average', '10 min, 100% intensity, average', '20 min, 100% intensity, average', '1 min, 20% intensity, average', '20 min, 20% intensity, average'}, 'location', 'northeast')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
title('Untreated Averaged Normalized Traces')
saveas(gcf, 'untreated_normAvg.png');
saveas(gcf, 'untreated_normAvg.fig');

%compare corrected traces
figure, hold on
%ciplot(mean(untreated_100{1,7}, 1, 'omitnan')-std(untreated_100{1,7}, 0, 1, 'omitnan'), mean(untreated_100{1,7}, 1, 'omitnan')+std(untreated_100{1,7}, 0, 1, 'omitnan'), untreated_100{1,2}, colorcode2{1}, transparency)
plot(untreated_100{1,2}, mean(untreated_100{1,7}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(untreated_100{2,7}, 1, 'omitnan')-std(untreated_100{2,7}, 0, 1, 'omitnan'), mean(untreated_100{2,7}, 1, 'omitnan')+std(untreated_100{2,7}, 0, 1, 'omitnan'), untreated_100{2,2}, colorcode2{2}, transparency)
plot(untreated_100{2,2}, mean(untreated_100{2,7}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)

%ciplot(mean(untreated_100{3,7}, 1, 'omitnan')-std(untreated_100{3,7}, 0, 1, 'omitnan'), mean(untreated_100{3,7}, 1, 'omitnan')+std(untreated_100{3,7}, 0, 1, 'omitnan'), untreated_100{3,2}, colorcode2{4}, transparency)
plot(untreated_100{3,2}, mean(untreated_100{3,7}, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

%ciplot(mean(untreated_100{4,7}, 1, 'omitnan')-std(untreated_100{4,7}, 0, 1, 'omitnan'), mean(untreated_100{4,7}, 1, 'omitnan')+std(untreated_100{4,7}, 0, 1, 'omitnan'), untreated_100{4,2}, colorcode2{5}, transparency)
plot(untreated_100{4,2}, mean(untreated_100{4,7}, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(untreated_100{5,7}, 1, 'omitnan')-std(untreated_100{5,7}, 0, 1, 'omitnan'), mean(untreated_100{5,7}, 1, 'omitnan')+std(untreated_100{5,7}, 0, 1, 'omitnan'), untreated_100{5,2}, colorcode2{6}, transparency)
plot(untreated_100{5,2}, mean(untreated_100{5,7}, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

%ciplot(mean(untreated_20{1,7}, 1, 'omitnan')-std(untreated_20{1,7}, 0, 1, 'omitnan'), mean(untreated_20{1,7}, 1, 'omitnan')+std(untreated_20{1,7}, 0, 1, 'omitnan'), untreated_20{1,2}, colorcode2{7}, transparency)
plot(untreated_20{1,2}, mean(untreated_20{1,7}, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

%ciplot(mean(untreated_20{2,7}, 1, 'omitnan')-std(untreated_20{2,7}, 0, 1, 'omitnan'), mean(untreated_20{2,7}, 1, 'omitnan')+std(untreated_20{2,7}, 0, 1, 'omitnan'), untreated_20{2,2}, colorcode2{8}, transparency)
plot(untreated_20{2,2}, mean(untreated_20{2,7}, 1, 'omitnan'), 'Color', colorcode{8}, 'LineWidth', 1)

%legend({'1 s, 100% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', '', '1 s, 20% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', ''})
legend({'1 min, 100% intensity, average', '1 min, 100% intensity, average', '5 min, 100% intensity, average', '10 min, 100% intensity, average', '20 min, 100% intensity, average', '1 min, 20% intensity, average', '20 min, 20% intensity, average'}, 'location', 'northeast')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
title('Untreated Averaged Corrected Traces')
saveas(gcf, 'untreated_correctedAvg.png');
saveas(gcf, 'untreated_correctedAvg.fig');

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

%% plot to compare the PBS experiments
cd(dirsave)

figure, hold on
%ciplot(mean(PBS_100{1,7}, 1, 'omitnan')-std(PBS_100{1,7}, 0, 1, 'omitnan'), mean(PBS_100{1,7}, 1, 'omitnan')+std(PBS_100{1,7}, 0, 1, 'omitnan'), PBS_100{1,2}, colorcode2{1}, transparency)
plot(PBS_100{1,2}, mean(PBS_100{1,7}, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(PBS_100{2,7}, 1, 'omitnan')-std(PBS_100{2,7}, 0, 1, 'omitnan'), mean(PBS_100{2,7}, 1, 'omitnan')+std(PBS_100{2,7}, 0, 1, 'omitnan'), PBS_100{2,2}, colorcode2{2}, transparency)
plot(PBS_100{2,2}, mean(PBS_100{2,7}, 1, 'omitnan'), 'Color', colorcode{2}, 'LineWidth', 1)

%ciplot(mean(PBS_100{3,7}, 1, 'omitnan')-std(PBS_100{3,7}, 0, 1, 'omitnan'), mean(PBS_100{3,7}, 1, 'omitnan')+std(PBS_100{3,7}, 0, 1, 'omitnan'), PBS_100{3,2}, colorcode2{4}, transparency)
plot(PBS_100{3,2}, mean(PBS_100{3,7}, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

%ciplot(mean(PBS_100{4,7}, 1, 'omitnan')-std(PBS_100{4,7}, 0, 1, 'omitnan'), mean(PBS_100{4,7}, 1, 'omitnan')+std(PBS_100{4,7}, 0, 1, 'omitnan'), PBS_100{4,2}, colorcode2{5}, transparency)
plot(PBS_100{4,2}(:, 1:60), mean(PBS_100{4,7}(:, 1:60), 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(PBS_100{5,7}, 1, 'omitnan')-std(PBS_100{5,7}, 0, 1, 'omitnan'), mean(PBS_100{5,7}, 1, 'omitnan')+std(PBS_100{5,7}, 0, 1, 'omitnan'), PBS_100{5,2}, colorcode2{6}, transparency)
plot(PBS_100{5,2}, mean(PBS_100{5,7}, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

%ciplot(mean(PBS_100{6,7}, 1, 'omitnan')-std(PBS_100{6,7}, 0, 1, 'omitnan'), mean(PBS_100{6,7}, 1, 'omitnan')+std(PBS_100{6,7}, 0, 1, 'omitnan'), PBS_100{6,2}, colorcode2{7}, transparency)
plot(PBS_100{6,2}, mean(PBS_100{6,7}, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

%ciplot(mean(PBS_100{7,7}, 1, 'omitnan')-std(PBS_100{7,7}, 0, 1, 'omitnan'), mean(PBS_100{7,7}, 1, 'omitnan')+std(PBS_100{7,7}, 0, 1, 'omitnan'), PBS_100{7,2}, colorcode2{8}, transparency)
plot(PBS_100{7,2}, mean(PBS_100{7,7}, 1, 'omitnan'), 'Color', colorcode{8}, 'LineWidth', 1)

%legend({'1 s, 100% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', '', '1 s, 20% intensity, average', '', '2 s, 100% intensity, average', '', '3 s, 100% intensity, average', ''})
legend({'PBS 2 min, average', 'PBS 2 min, average', 'PBS 20 min, average', 'PBS 20 min, average', 'PBS 60 min, average', 'PBS 60 min, average', 'PBS 120 min, average'}, 'location', 'northeast')
ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
title('PBS-incubated Averaged Corrected Traces')
saveas(gcf, 'PBS_correctedAvg.png');
saveas(gcf, 'PBS_correctedAvg.fig');

%% close all plots
close all

%% find where the mean traces plateau
avg_100=cell(height(untreated_100),1);
imend_100=nan(height(untreated_100),1);

avg_20=cell(height(untreated_20),1);
imend_20=nan(height(untreated_20),1);

tme_100 = [];
tme_20 = [];

for i=1:height(untreated_100)
    [avg_100{i}, imend_100(i)]=traceAvg(untreated_100{i, 7});
    idx=[1:imend_100(i)];
    tme_100=union(tme_100, untreated_100{i,2}(idx));
end

for i=1:height(untreated_20)
    [avg_20{i}, imend_20(i)]=traceAvg(untreated_20{i, 7});
    idx=[1:imend_20(i)];
    tme_20=union(tme_20, untreated_20{i,2}(idx));
end

%create a new matrix to combine the mean traces
ftrace_100 = nan(height(untreated_100), length(tme_100));
ftrace_20 = nan(height(untreated_20), length(tme_20));

for i=1:height(untreated_100)
    [~, ia, ib] = intersect(untreated_100{i,2}(1:imend_100(i)), tme_100);
    ftrace_100(i, ib)=avg_100{i}(1, ia);
end

for i=1:height(untreated_20)
    [~, ia, ib] = intersect(untreated_20{i,2}(1:imend_20(i)), tme_20);
    ftrace_20(i, ib)=avg_20{i}(1, ia);
end

%now calculate the standard error of the mean
smplesz_100 =  sum(~isnan(ftrace_100), 1);
ste_100 = std(ftrace_100, 0, 1, 'omitnan')./sqrt(smplesz_100);

smplesz_20 =  sum(~isnan(ftrace_20), 1);
ste_20 = std(ftrace_20, 0, 1, 'omitnan')./sqrt(smplesz_20);

%% plot the mean and ste of the mean for the corrected untreated LB trace
figure, hold on
ciplot(mean(ftrace_100, 1, 'omitnan')-ste_100, mean(ftrace_100, 1, 'omitnan')+ste_100, tme_100, colorcode2{1}, transparency)
plot(tme_100, mean(ftrace_100, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(ftrace_20, 1, 'omitnan')-ste_20, mean(ftrace_20, 1, 'omitnan')+ste_20, tme_20, colorcode2{7}, transparency)
plot(tme_20, mean(ftrace_20, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)

ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

legend('', '100% intensity', '', '20% intensity')
saveas(gcf, 'fullTrace.png');
saveas(gcf, 'fullTrace.fig');

%% fit the corrected untreated average
[tau_100, yhat_100]=tauCalc(tme_100', mean(ftrace_100, 1, 'omitnan'), 1);
[tau_20, yhat_20]=tauCalc(tme_20', mean(ftrace_20, 1, 'omitnan'), 1);

%find when fluor goes to zero
xpred = 0:420;
modelfun = @(tau, x)exp(-x./tau);
ypred_100 = modelfun(tau_100, xpred);
ypred_20 = modelfun(tau_20, xpred);

%% plot to see the fit
figure, hold on
plot(tme_100, mean(ftrace_100, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
plot(xpred, ypred_100, '--m', 'LineWidth', 1)

plot(tme_20, mean(ftrace_20, 1, 'omitnan'), 'Color', colorcode{7}, 'LineWidth', 1)
plot(xpred, ypred_20, '--k', 'LineWidth', 1)

ylim([0 1.1])
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
legend('100% intensity', '100% intensity fit', '20% intensity', '20% intensity fit')
saveas(gcf, 'tauFit.png')
saveas(gcf, 'tauFit.fig')

%% Functions
%to aggregate data and normalize
function [intensity, adjintensity, normintensity, lCell, time, tme, imstart]=dataNormalize(datadir, imstart)
        
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
      
    %remove cells without a value for imstart
%     idx=find(~isnan(intensity(:, imstart)));
%     intensity=intensity(idx,:);
%     lCell=lCell(idx,:);    
       
    %subtract final fluor. value from the trace
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<0)=NaN;      

%     %interpolate the fluor values during detergent perfusion
%     for n=1:height(adjintensity)
%         idx=find(adjintensity(n,:)>adjintensity(n, imstart));
%             if ~isempty(idx)
%                 [nrow, ncol]=size(adjintensity);
%                 x=setdiff(1:ncol, idx); %x=time
%                 v=adjintensity(n, x); %v=intensity
%                 vq=interp1(x, v, idx); %vq=interpolation at query time pts
%                 adjintensity(n,idx)=vq;
%             end
%     end
        
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:end)./adjintensity(:, imstart);
    [nrow, ncol]=size(normintensity);
    
    %interpolate the fluor values during detergent perfusion (& other noisy
    %spike in fluor.)
    for n=1:height(adjintensity)

        dv=diff(normintensity(n,:), 1, 2);
        idx1=find(dv>0)+1;
        idx2=find(normintensity(n,:)>1);
        idx=union(idx1, idx2);
        
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

%to aggregate control data and normalize
function [intensity, bgintensity, adjintensity, normintensity, lCell, time, tme, imstart]=controlNormalize(datadir)

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

%to find the time point at which the corrected trace plateaus and calculate
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