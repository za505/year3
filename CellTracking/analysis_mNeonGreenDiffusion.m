%Author: Zarina Akbary
%Date: 02/12/2022
%Purpose: To combine and analyze decayMeasure (dm.mat) data from mNeonGreen diffusion
%experiments. 

clear, close all

%Inputs
%dirsave = directory where all the dm.mat files are located and where plots
%are saved

%Outputs
%intensity = cell x time matrix of raw intensity values; concatenated
%icell_intensity variable from dm.mat files

%time = 1 x time matrix of time (minutes) values that correspond to
%intensity traces

%bgintensity = m x time matrix of background intensity values; concatenated
%bg_intensity variable from dm.mat files

%tme = 1 x n matrix of time (minutes) values that correspond to
%the pre- or post-lysis time point (imstart) and imend

%adjintensity = cell x tme maxtrix of intensity values minus bg intensity
%values from imstart to imend

%normintensity = cell x tme matrix of adjintensity values normalized to the
%initial value

%% User Input 

%color codes must be in RGB [0-1] format to be used in ciplot
colorcode={[204 0 0], [204 102 0], [204 204 0], [102 204 0], [0 204 204], [0 0 204], [102 0 204], [204 0 204], [204 0 102], [255 102 102], [255 178 102], [102 255 102], [102 255 255], [102 178 255], [178 102 255], [255 102 255],[255 102 178]};
colorcode2={[255 51 51], [255 153 51], [255 255 51], [153 255 51], [51 255 255], [51 51 255], [153 51 255], [255 51 255], [255 51 153], [255 204 204], [255 229 204], [204 255 204], [204 255 255], [204 229 255], [229 204 255], [255 204 255],[255 204 229]};

colorcode=cellfun(@(x)(x./255), colorcode, 'UniformOutput', false);
colorcode2=cellfun(@(x)(x./255), colorcode2, 'UniformOutput', false);

transparency = 0.3; %this is the alpha argument for the ciplot function

lengthCheck=1;

%location of the dm.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/02102022_reanalysis';
cd([dirsave '/MatFiles'])

LB1a=dir(['10232021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 1, full priming
LB1b=dir(['10262021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 2, full priming
LB20a=dir(['02092022_Exp1' '*dm.mat']); %LB, fr = 20 min, full priming

PBS2a=dir(['11192021_Exp2' '*dm.mat']); %PBS 2 min, rep 1
PBS2b=dir(['12082021_Exp3' '*dm.mat']); %PBS 2 min, rep 2
PBS20a=dir(['11192021_Exp1' '*dm.mat']); %PBS 20 min, rep 1
PBS20b=dir(['11302021_Exp1' '*dm.mat']); %PBS 20 min, rep 2
PBS60a=dir(['10232021_Exp2' '*dm.mat']); %PBS 60 min, rep 1
PBS60b=dir(['10262021_Exp2' '*dm.mat']); %PBS 60 min, rep 2
PBS120a=dir(['01142022_Exp1' '*dm.mat']); %PBS 120 min

LBMga=dir(['01172022_Exp1' '*dm.mat']); %LB, 20 mM Mg2+
LBEa=dir(['01172022_Exp2' '*dm.mat']); %LB, 10 mM EDTA
LBtuna=dir(['01242022_Exp1' '*dm.mat']); %LB, 0.5 ug/mL tunicamycin
LBvana=dir(['01262022_Exp1' '*dm.mat']); %LB, 1 ug/mL vancomycin

LB1s=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
LB2s=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
LB3s=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s

%% calculate normalized fluorescence traces
[intensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lcell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize(LB1a, 6); 
[intensity_LB1b, adjintensity_LB1b, normintensity_LB1b, lcell_LB1b, time_LB1b, tme_LB1b, imstart_LB1b]=dataNormalize(LB1b, 5); 
[intensity_LB20a, adjintensity_LB20a, normintensity_LB20a, lcell_LB20a, time_LB20a, tme_LB20a, imstart_LB20a]=dataNormalize(LB20a, 6); 

[intensity_PBS2a, adjintensity_PBS2a, normintensity_PBS2a, lcell_PBS2a, time_PBS2a, tme_PBS2a, imstart_PBS2a]=dataNormalize(PBS2a, 10); 
[intensity_PBS2b, adjintensity_PBS2b, normintensity_PBS2b, lcell_PBS2b, time_PBS2b, tme_PBS2b, imstart_PBS2b]=dataNormalize(PBS2b, 10); 
[intensity_PBS20a, adjintensity_PBS20a, normintensity_PBS20a, lcell_PBS20a, time_PBS20a, tme_PBS20a, imstart_PBS20a]=dataNormalize(PBS20a, 4); 
[intensity_PBS20b, adjintensity_PBS20b, normintensity_PBS20b, lcell_PBS20b, time_PBS20b, tme_PBS20b, imstart_PBS20b]=dataNormalize(PBS20b, 4); 
[intensity_PBS60a, adjintensity_PBS60a, normintensity_PBS60a, lcell_PBS60a, time_PBS60a, tme_PBS60a, imstart_PBS60a]=dataNormalize(PBS60a, 3); 
[intensity_PBS60b, adjintensity_PBS60b, normintensity_PBS60b, lcell_PBS60b, time_PBS60b, tme_PBS60b, imstart_PBS60b]=dataNormalize(PBS60b, 5); 
[intensity_PBS120a, adjintensity_PBS120a, normintensity_PBS120a, lcell_PBS120a, time_PBS120a, tme_PBS120a, imstart_PBS120a]=dataNormalize(PBS120a, 9); 

[intensity_LBMga, adjintensity_LBMga, normintensity_LBMga, lcell_LBMga, time_LBMga, tme_LBMga, imstart_LBMga]=dataNormalize(LBMga, 6); 
[intensity_LBEa, adjintensity_LBEa, normintensity_LBEa, lcell_LBEa, time_LBEa, tme_LBEa, imstart_LBEa]=dataNormalize(LBEa, 6); 
[intensity_LBtuna, adjintensity_LBtuna, normintensity_LBtuna, lcell_LBtuna, time_LBtuna, tme_LBtuna, imstart_LBtuna]=dataNormalize(LBtuna, 3); 
[intensity_LBvana, adjintensity_LBvana, normintensity_LBvana, lcell_LBvana, time_LBvana, tme_LBvana, imstart_LBvana]=dataNormalize(LBvana, 4); 

[intensity_LB1s, bgintensity_LB1s, adjintensity_LB1s, normintensity_LB1s, lcell_LB1s, time_LB1s, tme_LB1s, imstart_LB1s]=controlNormalize(LB1s); 
[intensity_LB2s, bgintensity_LB2s, adjintensity_LB2s, normintensity_LB2s, lcell_LB2s, time_LB2s, tme_LB2s, imstart_LB2s]=controlNormalize(LB2s); 
[intensity_LB3s, bgintensity_LB3s, adjintensity_LB3s, normintensity_LB3s, lcell_LB3s, time_LB3s, tme_LB3s, imstart_LB3s]=controlNormalize(LB3s); 

%% plot length traces, as a sanity check for imstart
if lengthCheck==1
    figure, hold on
    ciplot(mean(lcell_LB1a, 1, 'omitnan')-std(lcell_LB1a, 0, 1, 'omitnan'), mean(lcell_LB1a, 1, 'omitnan')+std(lcell_LB1a, 0, 1, 'omitnan'), time_LB1a, colorcode2{1}, transparency)
    plot(time_LB1a, mean(lcell_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_LB1b, 1, 'omitnan')-std(lcell_LB1b, 0, 1, 'omitnan'), mean(lcell_LB1b, 1, 'omitnan')+std(lcell_LB1b, 0, 1, 'omitnan'), time_LB1b, colorcode2{5}, transparency)
    plot(time_LB1b, mean(lcell_LB1b, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_LB20a, 1, 'omitnan')-std(lcell_LB20a, 0, 1, 'omitnan'), mean(lcell_LB20a, 1, 'omitnan')+std(lcell_LB20a, 0, 1, 'omitnan'), time_LB20a, colorcode2{4}, transparency)
    plot(time_LB20a, mean(lcell_LB20a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_PBS2a, 1, 'omitnan')-std(lcell_PBS2a, 0, 1, 'omitnan'), mean(lcell_PBS2a, 1, 'omitnan')+std(lcell_PBS2a, 0, 1, 'omitnan'), time_PBS2a, colorcode2{6}, transparency)
    plot(time_PBS2a, mean(lcell_PBS2a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_PBS2b, 1, 'omitnan')-std(lcell_PBS2b, 0, 1, 'omitnan'), mean(lcell_PBS2b, 1, 'omitnan')+std(lcell_PBS2b, 0, 1, 'omitnan'), time_PBS2b, colorcode2{6}, transparency)
    plot(time_PBS2b, mean(lcell_PBS2b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_PBS20a, 1, 'omitnan')-std(lcell_PBS20a, 0, 1, 'omitnan'), mean(lcell_PBS20a, 1, 'omitnan')+std(lcell_PBS20a, 0, 1, 'omitnan'), time_PBS20a, colorcode2{6}, transparency)
    plot(time_PBS20a, mean(lcell_PBS20a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_PBS20b, 1, 'omitnan')-std(lcell_PBS20b, 0, 1, 'omitnan'), mean(lcell_PBS20b, 1, 'omitnan')+std(lcell_PBS20b, 0, 1, 'omitnan'), time_PBS20b, colorcode2{6}, transparency)
    plot(time_PBS20b, mean(lcell_PBS20b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_PBS60a, 1, 'omitnan')-std(lcell_PBS60a, 0, 1, 'omitnan'), mean(lcell_PBS60a, 1, 'omitnan')+std(lcell_PBS60a, 0, 1, 'omitnan'), time_PBS60a, colorcode2{6}, transparency)
    plot(time_PBS60a, mean(lcell_PBS60a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_PBS60b, 1, 'omitnan')-std(lcell_PBS60b, 0, 1, 'omitnan'), mean(lcell_PBS60b, 1, 'omitnan')+std(lcell_PBS60b, 0, 1, 'omitnan'), time_PBS60b, colorcode2{6}, transparency)
    plot(time_PBS60b, mean(lcell_PBS60b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_PBS120a, 1, 'omitnan')-std(lcell_PBS120a, 0, 1, 'omitnan'), mean(lcell_PBS120a, 1, 'omitnan')+std(lcell_PBS120a, 0, 1, 'omitnan'), time_PBS120a, colorcode2{6}, transparency)
    plot(time_PBS120a, mean(lcell_PBS120a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_LBMga, 1, 'omitnan')-std(lcell_LBMga, 0, 1, 'omitnan'), mean(lcell_LBMga, 1, 'omitnan')+std(lcell_LBMga, 0, 1, 'omitnan'), time_LBMga, colorcode2{6}, transparency)
    plot(time_LBMga, mean(lcell_LBMga, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_LBEa, 1, 'omitnan')-std(lcell_LBEa, 0, 1, 'omitnan'), mean(lcell_LBEa, 1, 'omitnan')+std(lcell_LBEa, 0, 1, 'omitnan'), time_LBEa, colorcode2{3}, transparency)
    plot(time_LBEa, mean(lcell_LBEa, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_LBtuna, 1, 'omitnan')-std(lcell_LBtuna, 0, 1, 'omitnan'), mean(lcell_LBtuna, 1, 'omitnan')+std(lcell_LBtuna, 0, 1, 'omitnan'), time_LBtuna, colorcode2{6}, transparency)
    plot(time_LBtuna, mean(lcell_LBtuna, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    figure, hold on
    ciplot(mean(lcell_LBvana, 1, 'omitnan')-std(lcell_LBvana, 0, 1, 'omitnan'), mean(lcell_LBvana, 1, 'omitnan')+std(lcell_LBvana, 0, 1, 'omitnan'), time_LBvana, colorcode2{6}, transparency)
    plot(time_LBvana, mean(lcell_LBvana, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)

    xlabel('Time (minutes)')
    ylabel('Length')

end
%% plot the raw fluorescence intensity
figure, hold on
plot(time_LB1a, intensity_LB1a, 'Color', colorcode{1});
plot(time_LB1a, bgintensity_LB1a, '--k');

figure, hold on
plot(time_LB1b, intensity_LB1b, 'Color', colorcode{5});
plot(time_LB1b, bgintensity_LB1b, '--k');

figure, hold on
plot(time_LB20a, intensity_LB20a, 'Color', colorcode{3});
plot(time_LB20a, bgintensity_LB20a, '--k');
 
% figure, hold on
% plot(time_LBEa, intensity_LBc, 'Color', colorcode{3});
% plot(time_LBEa, bgintensity_LBEa, '--k');

% figure, hold on
% plot(time_LBEa, intensity_LBEa, 'Color', colorcode{3});
% plot(time_LB20a, intensity_LB20a, 'Color', colorcode{4});
% plot(time_PBS2a, intensity_PBS2a, 'Color', colorcode{6});
% 
% figure, hold on
% plot(time_LB1a, intensity_LB1a, 'Color', colorcode{1});
% plot(time_LB1b, intensity_LB1b, 'Color', colorcode{5});
% plot(time_LB20a, intensity_LB20a, 'Color', colorcode{4});
% plot(time_PBS2a, intensity_PBS2a, 'Color', colorcode{6});

figure, hold on
plot(time_LB1s, intensity_LB1s, '-r');
plot(time_LB1s, bgintensity_LB1s, '--r')
plot(time_LB2s, intensity_LB2s, '-b');
plot(time_LB2s, bgintensity_LB2s, '--b')
plot(time_LB3s, intensity_LB3s, '-g'); 
plot(time_LB3s, bgintensity_LB3s, '--g')

% figure, hold on
% ciplot(mean(intensity_LB1a, 1, 'omitnan')-std(intensity_LB1a, 0, 1, 'omitnan'), mean(intensity_LB1a, 1, 'omitnan')+std(intensity_LB1a, 0, 1, 'omitnan'), time_LB1a, colorcode2{1}, transparency)
% plot(time_LB1a, mean(intensity_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% ciplot(mean(intensity_LB1b, 1, 'omitnan')-std(intensity_LB1b, 0, 1, 'omitnan'), mean(intensity_LB1b, 1, 'omitnan')+std(intensity_LB1b, 0, 1, 'omitnan'), time_LB1b, colorcode2{5}, transparency)
% plot(time_LB1b, mean(intensity_LB1b, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% 
% ciplot(mean(intensity_LBEa, 1, 'omitnan')-std(intensity_LBEa, 0, 1, 'omitnan'), mean(intensity_LBEa, 1, 'omitnan')+std(intensity_LBEa, 0, 1, 'omitnan'), time_LBEa, colorcode2{3}, transparency)
% plot(time_LBEa, mean(intensity_LBEa, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% 
% ciplot(mean(intensity_LB20a, 1, 'omitnan')-std(intensity_LB20a, 0, 1, 'omitnan'), mean(intensity_LB20a, 1, 'omitnan')+std(intensity_LB20a, 0, 1, 'omitnan'), time_LB20a, colorcode2{4}, transparency)
% plot(time_LB20a, mean(intensity_LB20a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)
% 
% ciplot(mean(intensity_PBS2a, 1, 'omitnan')-std(intensity_PBS2a, 0, 1, 'omitnan'), mean(intensity_PBS2a, 1, 'omitnan')+std(intensity_PBS2a, 0, 1, 'omitnan'), time_PBS2a, colorcode2{6}, transparency)
% plot(time_PBS2a, mean(intensity_PBS2a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
% xlabel('Time (minutes)')
% ylabel('Intensity')

% figure, hold on
% ciplot(mean(intensity_LB1s, 1, 'omitnan')-std(intensity_LB1s, 0, 1, 'omitnan'), mean(intensity_LB1s, 1, 'omitnan')+std(intensity_LB1s, 0, 1, 'omitnan'), time_LB1s, colorcode2{1}, transparency)
% plot(time_LB1s, mean(intensity_LB1s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% ciplot(mean(intensity_LB2s, 1, 'omitnan')-std(intensity_LB2s, 0, 1, 'omitnan'), mean(intensity_LB2s, 1, 'omitnan')+std(intensity_LB2s, 0, 1, 'omitnan'), time_LB2s, colorcode2{5}, transparency)
% plot(time_LB2s, mean(intensity_LB2s, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% 
% ciplot(mean(intensity_LB3s, 1, 'omitnan')-std(intensity_LB3s, 0, 1, 'omitnan'), mean(intensity_LB3s, 1, 'omitnan')+std(intensity_LB3s, 0, 1, 'omitnan'), time_LB3s, colorcode2{3}, transparency)
% plot(time_LB3s, mean(intensity_LB3s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% xlabel('Time (minutes)')
% ylabel('Intensity')

%Conclusion: LB3s was cut short (before it hit the plateau) and LB1a
%started with the highest initial value. The rest all start at similar
%intiial values but end at different final values. I wonder if the LB1s and
%LB2s movies bleach the background faster than the flow of media through
%the chip?

%% plot the adjusted intensity
%figure, hold on
figure, plot(time_LB1a, adjintensity_LB1a, 'Color', colorcode{1});
figure, plot(time_LB1b, adjintensity_LB1b, 'Color', colorcode{5});
figure, plot(time_LB20a, adjintensity_LB20a, 'Color', colorcode{3});

% figure, hold on
% plot(time_LBEa, adjintensity_LBEa, 'Color', colorcode{3});
% plot(time_LB20a, adjintensity_LB20a, 'Color', colorcode{4});
% plot(time_PBS2a, adjintensity_PBS2a, 'Color', colorcode{6});
% 
% figure, hold on
% plot(time_LB1a, adjintensity_LB1a, 'Color', colorcode{1});
% plot(time_LB1b, adjintensity_LB1b, 'Color', colorcode{5});
% plot(time_LB20a, adjintensity_LB20a, 'Color', colorcode{4});
% plot(time_PBS2a, adjintensity_PBS2a, 'Color', colorcode{6});

figure, hold on
plot(time_LB1s, adjintensity_LB1s, '-r');
plot(time_LB2s, adjintensity_LB2s, '-b');
plot(time_LB3s, adjintensity_LB3s, '-g'); 

% figure, hold on
% ciplot(mean(adjintensity_LB1a, 1, 'omitnan')-std(adjintensity_LB1a, 0, 1, 'omitnan'), mean(adjintensity_LB1a, 1, 'omitnan')+std(adjintensity_LB1a, 0, 1, 'omitnan'), time_LB1a, colorcode2{1}, transparency)
% plot(time_LB1a, mean(adjintensity_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LB1b, 1, 'omitnan')-std(adjintensity_LB1b, 0, 1, 'omitnan'), mean(adjintensity_LB1b, 1, 'omitnan')+std(adjintensity_LB1b, 0, 1, 'omitnan'), time_LB1b, colorcode2{5}, transparency)
% plot(time_LB1b, mean(adjintensity_LB1b, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LBEa, 1, 'omitnan')-std(adjintensity_LBEa, 0, 1, 'omitnan'), mean(adjintensity_LBEa, 1, 'omitnan')+std(adjintensity_LBEa, 0, 1, 'omitnan'), time_LBEa, colorcode2{3}, transparency)
% plot(time_LBEa, mean(adjintensity_LBEa, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LB20a, 1, 'omitnan')-std(adjintensity_LB20a, 0, 1, 'omitnan'), mean(adjintensity_LB20a, 1, 'omitnan')+std(adjintensity_LB20a, 0, 1, 'omitnan'), time_LB20a, colorcode2{4}, transparency)
% plot(time_LB20a, mean(adjintensity_LB20a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_PBS2a, 1, 'omitnan')-std(adjintensity_PBS2a, 0, 1, 'omitnan'), mean(adjintensity_PBS2a, 1, 'omitnan')+std(adjintensity_PBS2a, 0, 1, 'omitnan'), time_PBS2a, colorcode2{6}, transparency)
% plot(time_PBS2a, mean(adjintensity_PBS2a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
% xlabel('Time (minutes)')
% ylabel('Intensity')

% figure, hold on
% ciplot(mean(adjintensity_LB1s, 1, 'omitnan')-std(adjintensity_LB1s, 0, 1, 'omitnan'), mean(adjintensity_LB1s, 1, 'omitnan')+std(adjintensity_LB1s, 0, 1, 'omitnan'), time_LB1s, colorcode2{1}, transparency)
% plot(time_LB1s, mean(adjintensity_LB1s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LB2s, 1, 'omitnan')-std(adjintensity_LB2s, 0, 1, 'omitnan'), mean(adjintensity_LB2s, 1, 'omitnan')+std(adjintensity_LB2s, 0, 1, 'omitnan'), time_LB2s, colorcode2{5}, transparency)
% plot(time_LB2s, mean(adjintensity_LB2s, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LB3s, 1, 'omitnan')-std(adjintensity_LB3s, 0, 1, 'omitnan'), mean(adjintensity_LB3s, 1, 'omitnan')+std(adjintensity_LB3s, 0, 1, 'omitnan'), time_LB3s, colorcode2{3}, transparency)
% plot(time_LB3s, mean(adjintensity_LB3s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% xlabel('Time (minutes)')
% ylabel('Intensity')

%% plot normalized fluorescence traces
% figure, hold on
% plot(tme_LB1a, normintensity_LB1a, 'Color', colorcode{1});
% plot(tme_LB1b, normintensity_LB1b, 'Color', colorcode{5});
% plot(tme_LBEa, normintensity_LBEa, 'Color', colorcode{3});
% 
% figure, hold on
% plot(tme_LBEa, normintensity_LBEa, 'Color', colorcode{3});
% plot(tme_LB20a, normintensity_LB20a, 'Color', colorcode{4});
% plot(tme_PBS2a, normintensity_PBS2a, 'Color', colorcode{6});
% 
% figure, hold on
% plot(tme_LB1a, normintensity_LB1a, 'Color', colorcode{1});
% plot(tme_LB1b, normintensity_LB1b, 'Color', colorcode{5});
% plot(tme_LB20a, normintensity_LB20a, 'Color', colorcode{4});
% plot(tme_PBS2a, normintensity_PBS2a, 'Color', colorcode{6});

figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r');
plot(tme_LB2s, normintensity_LB2s, '-b');
plot(tme_LB3s, normintensity_LB3s, '-g'); 

figure, hold on
ciplot(mean(normintensity_LB1a, 1, 'omitnan')-std(normintensity_LB1a, 0, 1, 'omitnan'), mean(normintensity_LB1a, 1, 'omitnan')+std(normintensity_LB1a, 0, 1, 'omitnan'), tme_LB1a, colorcode2{1}, transparency)
plot(tme_LB1a, mean(normintensity_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(normintensity_LB1b, 1, 'omitnan')-std(normintensity_LB1b, 0, 1, 'omitnan'), mean(normintensity_LB1b, 1, 'omitnan')+std(normintensity_LB1b, 0, 1, 'omitnan'), tme_LB1b, colorcode2{5}, transparency)
plot(tme_LB1b, mean(normintensity_LB1b, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

% ciplot(mean(normintensity_LBEa, 1, 'omitnan')-std(normintensity_LBEa, 0, 1, 'omitnan'), mean(normintensity_LBEa, 1, 'omitnan')+std(normintensity_LBEa, 0, 1, 'omitnan'), tme_LBEa, colorcode2{3}, transparency)
% plot(tme_LBEa, mean(normintensity_LBEa, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

% ciplot(mean(normintensity_LB20a, 1, 'omitnan')-std(normintensity_LB20a, 0, 1, 'omitnan'), mean(normintensity_LB20a, 1, 'omitnan')+std(normintensity_LB20a, 0, 1, 'omitnan'), tme_LB20a, colorcode2{4}, transparency)
% plot(tme_LB20a, mean(normintensity_LB20a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

ciplot(mean(normintensity_PBS2a, 1, 'omitnan')-std(normintensity_PBS2a, 0, 1, 'omitnan'), mean(normintensity_PBS2a, 1, 'omitnan')+std(normintensity_PBS2a, 0, 1, 'omitnan'), tme_PBS2a, colorcode2{6}, transparency)
plot(tme_PBS2a, mean(normintensity_PBS2a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
legend({'LB 1 min, rep 1, std', 'LB 1 min, rep 1, avg', 'LB 1 min, rep 2, std', 'LB 1 min, rep 2, avg', 'LB 20 min, std', 'LB 20 min, avg'})
xlabel('Time (minutes)')
ylabel('Normalized Intensity')

figure, hold on
ciplot(mean(normintensity_LB1s, 1, 'omitnan')-std(normintensity_LB1s, 0, 1, 'omitnan'), mean(normintensity_LB1s, 1, 'omitnan')+std(normintensity_LB1s, 0, 1, 'omitnan'), tme_LB1s, colorcode2{1}, transparency)
plot(tme_LB1s, mean(normintensity_LB1s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(normintensity_LB2s, 1, 'omitnan')-std(normintensity_LB2s, 0, 1, 'omitnan'), mean(normintensity_LB2s, 1, 'omitnan')+std(normintensity_LB2s, 0, 1, 'omitnan'), tme_LB2s, colorcode2{5}, transparency)
plot(tme_LB2s, mean(normintensity_LB2s, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(normintensity_LB3s, 1, 'omitnan')-std(normintensity_LB3s, 0, 1, 'omitnan'), mean(normintensity_LB3s, 1, 'omitnan')+std(normintensity_LB3s, 0, 1, 'omitnan'), tme_LB3s, colorcode2{3}, transparency)
plot(tme_LB3s, mean(normintensity_LB3s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Intensity')

xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

%% calculate alpha
[tau1, yhat_LB1s]=tauCalc(tme_LB1s, normintensity_LB1s);
[tau2, yhat_LB2s]=tauCalc(tme_LB2s, normintensity_LB2s);
[tau3, yhat_LB3s]=tauCalc(tme_LB3s, normintensity_LB3s);

figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r')
plot(tme_LB1s, yhat_LB1s, '--k')
%saveas(gcf, 'tauFit_LB1s.png')

figure, hold on
plot(tme_LB2s, normintensity_LB2s, '-b')
plot(tme_LB2s, yhat_LB2s, '--k')
%saveas(gcf, 'tauFit_LB2s.png')

figure, hold on
plot(tme_LB3s, normintensity_LB3s, '-g')
plot(tme_LB3s, yhat_LB3s, '--k')
%saveas(gcf, 'tauFit_LB3s.png')

tau_means = [mean(tau1, 'omitnan'), mean(tau2, 'omitnan'), mean(tau3, 'omitnan')];
tau_std = [std(tau1, 0, 'omitnan'), std(tau2, 0, 'omitnan'), std(tau3, 0, 'omitnan')];

dt1=tme_LB1s(end)-tme_LB1s(end-1);
dt2=tme_LB2s(end)-tme_LB2s(end-1);
dt3=tme_LB3s(end)-tme_LB3s(end-1);

linearCoef1 = polyfit([repelem(dt1, length(tau1)), repelem(dt2, length(tau2)), repelem(dt3, length(tau3))],[tau1', tau2', tau3'],1);
linearFit1= polyval(linearCoef1,[0 dt1 dt2 dt3]);

figure, hold on
scatter(repelem(dt1, length(tau1)), tau1, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(dt2, length(tau2)), tau2, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
scatter(repelem(dt3, length(tau3)), tau3, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')

scatter([dt1, dt2, dt3], tau_means, 'MarkerFaceColor', 'black')
errorbar(dt1, tau_means(1), tau_std(1), 'Color', 'black')
errorbar(dt2, tau_means(2), tau_std(2), 'Color', 'black')
errorbar(dt3, tau_means(3), tau_std(3), 'Color', 'black')

plot([0 dt1 dt2 dt3], linearFit1, '--b')
xlim([0, 0.06])
ylabel('\tau (min^{-1})')
xticks([0 dt1 dt2 dt3])
xticklabels({'0', '1.2 s', '2 s', '3 s'})
title('Tau vs Frame Rate')

%% correct for photobleaching
% alpha=38.898;
% intercept=-0.0765;

% alpha=28.9210;
% intercept=0.2482;

% alpha=39.1226;
% intercept=-0.0827;

%02/10/2022
alpha=32.2114;
intercept=0.1614;

%02/11/2022
% alpha=53.1363;
% intercept=-0.2943;

[Cnew_LB1a, dCB_LB1a, dCT_LB1a, dCP_LB1a, CblExp_LB1a, unbFrac_LB1a]=photoCorrect(tme_LB1a, normintensity_LB1a, alpha, intercept);
[Cnew_LB1b, dCB_LB1b, dCT_LB1b, dCP_LB1b, CblExp_LB1b, unbFrac_LB1b]=photoCorrect(tme_LB1b, normintensity_LB1b, alpha, intercept);
[Cnew_LB20a, dCB_LB20a, dCT_LB20a, dCP_LB20a, CblExp_LB20a, unbFrac_LB20a]=photoCorrect(tme_LB20a, normintensity_LB20a, alpha, intercept);

[Cnew_LB1s, dCB_LB1s, dCT_LB1s, dCP_LB1s, CblExp_LB1s, unbFrac_LB1s]=photoCorrect(tme_LB1s, normintensity_LB1s, alpha, intercept);
[Cnew_LB2s, dCB_LB2s, dCT_LB2s, dCP_LB2s, CblExp_LB2s, unbFrac_LB2s]=photoCorrect(tme_LB2s, normintensity_LB2s, alpha, intercept);
[Cnew_LB3s, dCB_LB3s, dCT_LB3s, dCP_LB3s, CblExp_LB3s, unbFrac_LB3s]=photoCorrect(tme_LB3s, normintensity_LB3s, alpha, intercept);

%% plot corrected traces
figure, plot(tme_LB1a, Cnew_LB1a, '-r');
figure, plot(tme_LB1b, Cnew_LB1b, '-b');
figure, plot(tme_LB20a, Cnew_LB20a, '-g'); 

figure, hold on
%ciplot(mean(Cnew_LB, 1, 'omitnan')-std(Cnew_LB, 0, 1, 'omitnan'), mean(Cnew_LB, 1, 'omitnan')+std(Cnew_LB, 0, 1, 'omitnan'), tme_LB, colorcode2{1}, transparency)
plot(tme_LB1a, mean(Cnew_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(Cnew_LBEa, 1, 'omitnan')-std(Cnew_LBEa, 0, 1, 'omitnan'), mean(Cnew_LBEa, 1, 'omitnan')+std(Cnew_LBEa, 0, 1, 'omitnan'), tme_LBEa, colorcode2{5}, transparency)
plot(tme_LB1b, mean(Cnew_LB1b, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(Cnew_LB20, 1, 'omitnan')-std(Cnew_LB20, 0, 1, 'omitnan'), mean(Cnew_LB20, 1, 'omitnan')+std(Cnew_LB20, 0, 1, 'omitnan'), tme_LB20, colorcode2{3}, transparency)
plot(tme_LB20a, mean(Cnew_LB20a, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
legend({'LB, 1 min, rep 1, avg', 'LB, 1 min, rep 2, avg', 'LB, 20 min, avg'})
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

figure, hold on
plot(tme_LB1s, Cnew_LB1s, '-r');
plot(tme_LB2s, Cnew_LB2s, '-b');
plot(tme_LB3s, Cnew_LB3s, '-g'); 
ylim([0 1.1])

figure, hold on
plot(tme_LB1s, mean(Cnew_LB1s, 1, 'omitnan'), '-r');
plot(tme_LB2s,mean(Cnew_LB2s, 1, 'omitnan'), '-b');
plot(tme_LB3s, mean(Cnew_LB3s, 1, 'omitnan'), '-g'); 

%% sanity check
time=tme_LB1s(1:end-1);
dCB_LB1s=dCB_LB1s.*-1;

figure, hold on
%ciplot(mean(dCT_LB1s, 1, 'omitnan')-std(dCT_LB1s, 0, 1, 'omitnan'), mean(dCT_LB1s, 1, 'omitnan')+std(dCT_LB1s, 0, 1, 'omitnan'), time, colorcode2{1}, transparency)
plot(time, mean(dCT_LB1s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(dCB_LB1s, 1, 'omitnan')-std(dCB_LB1s, 0, 1, 'omitnan'), mean(dCB_LB1s, 1, 'omitnan')+std(dCB_LB1s, 0, 1, 'omitnan'), time, colorcode2{5}, transparency)
plot(time, mean(dCB_LB1s, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(dCP_LB1s, 1, 'omitnan')-std(dCP_LB1s, 0, 1, 'omitnan'), mean(dCP_LB1s, 1, 'omitnan')+std(dCP_LB1s, 0, 1, 'omitnan'), time, colorcode2{3}, transparency)
plot(time, mean(dCP_LB1s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
legend({'dCT', 'dCB', 'dCP'})
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

time=tme_LB1s(1:end-1);
dCB_LB1s=dCB_LB1s.*-1;

figure, hold on
%ciplot(mean(dCT_LB1s, 1, 'omitnan')-std(dCT_LB1s, 0, 1, 'omitnan'), mean(dCT_LB1s, 1, 'omitnan')+std(dCT_LB1s, 0, 1, 'omitnan'), time, colorcode2{1}, transparency)
plot(time, mean(dCT_LB1s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(dCB_LB1s, 1, 'omitnan')-std(dCB_LB1s, 0, 1, 'omitnan'), mean(dCB_LB1s, 1, 'omitnan')+std(dCB_LB1s, 0, 1, 'omitnan'), time, colorcode2{5}, transparency)
plot(time, mean(dCB_LB1s, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(dCP_LB1s, 1, 'omitnan')-std(dCP_LB1s, 0, 1, 'omitnan'), mean(dCP_LB1s, 1, 'omitnan')+std(dCP_LB1s, 0, 1, 'omitnan'), time, colorcode2{3}, transparency)
plot(time, mean(dCP_LB1s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
legend({'dCT', 'dCB', 'dCP'})
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

% sanity check for 2s frame rate
time=tme_LB2s(1:end-1);
dCB_LB2s=dCB_LB2s.*-1;

figure, hold on
%ciplot(mean(dCT_LB2s, 1, 'omitnan')-std(dCT_LB2s, 0, 1, 'omitnan'), mean(dCT_LB2s, 1, 'omitnan')+std(dCT_LB2s, 0, 1, 'omitnan'), time, colorcode2{1}, transparency)
plot(time, mean(dCT_LB2s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(dCB_LB2s, 1, 'omitnan')-std(dCB_LB2s, 0, 1, 'omitnan'), mean(dCB_LB2s, 1, 'omitnan')+std(dCB_LB2s, 0, 1, 'omitnan'), time, colorcode2{5}, transparency)
plot(time, mean(dCB_LB2s, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(dCP_LB2s, 1, 'omitnan')-std(dCP_LB2s, 0, 1, 'omitnan'), mean(dCP_LB2s, 1, 'omitnan')+std(dCP_LB2s, 0, 1, 'omitnan'), time, colorcode2{3}, transparency)
plot(time, mean(dCP_LB2s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
legend({'dCT', 'dCB', 'dCP'})
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

% last sanity check
time=tme_LB3s(1:end-1);
dCB_LB3s=dCB_LB3s.*-1;

figure, hold on
%ciplot(mean(dCT_LB3s, 1, 'omitnan')-std(dCT_LB3s, 0, 1, 'omitnan'), mean(dCT_LB3s, 1, 'omitnan')+std(dCT_LB3s, 0, 1, 'omitnan'), time, colorcode2{1}, transparency)
plot(time, mean(dCT_LB3s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(dCB_LB3s, 1, 'omitnan')-std(dCB_LB3s, 0, 1, 'omitnan'), mean(dCB_LB3s, 1, 'omitnan')+std(dCB_LB3s, 0, 1, 'omitnan'), time, colorcode2{5}, transparency)
plot(time, mean(dCB_LB3s, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

%ciplot(mean(dCP_LB3s, 1, 'omitnan')-std(dCP_LB3s, 0, 1, 'omitnan'), mean(dCP_LB3s, 1, 'omitnan')+std(dCP_LB3s, 0, 1, 'omitnan'), time, colorcode2{3}, transparency)
plot(time, mean(dCP_LB3s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
legend({'dCT', 'dCB', 'dCP'})
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')

%% Functions
%to aggregate data and normalize
function [intensity, adjintensity, normintensity, lCell, time, tme, imstart]=dataNormalize(datadir)
        
    %pre-allocate variables
    intensity=[];
    lCell=[];
        
    %go through the data for each position
    for i=1:length(datadir)

        %load decayMeasure .mat file
        cd(datadir(i).folder)
        load(datadir(i).name, 'icell_intensity', 'time', 'lcell')

        intensity=[intensity; icell_intensity];
        lCell=[lCell; lcell];

        if i==1
            tme=time; %pre-set new time vector
        end

    end
        
    %find the final pre-lysis frame by identifying the biggest
    %change in length
    dl=diff(lCell, 1);
    lvg=mean(dl, 1, 'omitnan');
    imstart=find(lvg<0, 1, 'first');
    %[~, imstart]=min(lvg);
      
    %remove cells without a value for imstart
    idx=find(~isnan(intensity(:, imstart)));
    intensity=intensity(idx,:);
    lCell=lCell(idx,:);    
        
    %subtract final fluor. value from the trace
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<0)=NaN;      

    %interpolate the fluor values during detergent perfusion
    for n=1:height(adjintensity)
        idx=find(adjintensity(n,:)>adjintensity(n, imstart));
            if ~isempty(idx)
                [nrow, ncol]=size(adjintensity);
                x=setdiff(1:ncol, idx); %x=time
                v=adjintensity(n, x); %v=intensity
                vq=interp1(x, v, idx); %vq=interpolation at query time pts
                adjintensity(n,idx)=vq;
            end
    end
        
    %adjust the time points in the fluor matrix and initialize to first
    %frame
    normintensity=adjintensity(:, imstart:end)./adjintensity(:, imstart);

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
        dl=diff(lCell, 1);
        lvg=mean(dl, 1, 'omitnan');
        imstart=find(lvg<0, 1, 'first')+2;
    else
        if dt(end)<1 %discrepancies in dt might be found in the meta data, this is the ad hoc fix
            imstart=find(dt==dt(end), 1, 'first')+2;
        else 
            imstart=find(dt==dt(end), 1, 'first');
        end
    end
  
    %remove cells without a value for imstart
    idx=find(~isnan(intensity(:, imstart)));
    intensity=intensity(idx,:);
    lCell=lCell(idx,:);    
        
    %subtract final fluor. value from the trace and set limit of
    %detection
    adjintensity=intensity-intensity(:, end); 
    adjintensity(adjintensity<200)=NaN; 

        
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
function [tau, yhat]=tauCalc(tme, normintensity)
    
    [nrow, ~]=size(normintensity);
    tau=nan(nrow);
    modelfun = @(tau, x)exp(-x./tau);
    for i=1:nrow
        tau(i, 1)=nlinfit(tme, normintensity(i,:), modelfun);
    end
    
    yhat=modelfun(tau, tme);
    
end

%to correct for photobleaching
function [Cnew, dCB, dCT, dCP, Cbl_exp, unb_frac]=photoCorrect(tme, normintensity, alpha, intercept)
        
        %pre-allocate variables
        %assume that the initial 'measured' fluorescence values and corrected
        %fluor. values will be equal. I prefer to pre-allocate with nan in case 
        %some values are missing in the raw data
        Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
        Cnew(:, 1)=normintensity(:);
             
        dCB=nan(height(normintensity), length(tme)-1); %change in fluor. due to photobleaching
        dCT=nan(height(normintensity), length(tme)-1); %total change in fluor.
        dCP=nan(height(normintensity), length(tme)-1); %this is the dCP, or loss attributable to permeability

        unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
        unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

        Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
        Cbl_exp(:, 1)=0;
        
        %calculate dt (the dt between frames may vary in a single run).
        %Note: this is also the frame rate
        dt=round(diff(tme));
        
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
