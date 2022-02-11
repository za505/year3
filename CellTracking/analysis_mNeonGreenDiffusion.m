%Author: Zarina Akbary
%Date: 02/11/2022
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

%location of the dm.mat files 
dirsave='/Users/zarina/Downloads/NYU/Year3_2022_Spring/02102022_reanalysis';
cd(dirsave)

LB1a=dir(['10232021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 1, full priming
LB1b=dir(['10262021_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 2, full priming
LB1c=dir(['02082022_Exp1' '*dm.mat']); %LB, fr = 1 min, rep 2, short priming
LB20a=dir(['02092022_Exp1' '*dm.mat']); %LB, fr = 20 min, full priming
LB20b=dir(['02042022_Exp1' '*dm.mat']); %LB, fr = 20 min, short priming

LB1s=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate = 1.2 s
LB2s=dir(['12082021_Exp1' '*dm.mat']); %LB, frame rate = 2 s
LB3s=dir(['12082021_Exp2' '*dm.mat']); %LB, frame rate = 3 s

%% calculate normalized fluorescence traces
[intensity_LB1a, bgintensity_LB1a, adjintensity_LB1a, normintensity_LB1a, lcell_LB1a, time_LB1a, tme_LB1a, imstart_LB1a]=dataNormalize1(LB1a, 1); 
[intensity_LB1b, bgintensity_LB1b, adjintensity_LB1b, normintensity_LB1b, lcell_LB1b, time_LB1b, tme_LB1b, imstart_LB1b]=dataNormalize1(LB1b, 1); 
[intensity_LB1c, bgintensity_LB1c, adjintensity_LB1c, normintensity_LB1c, lcell_LB1c, time_LB1c, tme_LB1c, imstart_LB1c]=dataNormalize1(LB1c, 1); 
[intensity_LB20a, bgintensity_LB20a, adjintensity_LB20a, normintensity_LB20a, lcell_LB20a, time_LB20a, tme_LB20a, imstart_LB20a]=dataNormalize1(LB20a, 1); 
[intensity_LB20b, bgintensity_LB20b, adjintensity_LB20b, normintensity_LB20b, lcell_LB20b, time_LB20b, tme_LB20b, imstart_LB20b]=dataNormalize1(LB20b, 1); 

[intensity_LB1s, bgintensity_LB1s, adjintensity_LB1s, normintensity_LB1s, lcell_LB1s, time_LB1s, tme_LB1s, imstart_LB1s]=dataNormalize1(LB1s, 2); 
[intensity_LB2s, bgintensity_LB2s, adjintensity_LB2s, normintensity_LB2s, lcell_LB2s, time_LB2s, tme_LB2s, imstart_LB2s]=dataNormalize1(LB2s, 2); 
[intensity_LB3s, bgintensity_LB3s, adjintensity_LB3s, normintensity_LB3s, lcell_LB3s, time_LB3s, tme_LB3s, imstart_LB3s]=dataNormalize1(LB3s, 2); 

%% plot length traces
figure, hold on
ciplot(mean(lcell_LB1a, 1, 'omitnan')-std(lcell_LB1a, 0, 1, 'omitnan'), mean(lcell_LB1a, 1, 'omitnan')+std(lcell_LB1a, 0, 1, 'omitnan'), time_LB1a, colorcode2{1}, transparency)
plot(time_LB1a, mean(lcell_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(lcell_LB1b, 1, 'omitnan')-std(lcell_LB1b, 0, 1, 'omitnan'), mean(lcell_LB1b, 1, 'omitnan')+std(lcell_LB1b, 0, 1, 'omitnan'), time_LB1b, colorcode2{5}, transparency)
plot(time_LB1b, mean(lcell_LB1b, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

% ciplot(mean(lcell_LB1c, 1, 'omitnan')-std(lcell_LB1c, 0, 1, 'omitnan'), mean(lcell_LB1c, 1, 'omitnan')+std(lcell_LB1c, 0, 1, 'omitnan'), time_LB1c, colorcode2{3}, transparency)
% plot(time_LB1c, mean(lcell_LB1c, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

ciplot(mean(lcell_LB20a, 1, 'omitnan')-std(lcell_LB20a, 0, 1, 'omitnan'), mean(lcell_LB20a, 1, 'omitnan')+std(lcell_LB20a, 0, 1, 'omitnan'), time_LB20a, colorcode2{4}, transparency)
plot(time_LB20a, mean(lcell_LB20a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)
% 
% ciplot(mean(lcell_LB20b, 1, 'omitnan')-std(lcell_LB20b, 0, 1, 'omitnan'), mean(lcell_LB20b, 1, 'omitnan')+std(lcell_LB20b, 0, 1, 'omitnan'), time_LB20b, colorcode2{6}, transparency)
% plot(time_LB20b, mean(lcell_LB20b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Length')

figure, hold on
ciplot(mean(lcell_LB1s, 1, 'omitnan')-std(lcell_LB1s, 0, 1, 'omitnan'), mean(lcell_LB1s, 1, 'omitnan')+std(lcell_LB1s, 0, 1, 'omitnan'), time_LB1s, colorcode2{1}, transparency)
plot(time_LB1s, mean(lcell_LB1s, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(lcell_LB2s, 1, 'omitnan')-std(lcell_LB2s, 0, 1, 'omitnan'), mean(lcell_LB2s, 1, 'omitnan')+std(lcell_LB2s, 0, 1, 'omitnan'), time_LB2s, colorcode2{5}, transparency)
plot(time_LB2s, mean(lcell_LB2s, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

ciplot(mean(lcell_LB3s, 1, 'omitnan')-std(lcell_LB3s, 0, 1, 'omitnan'), mean(lcell_LB3s, 1, 'omitnan')+std(lcell_LB3s, 0, 1, 'omitnan'), time_LB3s, colorcode2{3}, transparency)
plot(time_LB3s, mean(lcell_LB3s, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Length')

%Conclusion: lengths are the same for all full prime experiments. The cells
%in the shortened priming experiments are shorter in length but,
%percuilarly, they lengthen over time post-lysis.

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
% plot(time_LB1c, intensity_LBc, 'Color', colorcode{3});
% plot(time_LB1c, bgintensity_LB1c, '--k');

% figure, hold on
% plot(time_LB1c, intensity_LB1c, 'Color', colorcode{3});
% plot(time_LB20a, intensity_LB20a, 'Color', colorcode{4});
% plot(time_LB20b, intensity_LB20b, 'Color', colorcode{6});
% 
% figure, hold on
% plot(time_LB1a, intensity_LB1a, 'Color', colorcode{1});
% plot(time_LB1b, intensity_LB1b, 'Color', colorcode{5});
% plot(time_LB20a, intensity_LB20a, 'Color', colorcode{4});
% plot(time_LB20b, intensity_LB20b, 'Color', colorcode{6});

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
% ciplot(mean(intensity_LB1c, 1, 'omitnan')-std(intensity_LB1c, 0, 1, 'omitnan'), mean(intensity_LB1c, 1, 'omitnan')+std(intensity_LB1c, 0, 1, 'omitnan'), time_LB1c, colorcode2{3}, transparency)
% plot(time_LB1c, mean(intensity_LB1c, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% 
% ciplot(mean(intensity_LB20a, 1, 'omitnan')-std(intensity_LB20a, 0, 1, 'omitnan'), mean(intensity_LB20a, 1, 'omitnan')+std(intensity_LB20a, 0, 1, 'omitnan'), time_LB20a, colorcode2{4}, transparency)
% plot(time_LB20a, mean(intensity_LB20a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)
% 
% ciplot(mean(intensity_LB20b, 1, 'omitnan')-std(intensity_LB20b, 0, 1, 'omitnan'), mean(intensity_LB20b, 1, 'omitnan')+std(intensity_LB20b, 0, 1, 'omitnan'), time_LB20b, colorcode2{6}, transparency)
% plot(time_LB20b, mean(intensity_LB20b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
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
% figure, hold on
% plot(time_LB1a, adjintensity_LB1a, 'Color', colorcode{1});
% plot(time_LB1b, adjintensity_LB1b, 'Color', colorcode{5});
% plot(time_LB1c, adjintensity_LB1c, 'Color', colorcode{3});

% figure, hold on
% plot(time_LB1c, adjintensity_LB1c, 'Color', colorcode{3});
% plot(time_LB20a, adjintensity_LB20a, 'Color', colorcode{4});
% plot(time_LB20b, adjintensity_LB20b, 'Color', colorcode{6});
% 
% figure, hold on
% plot(time_LB1a, adjintensity_LB1a, 'Color', colorcode{1});
% plot(time_LB1b, adjintensity_LB1b, 'Color', colorcode{5});
% plot(time_LB20a, adjintensity_LB20a, 'Color', colorcode{4});
% plot(time_LB20b, adjintensity_LB20b, 'Color', colorcode{6});

% figure, hold on
% plot(time_LB1s, adjintensity_LB1s, '-r');
% plot(time_LB2s, adjintensity_LB2s, '-b');
% plot(time_LB3s, adjintensity_LB3s, '-g'); 

% figure, hold on
% ciplot(mean(adjintensity_LB1a, 1, 'omitnan')-std(adjintensity_LB1a, 0, 1, 'omitnan'), mean(adjintensity_LB1a, 1, 'omitnan')+std(adjintensity_LB1a, 0, 1, 'omitnan'), time_LB1a, colorcode2{1}, transparency)
% plot(time_LB1a, mean(adjintensity_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LB1b, 1, 'omitnan')-std(adjintensity_LB1b, 0, 1, 'omitnan'), mean(adjintensity_LB1b, 1, 'omitnan')+std(adjintensity_LB1b, 0, 1, 'omitnan'), time_LB1b, colorcode2{5}, transparency)
% plot(time_LB1b, mean(adjintensity_LB1b, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LB1c, 1, 'omitnan')-std(adjintensity_LB1c, 0, 1, 'omitnan'), mean(adjintensity_LB1c, 1, 'omitnan')+std(adjintensity_LB1c, 0, 1, 'omitnan'), time_LB1c, colorcode2{3}, transparency)
% plot(time_LB1c, mean(adjintensity_LB1c, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LB20a, 1, 'omitnan')-std(adjintensity_LB20a, 0, 1, 'omitnan'), mean(adjintensity_LB20a, 1, 'omitnan')+std(adjintensity_LB20a, 0, 1, 'omitnan'), time_LB20a, colorcode2{4}, transparency)
% plot(time_LB20a, mean(adjintensity_LB20a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)
% 
% ciplot(mean(adjintensity_LB20b, 1, 'omitnan')-std(adjintensity_LB20b, 0, 1, 'omitnan'), mean(adjintensity_LB20b, 1, 'omitnan')+std(adjintensity_LB20b, 0, 1, 'omitnan'), time_LB20b, colorcode2{6}, transparency)
% plot(time_LB20b, mean(adjintensity_LB20b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
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
% plot(tme_LB1c, normintensity_LB1c, 'Color', colorcode{3});
% 
% figure, hold on
% plot(tme_LB1c, normintensity_LB1c, 'Color', colorcode{3});
% plot(tme_LB20a, normintensity_LB20a, 'Color', colorcode{4});
% plot(tme_LB20b, normintensity_LB20b, 'Color', colorcode{6});
% 
% figure, hold on
% plot(tme_LB1a, normintensity_LB1a, 'Color', colorcode{1});
% plot(tme_LB1b, normintensity_LB1b, 'Color', colorcode{5});
% plot(tme_LB20a, normintensity_LB20a, 'Color', colorcode{4});
% plot(tme_LB20b, normintensity_LB20b, 'Color', colorcode{6});

figure, hold on
plot(tme_LB1s, normintensity_LB1s, '-r');
plot(tme_LB2s, normintensity_LB2s, '-b');
plot(tme_LB3s, normintensity_LB3s, '-g'); 

figure, hold on
ciplot(mean(normintensity_LB1a, 1, 'omitnan')-std(normintensity_LB1a, 0, 1, 'omitnan'), mean(normintensity_LB1a, 1, 'omitnan')+std(normintensity_LB1a, 0, 1, 'omitnan'), tme_LB1a, colorcode2{1}, transparency)
plot(tme_LB1a, mean(normintensity_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

ciplot(mean(normintensity_LB1b, 1, 'omitnan')-std(normintensity_LB1b, 0, 1, 'omitnan'), mean(normintensity_LB1b, 1, 'omitnan')+std(normintensity_LB1b, 0, 1, 'omitnan'), tme_LB1b, colorcode2{5}, transparency)
plot(tme_LB1b, mean(normintensity_LB1b, 1, 'omitnan'), 'Color', colorcode{5}, 'LineWidth', 1)

% ciplot(mean(normintensity_LB1c, 1, 'omitnan')-std(normintensity_LB1c, 0, 1, 'omitnan'), mean(normintensity_LB1c, 1, 'omitnan')+std(normintensity_LB1c, 0, 1, 'omitnan'), tme_LB1c, colorcode2{3}, transparency)
% plot(tme_LB1c, mean(normintensity_LB1c, 1, 'omitnan'), 'Color', colorcode{3}, 'LineWidth', 1)

% ciplot(mean(normintensity_LB20a, 1, 'omitnan')-std(normintensity_LB20a, 0, 1, 'omitnan'), mean(normintensity_LB20a, 1, 'omitnan')+std(normintensity_LB20a, 0, 1, 'omitnan'), tme_LB20a, colorcode2{4}, transparency)
% plot(tme_LB20a, mean(normintensity_LB20a, 1, 'omitnan'), 'Color', colorcode{4}, 'LineWidth', 1)

ciplot(mean(normintensity_LB20b, 1, 'omitnan')-std(normintensity_LB20b, 0, 1, 'omitnan'), mean(normintensity_LB20b, 1, 'omitnan')+std(normintensity_LB20b, 0, 1, 'omitnan'), tme_LB20b, colorcode2{6}, transparency)
plot(tme_LB20b, mean(normintensity_LB20b, 1, 'omitnan'), 'Color', colorcode{6}, 'LineWidth', 1)
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
saveas(gcf, 'tauFit_LB1s.png')

figure, hold on
plot(tme_LB2s, normintensity_LB2s, '-b')
plot(tme_LB2s, yhat_LB2s, '--k')
saveas(gcf, 'tauFit_LB2s.png')

figure, hold on
plot(tme_LB3s, normintensity_LB3s, '-g')
plot(tme_LB3s, yhat_LB3s, '--k')
saveas(gcf, 'tauFit_LB3s.png')

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

alpha=32.2114;
intercept=0.1614;

[Cnew_LB1a, dCB_LB1a, dCT_LB1a, dCP_LB1a, CblExp_LB1a, unbFrac_LB1a]=photoCorrect(tme_LB1a, normintensity_LB1a, alpha, intercept);
[Cnew_LB1b, dCB_LB1b, dCT_LB1b, dCP_LB1b, CblExp_LB1b, unbFrac_LB1b]=photoCorrect(tme_LB1b, normintensity_LB1b, alpha, intercept);
[Cnew_LB20a, dCB_LB20a, dCT_LB20a, dCP_LB20a, CblExp_LB20a, unbFrac_LB20a]=photoCorrect(tme_LB20a, normintensity_LB20a, alpha, intercept);

[Cnew_LB1s, dCB_LB1s, dCT_LB1s, dCP_LB1s, CblExp_LB1s, unbFrac_LB1s]=photoCorrect(tme_LB1s, normintensity_LB1s, alpha, intercept);
[Cnew_LB2s, dCB_LB2s, dCT_LB2s, dCP_LB2s, CblExp_LB2s, unbFrac_LB2s]=photoCorrect(tme_LB2s, normintensity_LB2s, alpha, intercept);
[Cnew_LB3s, dCB_LB3s, dCT_LB3s, dCP_LB3s, CblExp_LB3s, unbFrac_LB3s]=photoCorrect(tme_LB3s, normintensity_LB3s, alpha, intercept);

%% plot corrected traces
figure, plot(tme_LB1a, Cnew_LB1a, '-r');
figure, plot(tme_LB1b, Cnew_LB1b, '-r');
figure, plot(tme_LB20a, Cnew_LB20a, '-r'); 

figure, hold on
%ciplot(mean(Cnew_LB, 1, 'omitnan')-std(Cnew_LB, 0, 1, 'omitnan'), mean(Cnew_LB, 1, 'omitnan')+std(Cnew_LB, 0, 1, 'omitnan'), tme_LB, colorcode2{1}, transparency)
plot(tme_LB1a, mean(Cnew_LB1a, 1, 'omitnan'), 'Color', colorcode{1}, 'LineWidth', 1)

%ciplot(mean(Cnew_LB1c, 1, 'omitnan')-std(Cnew_LB1c, 0, 1, 'omitnan'), mean(Cnew_LB1c, 1, 'omitnan')+std(Cnew_LB1c, 0, 1, 'omitnan'), tme_LB1c, colorcode2{5}, transparency)
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
function [intensity, bgintensity, adjintensity, normintensity, lCell, time, tme, imstart]=dataNormalize1(datadir, set)
        
        %set 1=pre-lysis
        %set 2=post-lysis 
        
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
        
        if set==1
            %find the final pre-lysis frame by identifying the biggest
            %change in length
            dl=diff(lCell, 1, 2);
            lvg=mean(dl, 1, 'omitnan');
            [~, imstart]=min(lvg);
        else
            %find the initial post-lysis frame by identifying the first
            %time point where dt matches the final dt
            dt=round(diff(time), 2);
            if dt(1)==dt(end)
                dl=diff(lCell, 1, 2); 
                dlvg=mean(dl, 1, 'omitnan');
                [~, imstart] = min(dlvg);
                imstart=imstart+2;
            else
                if dt(end)<1 %discrepancies in dt might be found in the meta data, this is the ad hoc fix
                    imstart=min(find(dt==dt(end)))+2;
                else 
                    imstart=min(find(dt==dt(end)));
                end
            end
        end
        
        %remove cells without a value for imstart
        idx=find(~isnan(intensity(:, imstart)));
        intensity=intensity(idx,:);
        lCell=lCell(idx,:);    
        
        %subtract final fluor. value from the trace
        adjintensity=intensity-intensity(:, end); 
        adjintensity(adjintensity<200)=NaN;
            
        %find when all values are below the limit of detection
        [nrow, ncol]=size(adjintensity);
        nsum=sum(isnan(adjintensity), 1);
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
        
        %interpolate the fluor values during detergent perfusion
        for n=1:height(normintensity)
            idx=find(normintensity(n,:)>1);
            if ~isempty(idx)
                [nrow, ncol]=size(normintensity);
                x=setdiff(1:ncol, idx); %x=time
                v=normintensity(n, x); %v=intensity
                vq=interp1(x, v, idx); %vq=interpolation at query time pts
                normintensity(n,idx)=vq;
            end
        end
        
end

%do not fit to y=(1-beta)*exp(-t./tau)+beta (beta=normalized limit of detection)
%fit to y=Ae^(-t/tau)
function [tau, yhat]=tauCalc(tme, normintensity)
    
    [nrow, ~]=size(normintensity);
    tau=nan(nrow, 1);
    modelfun = @(tau, x)exp(-x./tau);
    for i=1:nrow
        tau(i, 1)=nlinfit(tme, normintensity(i,:), modelfun, 1);
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
        Cnew(:, 1)=normintensity(:, 1);
        
        
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
