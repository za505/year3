%analysis_FSS.m

clear, close all

%% user input
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/FSS_diffusion_analysis";
dirsave="/Users/zarina/Downloads/NYU/Year3_2021_Fall/presentations/09062021";

%diffC=Cin-Cout
data = struct('experiment', [], 'length', [], 'ncells', [], 'diffC', [], 'time', []);
cd(datadir)
cutoff=400;
pulse=60;

%% 04202021 Mg2+ pulse
data(1).experiment = '04202021_Exp1';
load("04202021_Exp1_colony1_BTphase.mat", 'l', 'ncells', 'time')
load("04202021_Exp1_colony1_BTfluo.mat", 'diffC')

data(1).length=l;
data(1).ncells=ncells;
data(1).diffC=diffC;
data(1).time=time;

%% 04222021 FSS pulse
data(2).experiment = '04222021_Exp1';
load("04222021_Exp1_colony1_BTphase.mat", 'l', 'ncells', 'time')
load("04222021_Exp1_colony1_BTfluo.mat", 'diffC')

data(2).length=l;
data(2).ncells=ncells;
data(2).diffC=diffC;
data(2).time=time;

%% compare the difference in length contractions
% for i=1:length(data)
%     idx = find(data(i).time >= 400);
%     figure, hold on
%     if i==1
%         plot()
% end

figure(1), hold on
for n=1:data(1).ncells
    idx = find(data(1).time >= cutoff);
    %plot(data(1).diffC(n, idx), data(1).length(n, idx))
    plot(data(1).time(idx), data(1).length(n, idx))
end

figure(2), hold on
for n=1:data(2).ncells
    idx = find(data(2).time >= cutoff);
   plot(data(2).diffC(n, idx), data(2).length(n, idx))
end