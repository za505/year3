%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all
%re-analysis of 06/02/2021 data
%% first, I want to plot how intensity changes overtime
%the change when switching from PBS + 5% NLS to *PBS [+ Mg2]
dirname=['/Users/zarina/Documents/MATLAB/MatlabReady/06022021_analysis'];
cd(dirname)

datadir=dir('*.mat');

experiment=[];
exp=[];
flowR=[];
frameR=[];
intensity=[];

for i=1:height(datadir)
    cd(datadir(i).folder);
    filename=datadir(i).name;
    load(filename)
    experiment=[experiment; basename];
    exp=[exp; exposure];
    flowR=[flowR; flowRate];
    frameR=[frameR; frameRate];
    intensity=[intensity; intensityAvg];
end

data=table('experiment', experiment, 'exposure', exp,'flowRate', flowR,'frameRate', frameR, 'intensity', intensity);
data={experiment, exp, flowR, frameR, intensity};
figure, hold on
for i=1:height(datadir)
    plot(time, intensity(i,:))
end
title('Background Fluor. Intensity')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
legend({'exp: 10 ms, flow: 5 psi', 'exp: 10 ms, flow: 8 psi', 'exp: 10 ms, flow: 10 psi', 'exp: 20 ms, flow: 5 psi', 'exp: 20 ms, flow: 8 psi', 'exp: 20 ms, flow: 10 psi', 'exp: 40 ms, flow: 5 psi', 'exp: 40 ms, flow: 8 psi', 'exp: 40 ms, flow: 10 psi'})