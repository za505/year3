%BTfluoAVG.m Zarina Akbary, 03/25/21 based on BTfluo by Rico Rojas, updated
%1/21/19 Calculates the average cytoplasmic fluorescence intensity and
%intensity/background ratio from cell tracked with BacTrack.m and BTfluo.m

clear, close all

%INSTRUCTIONS FOR USE:
%Align phase and fluorescent stacks with imagealign.m, remove pillars from
%the phase image with eraseimage.m, and then analyze using BacTrack.m.
%Calculate fluorescent intensities using BTfluo.m. That data from that
%script will be load here.

%INPUT
%basename: name to obtain data from.
%channels: list of directories containing fluorescent image stacks.

%OUTPUT:



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename=["04062021_Exp4_colony1", "04062021_Exp4_colony2", "04062021_Exp4_colony3"];%Name of the image stack, used to save file.
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/04062021_analysis/04062021_Exp4'];
channel=['_FSS'];
recrunch=0;
B=length(basename);%number of main directories to analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==0
    
%load the time variables
base=char(basename(1))
cd([filename '/' base '/' base channel '/' base '_figures'])
filelist{1}=dir([base '_BTfluo.mat']);
T = cell2mat(struct2cell(load([filelist{1}.name],'T')));
time = cell2mat(struct2cell(load([filelist{1}.name],'time')));
time = cell2mat(struct2cell(load([filelist{1}.name],'tmid')));
frameAuto = cell2mat(struct2cell(load([filelist{1}.name],'frameAuto')));
frameInitial = cell2mat(struct2cell(load([filelist{1}.name],'frameInitial')));
xlabels = struct2cell(load([filelist{1}.name],'xlabels'));
xswitch = struct2cell(load([filelist{1}.name],'xswitch'));

xlabels = xlabels{1,1};
xswitch = xswitch{1,1};

%pre-allocate matricies
%there we just want to have in one place to see
cellAdj=zeros(B,1);
bgAdj=zeros(B,1);
bgIntensity=zeros(B,T);
cellIntensity=[];
ratio=[];
Iout=zeros(B, T);
Iin=[];
Irate=[];

%these are the ones we'll be re-calculating/plotting
avgIin=zeros(1, T);
stdIin=zeros(1, T);
avgRatio=zeros(1, T);
stdRatio=zeros(1, T);

for b=1:B
    
    base=char(basename(b))
    cd([filename '/' base '/' base channel '/' base '_figures'])
    filelist{b}=dir([base '_BTfluo.mat']);
    
    bgAdj(b,1) = cell2mat(struct2cell(load([filelist{b}.name],'bgAdj')));
    cellAdj(b,1) = cell2mat(struct2cell(load([filelist{b}.name],'cellAdj')));
    bgIntensity(b, :) = cell2mat(struct2cell(load([filelist{b}.name],'bgIntensity')));
    Iout(b, :) = cell2mat(struct2cell(load([filelist{b}.name],'Iout')));
    temp=struct2cell(load([filelist{b}.name],'cellIntensity'));
    ratio_temp=struct2cell(load([filelist{b}.name],'ratio'));
    Iin_temp=struct2cell(load([filelist{b}.name],'Iin'));
    Irate_temp=struct2cell(load([filelist{b}.name],'Irate'));
    
    %let's get all the raw intensities, adj intensities, and ratios in one place
    temp=temp{1,1};
    cellIntensity=[cellIntensity; temp];
    
    ratio_temp=ratio_temp{1,1};
    ratio=[ratio; ratio_temp];
    
    Iin_temp=Iin_temp{1,1};
    Iin=[Iin; Iin_temp];
    
    Irate_temp=Irate_temp{1,1};
    Irate=[Irate; Irate_temp];
  
end
    
%let's calulcate the population average intensity based on all the cells
avgIin = mean(Iin, 'omitnan');
stdIin = std(Iin, 'omitnan');

%let's calulcate the population average background based on all the cells
avgIout = mean(Iout, 'omitnan');
stdIout = std(Iout, 'omitnan');

%remove ratios from frames without dye
ratio(:,frameInitial:frameAuto)=NaN;

%let's calulcate the population average based on all the cells
avgRatio = mean(ratio, 'omitnan');
stdRatio = std(ratio, 'omitnan');

%let's calulcate the population average intensity rate 
avgRate = mean(Irate, 'omitnan');
stdRate = std(Irate, 'omitnan');

%change directory
cd(filename);
save(['BTfluoAVG' channel])

else
    %change directory
    cd(filename);
    load(['BTfluoAVG' channel])
end

%let's plot the population average intensity
figure, hold on
title('Population Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
for x=1:length(xlabels)
    xline(xswitch(x), '--k', xlabels(x)) 
end
ylim([-2 Inf])
ciplot(avgIin - stdIin, avgIin + stdIin, time, [0.75 0.75 1])
plot(time, avgIin, '-r')
saveas(gcf, ['avgIin' channel '.png'])

%let's plot the population average background intensity
figure, hold on
title('Average Background Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
for x=1:length(xlabels)
    xline(xswitch(x), '--k', xlabels(x)) 
end
ylim([-2 Inf])
ciplot(avgIout - stdIout, avgIout + stdIout, time, [0.75 0.75 1])
plot(time, avgIout, '-r')
saveas(gcf, ['avgIout' channel '.png'])

%let's plot the average ratio too
figure, hold on
title('Intensity/Background vs Time')
xlabel('Time (s)')
ylabel('Intensity/Background')
ylim([0 Inf])
fig2pretty
yline(1, '--k')
for x=1:length(xlabels)
    xline(xswitch(x), '--k', xlabels(x)) 
end
ciplot(avgRatio - stdRatio, avgRatio + stdRatio, time, [0.75 0.75 1])
plot(time, avgRatio, '-r') 
saveas(gcf, ['avgRatio' channel '.png'])

%and finally average intensity rate
figure, hold on
title('Intensity Rate vs Time')
xlabel('Time (s)')
ylabel('Intensity/Background')
ylim([0 Inf])
fig2pretty
yline(1, '--k')
for x=1:length(xlabels)
    xline(xswitch(x), '--k', xlabels(x)) 
end
ciplot(avgIrate - stdIrate, avgIrate + stdIrate, time, [0.75 0.75 1])
plot(tmid, avgIrate, '-r') 
saveas(gcf, ['avgRatio' channel '.png'])


