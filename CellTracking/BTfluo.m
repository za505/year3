%BTfluo.m
%Rico Rojas, updated 1/21/19
%Zarina Akbary, updated 3/11/21
%Calculates the average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack.m.  

clear, close all

%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save fluorescent image stacks
%directories by themselves. 

%INPUT
%basename: name to save the results to.
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%cellIntensity=vector with raw intensity values
%bgIntensity=vector with raw background intensity values
%Iin=cellular intensity values adjusted for autofluorescence
%Iout=background autointensity values adjusted for background

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='04222021_Exp1_colony1';
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/04222021_analysis/' basename];
savename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/04222021_analysis/' basename '/' basename '_FSS/'  basename '_figures'];;
channel=[filename '/' basename '_FSS/' basename '_aligned'];

frameSwitch=118; %this is the initial frame for the switch (after which we pulse w/ and w/out Mg2+)
frameInitial=8; %this is the FIRST frame that has dye without lysing the membrane
frameAuto=105; %this is the LAST frame that has dye without lysing the membrane
frameBg=105; %this is the frame that you'll pick the background area from
recrunch=0; %recrunch=1, just redo the plots, don't recalculate any values
%xlabels=["PBS + FSS" "PBS + FSS + 9 mM Mg" "PBS + FSS"];
%xswitch=[300 360 420];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    
    load ([filename '/' basename '_FSS/' basename '_figures/' basename '_BTfluo'])
    
%     xlabels=["PBS + 5% detergent" "PBS + FSS" "PBS + FSS + 9 mM Ca" "PBS + FSS"];
%     xswitch=[60 300 410 530];
    
elseif recrunch==0

%load fluorescent images
curdir=cd;
cd(channel); 
fluo_directory=dir('*.tif');

%load BacTrack data
load([filename '/' basename '_phase/' basename '_figures/' basename '_BTphase'], 'T', 'time', 'tmid', 'pixels', 'xlabels', 'xswitch')

%preallocate cells
cellIntensity=[]; %cellular intensities
avgIntensity=zeros(1,T); %population average cellular intensity
stdIntensity=zeros(1,T); %std of population cellular intensity

bgIntensity=zeros(1,T); %background intensity
cellAuto=zeros(height(pixels),frameAuto-frameInitial); %intensity of cells prior to dye perfusion, aka autofluorescence
bgAuto=zeros(1,frameSwitch-frameAuto-1); %intensity of background prior to dye perfusion, aka autofluorescence

%determine region where you'll measure background intensity
imagename=fluo_directory(frameBg).name;
im=imread(imagename);    
[p1, p2]=getBackground(imagename);
close all

%now let's track intensity over time
for t=1:T

    t

    %load the image
    imagename=fluo_directory(t).name;
    im=imread(imagename);
    
    %measure background level
    bglevel = measureBackground(imagename, p1, p2);
    bgIntensity(t)=bglevel; 
    
    %calculate background autofluorescence
    if t <= frameSwitch & t > frameAuto
        bgAuto(1, t)=measureBackground(imagename, p1, p2);
    else 
        bgAuto(1, t)=NaN;
    end
    
    %measure cellular intensity
    for j=1:height(pixels)
         
        %calculate intensity
        cellIntensity(j,t)=mean(im(pixels{j,t}));
        
        %calculate cellular autofluorescence
        if t <= frameAuto & t >= frameInitial
            cellAuto(j, t)=mean(im(pixels{j,t}));
        else
            cellAuto(j, t)=NaN;
        end
        
    end

end
    
%cellAuto(cellAuto==0)=NaN;
        
%take the population average of the autofluorescence 
cellAuto=mean(cellAuto, 'omitnan');

%now take the temporal average of the background and cellular
%autofluorescence
cellAdj=mean(cellAuto, 'all', 'omitnan');
bgAdj = mean(bgAuto, 'all', 'omitnan');

%now subtract the background autofluorescence from background intensity
Iout = bgIntensity - bgAdj;
Iout(Iout<=0)=0; 

%now subtract the background autofluorescence from background intensity
Iin = cellIntensity - cellAdj;
Iin(Iin<=0)=0;

%calculate average intensity and standard deviation
avgIntensity = mean(Iin, 'omitnan');
stdIntensity = std(Iin, 'omitnan');

%now, calculate the Iin/Iout ratio
ratio = Iin ./ Iout;
ratio(:,1:frameSwitch-1)=NaN; %note: 0/0 makes weird things happen

avgRatio = mean(ratio, 1, 'omitnan');
stdRatio = std(ratio, 0, 1, 'omitnan');

%calculate Intencity Rate
deltat=time(2:end)-time(1:end-1);
Irate=(ratio(:, 2:end)-ratio(:, 1:end-1))./((ratio(:, 1:end-1)+ratio(:, 2:end))/2);
for i=1:height(Irate)
    Irate(i, :)=Irate(i, :)./deltat;
end

avgIrate = mean(Irate, 1, 'omitnan');
stdIrate = std(Irate, 0, 1, 'omitnan');

end

%let's change folders to save the plots and variables
cd([filename '/' basename '_FSS/' basename '_figures'])
cd(savename)

%save the variables
save([basename '_BTfluo'])

%Plot data
%Let's plot cellular intensity first
figure, hold on, title('Intensity vs Time')
for i=1:height(pixels)
    plot(time,cellIntensity(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
for x=1:length(xlabels)
    xline(xswitch(x), '--k', xlabels(x)) 
end
ylim([-3 Inf])
saveas(gcf, [basename '_intensity.png'])

%now population average cell intensity
figure
ciplot(avgIntensity - stdIntensity, avgIntensity + stdIntensity, time, [0.75 0.75 1])
plot(time, avgIntensity,'-r')
title('Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
for x=1:length(xlabels)
    xline(xswitch(x), '--k', xlabels(x)) 
end
ylim([-3 Inf])
saveas(gcf, [basename,'_intensityAvg.png'])

%Plot average background fluorescence
figure, hold on 
title('Background Intensity vs Time')
plot(time,Iout)
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
for x=1:length(xlabels)
    xline(xswitch(x), '--k', xlabels(x)) 
end
ylim([-3 Inf])
hold off
saveas(gcf, [basename,'_background.png'])

%Plot the Iin/Iout ratio over time
figure, hold on
title('Intensity/Background vs Time')
ciplot(avgRatio - stdRatio, avgRatio + stdRatio, time, [0.75 0.75 1])
plot(time,avgRatio)
xlabel('Time (s)')
ylabel('Intensity/Background')
fig2pretty 
yline(1, '--k')
for x=1:length(xlabels)
    xline(xswitch(x), '--k', xlabels(x)) 
end
ylim([0 Inf])
hold off
saveas(gcf, [basename,'_ratioTime.png'])

%Plot Intensity Rate
figure, hold on
ciplot(avgIrate - stdIrate, avgIrate + stdIrate, tmid, [0.75 0.75 1])
plot(tmid,avgIrate)
xlabel('Time (s)')
ylabel('Intensity Rate (s^{-1})')
fig2pretty
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
ylim([0 Inf])
xlim([0 Inf])
hold off
saveas(gcf, [basename,'_intensityRate.png'])

%%%%%Functions
function [p1, p2]=getBackground(imagename)
        
        %Load last image
        %imagename=fluo_directory{i}(t).name;
        im2=imread(imagename);

        %Determine Background
        figure,imshow(im2,[]), hold on, title('Select Background')
        k=waitforbuttonpress;
        set(gcf,'Pointer')
        hold on
        axis manual
        point1=get(gca,'CurrentPoint');
        finalRect=rbbox;
        point2=get(gca,'CurrentPoint');
        point1=point1(1,1:2);
        point2=point2(1,1:2);
        point1(point1<1)=1;
        point2(point2<1)=1;
        p1=min(point1,point2);%Calculate locations
        p2=max(point1,point2);
        offset = abs(point1-point2);%And dimensions
        x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
        y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
        plot(x,y)
        p1=round(p1);
        p2=round(p2);  
end 

function bglevel = measureBackground(imagename, p1, p2)
        
        %Load last image
        %imagename=fluo_directory{i}(t).name;
        im2=imread(imagename);
        
        %Determine background
        backim=im2(p1(2):p2(2),p1(1):p2(1));
        [counts,bins]=imhist(backim);
        [~,binnum]=max(counts);
        maxpos=bins(binnum);
        bglevel=mean(mean(backim));
        
end 