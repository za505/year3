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
basename='05022021_Exp1_colony1';
filename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05022021_analysis/'];
savename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05022021_analysis/' basename '_FITC/'  basename '_figures'];;
channel=[filename '/' basename '_FITC/' basename '_aligned'];

frameSwitch=120; %this is the initial frame for the switch (after which we pulse w/ and w/out Mg2+)
frameInitial=14; %this is the FIRST frame that has dye without lysing the membrane
frameAuto=93; %this is the LAST frame that has dye without lysing the membrane
frameBg=90; %this is the frame that you'll pick the background area from
recrunch=0; %recrunch=1, just redo the plots, don't recalculate any values
vis=1;
outThresh=18000; %intensity cutoff; x boundary 
%xlabels=["PBS + FITC" "PBS + FITC + 9 mM Mg" "PBS + FITC"];
%xswitch=[300 360 420];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    
    load ([filename '/' basename '_FITC/' basename '_figures/' basename '_BTfluo'])
    
%     xlabels=["PBS + 5% detergent" "PBS + FITC" "PBS + FITC + 9 mM Ca" "PBS + FITC"];
%     xswitch=[60 300 410 530];
    
elseif recrunch==0

%load fluorescent images
curdir=cd;
cd(channel); 
fluo_directory=dir('*.tif');

%load BacTrack data
load([filename '/' basename '_phase/' basename '_figures/' basename '_BTphase'], 'time', 'T','directory', 'dirname', 'xlabels', 'xswitch')
load([filename '/' basename '_phase/' basename '_figures/' basename '_BTlab'],'labels', 'labels2')

%calculate tmid
tmid=(time(2:end)+time(1:end-1))/2;

%preallocate variables
cellnumber=zeros(1, T);

%Let's get the coordinates of the cells for each image
for t=1:T
    
    %Let's calculate the number of cells in each frame
    cellnumber(t)=max(max(labels(:,:,t)));
    
   %now let's get the row and column coordinates of the cells in each frame
    for n=1:cellnumber(t)
    
        [r1,c1] = find(labels(:,:,t)==n);   %Find coordinates of cells in frame T

        r1c1{n,t}=[r1 c1]; % concatenate XY coords frame T 
        
    end
    
end

%calculate max total cells
maxcellnumber=mode(cellnumber(frameInitial:frameAuto)); %if we can't correct for it, throw it out
    
%pre-allocate more variables
pixelIntensity=cell(maxcellnumber, T); %cellular intensities
avgIntensity=zeros(1,T); %population average cellular intensity
stdIntensity=zeros(1,T); %std of population cellular intensity

bgIntensity=zeros(1,T); %background intensity

Cin=nan(maxcellnumber, T);
Cexp=nan(maxcellnumber, T);
cellIntensity=nan(maxcellnumber, T);
intercept=nan(maxcellnumber, 1);
slope=nan(maxcellnumber, 1);
Cout=[]; %background autofluorescence
Cauto=[];

stats=cell(maxcellnumber, 1); %where all the models are stored

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
    
    %measure cellular intensity
    for n=1:maxcellnumber
        
        %make a separate variable to store the pixel intensities of 1 cell
        cellTemp=nan(height(r1c1{n,t}),1);
        
        %get the intensity of each pixel
        for p=1:height(cellTemp)
         cellTemp(p)=im(r1c1{n,t}(p,1),r1c1{n,t}(p,2));
        end
        
        %store this in the cellIntensity variable
        pixelIntensity{n,t}=cellTemp;
       
        %take the mean and calculate Cin
        Cin(n,t)=mean(pixelIntensity{n,t});
        
    end
    
    %index Cout and Cin
    if t>=frameInitial & t<=frameAuto & bgIntensity(t)<=outThresh
        Cout=[Cout bgIntensity(t)];
        Cauto=[Cauto Cin(:, t)];
    end
   
     %check that the coordinates are correct
    if vis==1 & (t<=6 | t>=T-6)
        figure
        imshow(im)
        hold on
        for n=1:cellnumber(t)
           plot(r1c1{n,t}(:, 2),r1c1{n,t}(:, 1),'-r')
        end
        pause
        close
    end

end
    
%now let's fit a line to the Cout and Cauto fluorescence data
for n=1:maxcellnumber
    stats{n,1}=fitlm(Cout, Cauto(n,:));
    arrayTemp=table2array(stats{n,1}.Coefficients);
    intercept(n,1)=arrayTemp(1,1);
    slope(n,1)=arrayTemp(2,1);
    
    for t=1:T
        Cexp(n,t)=slope(n,1)*bgIntensity(1,t)+intercept(n,1);
        cellIntensity(n,t)=Cin(n,t)-Cexp(n,t);
    end
    
end

%calculate average intensity and standard deviation
avgIntensity = mean(cellIntensity, 'omitnan');
stdIntensity = std(cellIntensity, 'omitnan');

%now, calculate the Iin/Iout ratio
ratio(:,frameInitial:frameAuto)=NaN; %note: 0/0 makes weird things happen
ratio = Cin ./ bgIntensity;

avgRatio = mean(ratio, 1, 'omitnan');
stdRatio = std(ratio, 0, 1, 'omitnan');

%calculate Ratio Rate
deltat=time(2:end)-time(1:end-1);
Irate=(ratio(:, 2:end)-ratio(:, 1:end-1))./((ratio(:, 1:end-1)+ratio(:, 2:end))/2);
for i=1:height(Irate)
    Irate(i, :)=Irate(i, :)./deltat;
end

avgIrate = mean(Irate, 1, 'omitnan');
stdIrate = std(Irate, 0, 1, 'omitnan');

end

%let's change folders to save the plots and variables
cd(savename)

%save the variables
save([basename '_BTfluo'])

%Plot data
%Let's plot cellular intensity first
figure, hold on, title('Intensity vs Time')
for i=1:maxcellnumber
    plot(time,cellIntensity(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
ylim([-3 Inf])
% saveas(gcf, [basename '_intensity.png'])

%now population average cell intensity
figure
ciplot(avgIntensity - stdIntensity, avgIntensity + stdIntensity, time, [0.75 0.75 1])
plot(time, avgIntensity,'-r')
title('Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
ylim([-3 Inf])
% saveas(gcf, [basename,'_intensityAvg.png'])

%Plot average background fluorescence
figure, hold on 
title('Background Intensity vs Time')
plot(time,bgIntensity)
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
ylim([-3 Inf])
hold off
% saveas(gcf, [basename,'_background.png'])

%Plot the Iin/Iout ratio over time
figure, hold on
title('Intensity/Background vs Time')
ciplot(avgRatio - stdRatio, avgRatio + stdRatio, time, [0.75 0.75 1])
plot(time,avgRatio)
xlabel('Time (s)')
ylabel('Intensity/Background')
fig2pretty 
yline(1, '--k')
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
ylim([0 Inf])
hold off
% saveas(gcf, [basename,'_ratioTime.png'])

%Plot Intensity Rate
figure, hold on
ciplot(avgIrate - stdIrate, avgIrate + stdIrate, tmid, [0.75 0.75 1])
plot(tmid,avgIrate)
xlabel('Time (s)')
ylabel('Ratio Rate (s^{-1})')
fig2pretty
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
ylim([0 Inf])
xlim([0 Inf])
hold off
% saveas(gcf, [basename,'_intensityRate.png'])

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