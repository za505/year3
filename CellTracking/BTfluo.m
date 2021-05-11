%BTfluo.m
%Rico Rojas, updated 1/21/19
%Zarina Akbary, updated 3/11/21
%Calculates the average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack.m.  

clear, close all

%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save fluorescent image stacks
%directories by themselves. Need the Signal Processing Toolbox

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
savename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/04222021_analysis/' basename '/' basename '_FSS/'  basename '_figures'];
labelname=[filename '/' basename '_phase/' basename '_figures/' basename '_BTlab'];
channel=[filename '/' basename '_FSS/' basename '_aligned'];

frameSwitch=202; %this is the initial frame for the switch (after which we pulse w/ and w/out Mg2+)
frameInitial=8; %this is the FIRST frame that has dye without lysing the membrane
frameAuto=105; %this is the LAST frame that has dye without lysing the membrane
frameBg=104; %this is the frame that you'll pick the background area from
recrunch=1; %recrunch=1, just redo the plots, don't recalculate any values
calibration=0; %this variable is to save time while I troubleshoot the code
vis=0;
xlabels=["PBS + FSS" "PBS + FSS + 9 mM Mg" "PBS + FSS"];
xswitch=[300 360 420];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    
    load ([savename '/' basename '_BTfluo'])
       
elseif recrunch==0

%load fluorescent images
curdir=cd;
cd(channel); 
fluo_directory=dir('*.tif');

%load BacTrack data
load([filename '/' basename '_phase/' basename '_figures/' basename '_BTphase'], 'time', 'T','directory', 'dirname', 'xlabels', 'xswitch')

%calculate tmid
tmid=(time(2:end)+time(1:end-1))/2;

%preallocate variables
cellnumber=zeros(1, T);

%%This code takes a bit of time, so we can make a short cut for
%%troubleshooting
if calibration==0
    load(labelname,'labels')
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
else
        load ([filename '/' basename '_FSS/' basename '_figures/' basename '_BTfluo'], 'cellnumber', 'r1', 'c1', 'r1c1')
end

%calculate max total cells
maxcellnumber=mode(cellnumber(frameInitial:frameAuto)); %we don't want to track cells that show up only once or twice
    
%pre-allocate more variables
pixelIntensity=cell(maxcellnumber, T); %pixel intensities for each cell in each frame

bgAuto=[]; %background autofluorescence
%cellAuto=[]; %cellular autofluorescence
%flag=1;
bgIntensity=nan(1,T); %background intensity for all the frames
cellIntensity=nan(maxcellnumber, T); %raw pixel intensity means
Cin=nan(maxcellnumber, T); %adjusted cell intensities for each cell in each frame
Cout=nan(1, T); %adjusted background intensity in each frame

avgIntensity=nan(1,T); %population average cellular intensity
stdIntensity=nan(1,T); %std of population cellular intensity
%%
%determine region where you'll measure background intensity
imagename=fluo_directory(frameBg).name;
im=imread(imagename);    
[p1, p2]=getBackground(imagename);
close all

 %prevent blips in cellular intensity
  %this smooths the data for a given cell at a given time point 
 for n=1:maxcellnumber
    
     %get the number of pixels in each time point for a given cell
     heightTemp=cellfun(@height, r1c1(n,1:T));

    %this smooths the data for a given cell at a given time point
    ping=[];
    for t=1:T-2
        lowThresh=heightTemp(1,t)*0.9;
        highThresh=heightTemp(1,t)*1.1;
        ind1=find(heightTemp(1, t:t+1)<lowThresh | heightTemp(1, t:t+1)>highThresh);
        ind2=find(heightTemp(1, t:t+2)<lowThresh | heightTemp(1, t:t+2)>highThresh);
            if length(ind1)==1 & length(ind2)==1
                ping=[ping t+1];
                r1c1{n,t+1}=r1c1{n,t};
            end
    end
    ping;
 end
 
 %make sure each cell is unique and not double counted
  for n=flip(2:maxcellnumber)
    for t=2:T
        match=setdiff(r1c1{n,t}, r1c1{n-1,t});
        if isempty(match)==1
            r1c1{n,t}=r1c1{n,t-1};
        end
    end 
  end
 

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
       
        %take the mean and raw cell. intensity for each cell at each frame
        cellIntensity(n,t)=mean(pixelIntensity{n,t});
        
    end
    
    %get background fluorescence
    if t>frameAuto & t<frameSwitch 
        bgAuto=[bgAuto bgIntensity(t)];
    end
    
     %check that the coordinates are correct
    if vis==1 & (t<=6 | t>=T-6)
        figure
        imshow(im, [])
        hold on
        for n=1:cellnumber(t)
           plot(r1c1{n,t}(:, 2),r1c1{n,t}(:, 1),'-r') %x is im 'column', y is im 'row'
        end
        pause
        close
    end

end
    
%%
%now take the mean of bgAuto
bgAuto=mean(bgAuto,2); %take the temporal average

%Subtract the autofluorescence from the raw intensity
Cout=bgIntensity-bgAuto;
for n=1:maxcellnumber
    Cin(n,:)=cellIntensity(n,:)-bgAuto;
end
%%
%calculate average intensity and standard deviation
avgIntensity = mean(Cin, 1, 'omitnan');
stdIntensity = std(Cin, 0, 1, 'omitnan');

%now, calculate the Iin/Iout ratio
ratio = Cin ./ Cout;
ratio(:,[1:frameInitial-1, frameAuto+1:frameSwitch-1])=NaN; %note: 0/0 makes weird things happen

avgRatio = mean(ratio, 1, 'omitnan');
stdRatio = std(ratio, 0, 1, 'omitnan');

%calculate the diffusion constant
deltat=[0 time(2:end)-time(1:end-1)];
dCin=[zeros(maxcellnumber,1) Cin(:,2:end)-Cin(:,1:end-1)];
diffC=nan(maxcellnumber,T);
D=nan(maxcellnumber,T);
for n=1:maxcellnumber
    dCin(n, :)=dCin(n, :)./deltat;
    diffC(n, :)=Cin(n, :)-Cout;
    D(n,:)=dCin(n,:)./diffC(n,:);
end

avgIrate = mean(dCin, 1, 'omitnan');
stdIrate = std(dCin, 0, 1, 'omitnan');

end

%let's change folders to save the plots and variables
cd(savename)

%save the variables (we are going to exclude saving labels, which requires
%a lot of storage
clear labels
save([basename '_BTfluo'])

%Plot data
%Let's plot cellular intensity first
figure, hold on, title('Adjusted Intensity vs Time')
for i=1:maxcellnumber
    plot(time,Cin(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
ylim([-200 Inf])
saveas(gcf, [basename '_intensityRaw.png'])

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
saveas(gcf, [basename,'_intensityAvg.png'])

%Plot adj background fluorescence
figure, hold on 
title('Adjusted Background Intensity vs Time')
plot(time,Cout)
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
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
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
%ylim([0 Inf])
hold off
saveas(gcf, [basename,'_ratioTime.png'])

% %Plot Intensity Rate
figure, hold on
title('Average Intensity Rate vs Time')
ciplot(avgIrate - stdIrate, avgIrate + stdIrate, time, [0.75 0.75 1])
plot(time,avgIrate)
xlabel('Time (s)')
ylabel('Intensity Rate')
fig2pretty 
% for x=1:length(xlabels)
%     xline(xswitch(x), '--k', xlabels(x)) 
% end
ylim([0 Inf])
hold off
saveas(gcf, [basename,'_intensityRate.png'])

%Let's plot Cin-Cout
figure, hold on, title('Cin-Cout vs Time')
for n=1:maxcellnumber
    plot(time,Cin(n,:)-Cout)
end
hold off
saveas(gcf, [basename,'_dConstant.png'])

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