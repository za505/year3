%BTfluo.m
%Rico Rojas, updated 1/21/19
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
%icell: Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a matrix (ncellxT) with the fluorescent intensities of each
        %cell, where rows are the cells and columns are time points.
%icell_av:  Cell array with length equal to the number of fluorescent
        %channels.  Each entry is a vector containing the population-
        %average of the single-cell fluorescent intensities.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='05082021_Exp5';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05082021_reanalysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory that the BTphase.mat file is saved in
savedir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05082021_reanalysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05082021_reanalysis/' basename '_colony1/' basename '_FITCK/' basename '_aligned']}; 
recrunch=0;
frameBg=20; %this is the frame that you'll pick the background area from
multiScale=0;
troubleshoot=0;
fluorFrames=[14:92,104:687]; %frames where FITC is perfused
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==0;
    
curdir=cd;
for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

cd(dirname)
load([basename '_BTphase'], 'pixels', 'ncells', 'lcell', 'B', 'time')

T=length(fluo_directory{1});
icell=cell(length(channels),1);
icell_av=cell(length(channels),1);
ratio=cell(length(channels),1);
ratio_av=cell(length(channels),1);
bgIntensity=nan(length(channels), T);

for i=1:length(channels)
    
    cd(channels{i}); 
    intensities_temp=zeros(ncells, T);
    imagename=fluo_directory{i}(frameBg).name;
    [p1, p2]=getBackground(imagename);
    
    for t=1:T
        t
        imagename=fluo_directory{i}(t).name;
        im=imread(imagename);
        for j=1:ncells
            intensities_temp(j,t)=mean(im(pixels{j,t}));    
        end
        
         %measure background level
         bglevel = measureBackground(imagename, p1, p2);
         bgIntensity(i,t)=bglevel; 
    end
    
    intensities_temp(intensities_temp==0)=NaN;
    icell{i}=intensities_temp;
end

for i=1:length(channels)
    for j=1:ncells
        ratio{i}(j,:)=icell{i}(j,:)./bgIntensity(i,:);
    end
end

for i=1:length(channels)
    icell_av{i}=nanmean(icell{i},1);
    ratio_av{i}=nanmean(ratio{i},1);
end

cd(savedir)
save([basename '_BTfluo'])

elseif recrunch==1
    load ([basename '_BTfluo'])
end

%% Troubleshooting
cd(channels{1}); 
 if troubleshoot==1
     for t=1:T
            t
            imagename=fluo_directory{1}(t).name;
            im=imread(imagename);


           figure
           imshow(im, [])
           hold on
           for k=1:ncells
                if isempty(B{k,t})==0
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
                else
                    continue
                end
           end
          pause
          close all
     end
 elseif troubleshoot==2
      for k=1:ncells
            k
            imagename=fluo_directory{1}(T).name;
            im=imread(imagename);


           figure
           imshow(im, [])
           hold on
           for t=1:T
                if isempty(B{k,t})==0
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
                else
                    continue
                end
           end
          pause
          close all
      end
 else
 end
%% Plot FITCK data
cd(savedir)

figure(1), hold on, 
for i=1:ncells
    plot(time, icell{1}(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
xlim([-0.2 Inf])
fig2pretty
saveas(gcf, [basename,'_intensity.fig'])
saveas(gcf, [basename,'_intensity.png'])

figure(2), hold on, 
plot(time,icell_av{1}, '-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
title('Average Intensity vs Time')
fig2pretty
xlim([-0.2 Inf])
saveas(gcf, [basename,'_intensityAvg.fig'])
saveas(gcf, [basename,'_intensityAvg.png'])

figure(3), hold on, 
for i=1:ncells
    plot(time, ratio{1}(i,:))
end
xlabel('Time (s)')
ylabel('Cellular Intensity/Background Intensity (A.U.)')
title('Intensity Ratio vs Time')
xlim([-0.2 Inf])
%ylim([0 1])
fig2pretty
saveas(gcf, [basename,'_ratio.fig'])
saveas(gcf, [basename,'_ratio.png'])

figure(4), hold on, 
plot(time,ratio_av{1}, '-r')
xlabel('Time (s)')
ylabel('Cellular Intensity/Background Intensity (A.U.)')
title('Average Intensity Ratio vs Time')
fig2pretty
xlim([-0.2 Inf])
%ylim([0 1])
saveas(gcf, [basename,'_ratioAvg.fig'])
saveas(gcf, [basename,'_ratioAvg.png'])

%phase1=[14:92];
phase2=[187:271];
phase3=[351:433];
phase4=[512:594];

f=cell(ncells,3);
for n=1:ncells
    yData=icell{1}(n, phase2);
    [xData, yData]=prepareCurveData(time(phase2), yData);
    f{n,1}=fit(xData, yData, 'exp1');
    plot(f{n,1}, xData, yData), pause, close
end

%% Functions
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