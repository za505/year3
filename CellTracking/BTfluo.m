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
basename='06162021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06162021_analysis/' basename '/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06162021_analysis/' basename '/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06162021_analysis/' basename '/' basename '_aligned']}; 
recrunch=0;
frameBg=4; %this is the frame that you'll pick the background area from
multiScale=1;
tscale1=[100 50 25 12.5 6 3 1 0.263];
T1={[0:6];[7:20];[21:46];[47:96];[97:198];[199:400];[401:1002];[1003:1244]};
tadd=[0 650 625 612.5 606 603 601 60.1];
t1=14;
refFrame=13;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==0;

%calculate time
tpoints=[];
if multiScale==0
     tpoints=[0:T-1]*tscale;
elseif multiScale==1
    for i=1:length(tscale1)
        if i==1
            t_temp=T1{i}*tscale1(i);
            tpoints=[tpoints t_temp];
        else
            t_temp=[t_temp(end)+tscale1(i):tscale1(i):t_temp(end)+tadd(i)+tscale1(i)]; 
            tpoints=[tpoints t_temp];
        end
    end
    
end

time=tpoints(1,:);
time2=tpoints(end,:);

curdir=cd;
for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

cd(dirname)
load([basename '_BTphase'], 'pixels', 'ncells', 'lcell')

T=length(fluo_directory{1})-refFrame;
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
        imagename=fluo_directory{i}(t+refFrame).name;
        im=imread(imagename);
        for j=1:ncells
            intensities_temp(j,t)=mean(im(pixels{j,refFrame-2}));    
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

%Plot data
figure, hold on, 
for i=1:ncells
    plot(time, icell{1}(i,:))
end
xlabel('Time (h)')
ylabel('Intensity (A.U.)')
xlim([-0.2 Inf])
fig2pretty
saveas(gcf, [basename,'_intensity.fig'])
saveas(gcf, [basename,'_intensity.png'])

figure, hold on, 
for i=1:ncells
    plot(time, ratio{1}(i,:))
end
xlabel('Time (h)')
ylabel('Cellular Intensity/Background Intensity (A.U.)')
xlim([-0.2 Inf])
fig2pretty
saveas(gcf, [basename,'_ratio.fig'])
saveas(gcf, [basename,'_ratio.png'])

figure, hold on, 
plot(time,icell_av{1}, '-r')
xlabel('Time (h)')
ylabel('Intensity (A.U.)')
fig2pretty
xlim([-0.2 Inf])
saveas(gcf, [basename,'_intensityAvg.fig'])
saveas(gcf, [basename,'_intensityAvg.png'])

figure, hold on, 
plot(time,ratio_av{1}, '-r')
xlabel('Time (h)')
ylabel('Cellular Intensity/Background Intensity (A.U.)')
fig2pretty
xlim([-0.2 Inf])
saveas(gcf, [basename,'_ratioAvg.fig'])
saveas(gcf, [basename,'_ratioAvg.png'])

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