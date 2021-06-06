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
basename='06062021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_aligned']; ['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_aligned']}; 
recrunch=0;
% frameSwitch=187; %this is the initial frame for the switch (after which we pulse w/ and w/out Mg2+)
% frameInitial=14; %this is the FIRST frame that has dye without lysing the membrane
% frameAuto=90; %this is the LAST frame that has dye without lysing the membrane
frameBg=85; %this is the frame that you'll pick the background area from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==0;

curdir=cd;
for i=1:length(channels)
    cd(channels{i}); 
    fluo_directory{i}=dir('*.tif');
end

cd(dirname)
load([basename '_BTphase'])
load([basename '_BTlab'])

intensities=cell(length(channels),1);
bgIntensity=nan(length(channels),T);
Cout=nan(length(channels),T);
Cadj=cell(length(channels),1);
bgAuto=nan(length(channels),1);
Cin{i}=cell(length(channels),1);
iDiff{i}=cell(length(channels),1);
ratio{i}=cell(length(channels),1);
    
for i=1:length(channels)
    cd(channels{i}); 
    intensities_temp=zeros(size(lcell));
    
    imagename=fluo_directory{i}(frameBg).name;
    [p1, p2]=getBackground(imagename);
    
    bgTemp=[];
    
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
        
%         %get auto fluorescence
%         if t>frameAuto+1 & t<frameSwitch-1 
%             bgTemp=[bgTemp bgIntensity(i,t)];
%         end
    
    end
    
    intensities_temp(intensities_temp==0)=NaN;
    icell{i}=intensities_temp;
   
    Cin{i}=icell{i};
    iDiff{i}=icell{i};
    ratio{i}=icell{i};
    
    %now take the mean of bgAuto
    %bgAuto(i)=mean(bgTemp,2); %take the temporal average
    
    %Subtract the autofluorescence from the raw intensity
    Cout(i,:)=bgIntensity(i,:); %-bgAuto(i);
    for n=1:ncells
        Cin{i}(n,:)=icell{i}(n,:); %-bgAuto;
        iDiff{i}(n,:)=Cin{i}(n,:)-Cout(i,:);
        ratio{i}(n,:)=Cin{i}(n,:)./Cout(i,:);
    end
end

icell_av=cell(length(channels),1);
for i=1:length(channels)
    icell_av{i}=nanmean(icell{i});
end

elseif recrunch==1
    cd(savedir)
    load ([basename '_BTfluo'])
end

%Plot data
cd(savedir)

figure, hold on, 
for i=1:ncells
    plot(time,icell{1}(i,:))
end
xlabel('Time (s)')
ylabel('FSS Intensity (A.U.)')
fig2pretty
saveas(gcf, [basename '_intensityFSS.fig'])
saveas(gcf, [basename '_intensityFSS.png'])

%Plot adj background fluorescence
figure, hold on 
title('Adjusted Background Intensity vs Time')
plot(time,Cout(1,:))
xlabel('Time (s)')
ylabel('FSS Intensity (A.U.)')
fig2pretty
saveas(gcf, [basename,'_CoutFSS.fig'])
saveas(gcf, [basename,'_CoutFSS.png'])

%Plot adj adjusted cellular fluorescence
figure, hold on, title('Phase Contrast vs Time')
for n=1:ncells
    plot(time,Cin{2}(n,:))
end
xlabel('Time (s)')
ylabel('Phase Contrast (A.U.)')
fig2pretty
saveas(gcf, [basename '_intensityPhase.fig'])
saveas(gcf, [basename '_intensityPhase.png'])

%Let's plot Cin-Cout
figure, hold on, title('Cin-Cout vs Time')
for n=1:ncells
    plot(time,iDiff{1}(n,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
saveas(gcf, [basename '_iDiffFSS.fig'])
saveas(gcf, [basename '_iDiffFSS.png'])

%Let's plot Cin/Cout
figure, hold on, title('Ratio vs Time')
for n=1:ncells
    plot(time,ratio{1}(n,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
fig2pretty
saveas(gcf, [basename '_ratioFSS.fig'])
saveas(gcf, [basename '_ratioFSS.png'])

figure, hold on, 
plot(time,icell_av{1},'-r')
xlabel('Time (s)')
ylabel('Average Intensity (A.U.)')
fig2pretty
saveas(gcf, [basename,'_intensityAvgFSS.fig'])
saveas(gcf, [basename,'_intensityAvgFSS.png'])

figure, hold on, 
plot(time,icell_av{2},'-r')
xlabel('Time (s)')
ylabel('Phase Contrast (A.U.)')
fig2pretty
saveas(gcf, [basename,'_intensityAvgPhase.fig'])
saveas(gcf, [basename,'_intensityAvgPhase.png'])

cd(savedir)
save([basename '_BTfluo'])

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