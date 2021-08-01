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
%fluordir: name of folder containing fluorescent image stacks to quantify.

%OUTPUT:
%icell: a matrix (ncellxT) with the fluorescent intensities of each
        %cell, where rows are the cells and columns are time points.
%icell_av: Each entry is a vector containing the population-
        %average of the single-cell fluorescent intensities.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='03292021_Exp3';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03292021_analysis/' basename '/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory that the BTphase.mat file is saved in
savedir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03292021_analysis/' basename '/' basename '_colony1/' basename '_AF/' basename '_figures'];%Directory to save the output .mat file to.
fluordir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/03292021_analysis/' basename '/' basename '_colony1/' basename '_AF/' basename '_aligned']; 
recrunch=1;
%frameAuto=[185:200]; %post lysis *PBS perfusion (to calculate autofluorescence)
frameBg=45;
troubleshoot=2;
fluorFrames=[38:97]; %frames where AF is perfused
%dotThresh=490; %this is the lightest pixel value of the dot during AF perfusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==0;
    
cd(fluordir); 
fluo_directory=dir('*.tif');

cd(dirname)
load([basename '_BTphase'], 'pixels', 'ncells', 'lcell', 'B', 'time')

T=length(fluo_directory);
icell=nan(ncells,T);
icell_av=nan(ncells,T);
ratio=nan(ncells,T);
ratio_av=nan(ncells,T);
bgIntensity=nan(1, T);
 
cd(fluordir); 
imagename=fluo_directory(frameBg).name;
[p1, p2]=getBackground(imagename);
    
    for t=1:T
        t
        imagename=fluo_directory(t).name;
        im=imread(imagename);
        for j=1:ncells
            icell(j,t)=mean(im(pixels{j,t}));    
        end
        
         %measure background level
         bglevel = measureBackground(imagename, p1, p2);
         bgIntensity(1,t)=bglevel; 
    end
    
%     intensities_temp(intensities_temp==0)=NaN;
%     icell{i}=intensities_temp;

    for j=1:ncells
        ratio(j,:)=icell(j,:)./bgIntensity(1,:);
    end

    icell_av=mean(icell, 1, 'omitnan');
    ratio_av=mean(ratio,1, 'omitnan');

cd(savedir)
save([basename '_BTfluo'])

elseif recrunch==1
    cd(savedir)
    load ([basename '_BTfluo'])
    troubleshoot=0;
end

%% Troubleshooting
cd(fluordir); 
 if troubleshoot==1
     for t=1:T
            t
            imagename=fluo_directory(t).name;
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
            imagename=fluo_directory(T).name;
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
    plot(time, icell(i,:))
end
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
title('Intensity vs Time')
xlim([-0.2 Inf])
fig2pretty
xline(250, '--', {'*PBS'})
xline(370, '--', {'*PBS + AF'})
saveas(gcf, [basename,'_intensity.fig'])
saveas(gcf, [basename,'_intensity.png'])

figure(2), hold on, 
plot(time,icell_av, '-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
title('Average Intensity vs Time')
fig2pretty
xlim([-0.2 Inf])
xline(250, '--', {'*PBS'})
xline(370, '--', {'*PBS + AF'})
saveas(gcf, [basename,'_intensityAvg.fig'])
saveas(gcf, [basename,'_intensityAvg.png'])

figure(3), hold on, 
for i=1:ncells
    plot(time, ratio(i,:))
end
xlabel('Time (s)')
ylabel('Cellular Intensity/Background Intensity (A.U.)')
title('Intensity Ratio vs Time')
xlim([-0.2 Inf])
%ylim([0 1])
xline(250, '--', {'*PBS'})
xline(370, '--', {'*PBS + AF'})
fig2pretty
saveas(gcf, [basename,'_ratio.fig'])
saveas(gcf, [basename,'_ratio.png'])

figure(4), hold on, 
plot(time,ratio_av, '-r')
xlabel('Time (s)')
ylabel('Cellular Intensity/Background Intensity (A.U.)')
title('Average Intensity Ratio vs Time')
fig2pretty
xlim([-0.2 Inf])
%ylim([0 1])
xline(250, '--', {'*PBS'})
xline(370, '--', {'*PBS + AF'})
saveas(gcf, [basename,'_ratioAvg.fig'])
saveas(gcf, [basename,'_ratioAvg.png'])

figure(5), hold on, 
plot(time,bgIntensity, '-r')
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
title('Background Intensity vs Time')
fig2pretty
xlim([-0.2 Inf])
xline(250, '--', {'*PBS'})
xline(370, '--', {'*PBS + AF'})
saveas(gcf, [basename,'_intensityBg.fig'])
saveas(gcf, [basename,'_intensityBg.png'])

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