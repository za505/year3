%photobleachMeasure.m
%Zarina Akbary, updated 05/19/21
%Calculates the rate of photobleaching

clear, close all

%INSTRUCTIONS FOR USE:
%save fluorescent image stacks directories by themselves

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
basename='05272021_FITCK_e100';
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05272021_analysis/' basename '/' basename '_aligned'];
conname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05272021_analysis/' basename];
savename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05272021_analysis/05272021_FITCK_figures'];
recrunch=0;
tscale=1;
vis=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savename)
    load([basename '_pb'])
else
curdir=cd;

cd(dirname);
directory=dir('*.tif');
T=length(directory);
path(dirname,path)

%LB perfusion
imagename=directory(1).name;
im1=imread(imagename);
[imM,imN]=size(im1);

[p1, p2]=selectRegion(imagename); %row=y, col=x

intensityAvg=[];
for t=1:T

    t

    %load the image
    imagename=directory(t).name;
    im=imread(imagename);
    
    %measure background intensity
    intensityAvg(1,t)=mean(mean(im(p1(2):p2(2),p1(1):p2(1))));
 
    
     %check that the coordinates are correct
    if vis==1 & (t<=6 | t>=T-6)
        figure
        imshow(im, [])
        hold on
        plot(p1(1),p1(2), 'r+', 'MarkerSize', 30, 'LineWidth', 2);
        %plot(p1(2):p2(2),p1(1):p2(1), '-r')
%         for n=1:height(coordinates)
%            plot(coordinates(n,2),coordinates(n,1),'-r') %x is im 'column', y is im 'row'
%         end
        pause
        close
    end
end
    
%now measure the intensity of the control
cd(conname);
imname=dir('*tif');
imname=imname(1).name;
imc=imread(imname);
conavg=mean(mean(imc(p1(2):p2(2),p1(1):p2(1))));
    
%Calculate time variable
tpoints=[0:T-1]*tscale;
%idx=find(tpoints>=264 & tpoints<265); %for this experiment specifically, there is a blip ~4.5 minutes
time=tpoints(1,:);

%index average intensity for relevant times
intensityAvg=intensityAvg(1,:);

end

cd(savename)
save([basename '_pb'])

%plot figures
figure(1)
plot(time, intensityAvg)
title('Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Average Intensity')
yline(conavg, '--k', '*PBS control') 
fig2pretty
saveas(gcf, [basename,'_intensityAvg.fig'])
saveas(gcf, [basename,'_intensityAvg.png'])

%%%%%Functions
function [p1, p2]=selectRegion(imagename)
        
        %Load last image
        %imagename=fluo_directory{i}(t).name;
        im2=imread(imagename);

        %Determine Background
        figure,imshow(im2,[]), hold on, title('Select Region')
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
        close
end 