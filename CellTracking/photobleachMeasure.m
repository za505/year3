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
basename='05182021_FITCK_010';
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05182021_analysis/05182021_control/' basename];
savename=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05182021_analysis/05182021_control/05182021_figures'];
recrunch=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curdir=cd;

cd(dirname);
directory=dir('*.tif');
T=5; %length(directory);
path(dirname,path)

%LB perfusion
imagename=directory(1).name;
im1=imread(imagename);
[imM,imN]=size(im1);

ppix=0.5;
im2=norm16bit(im1,ppix);

figure,imshow(im1,stretchlim(im1)*65000)
set(gcf,'Pointer','fullcross')
hold on
axis manual

count=0;
k=0;
while k~=1
    count=count+1;
    k=waitforbuttonpress;
    point1=get(gca,'CurrentPoint');   
    finalRect=rbbox;                   
    %pause
    point2=get(gca,'CurrentPoint');    
    point1=point1(1,1:2);              
    point2=point2(1,1:2);
    point1(point1<1)=1;
    point2(point2<1)=1;
    if point1(2)>imM
        point1(2)=imM;
    end
    if point1(1)>imN
        point1(1)=imN;
    end
    if point2(2)>imM
        point2(2)=imM;
    end
    if point2(1)>imN
        point2(1)=imN;
    end
    p1=min(point1,point2);%Calculate locations
    p2=max(point1,point2);
    offset = abs(point1-point2);%And dimensions
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    plot(x,y)
    
    rp1(count,:)=round(p1);
    rp2(count,:)=round(p2);
end

for t=1:T
    t
    
    %Load image
    imagename=directory(t).name;

    im=imread(imagename);
    [imM,imN]=size(im);
    
    
    figure
    imshow(im, [])
    hold on
        
    for n=1:count-1 %the last coordinates are a repeat
        
        w=rp2(n,1)-rp1(n,1);
        h=rp2(n,2)-rp1(n,2);
        rec=[rp1(n,1), rp1(n,2), w, h];

        rectangle('Position', rec, 'Edgecolor', 'r')
        
    end
    
    pause
end