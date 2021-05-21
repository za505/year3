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
tscale=10;
vis=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savename)
    load([basename '_pbC'])
else
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

exclude=[];
for n=1:count-1
    for row=rp1(n,2):rp2(n,2)
        exc=[];
        for col=rp1(n,1):rp2(n,1)
            omit=[col, row];
            exc=[exc; omit];
        end
        exclude=[exclude; exc];
    end
end

coordinates=[];
for row=1:imM
    for col=1:imN
        tally=0;
        
        for n=1:height(exclude)
            if row~=exclude(n,2) & col~=exclude(n,1)
                tally=tally+1;
            end
        end
        
        if tally==height(exclude)
            coor=[row col];
            coordinates=[coordinates; coor];
        end
    end
end

pixelIntensity=nan(height(coordinates), T);
for t=1:T

    t

    %load the image
    imagename=directory(t).name;
    im=imread(imagename);
    
    %measure background intensity
    for n=1:height(coordinates)
         pixelIntensity(n,t)=im(coordinates(n,1),coordinates(n,2));
    end
    
     %check that the coordinates are correct
    if vis==1 & (t<=6 | t>=T-6)
        figure
        imshow(im, [])
        hold on
        %plot(rp1(1,2),rp1(1,1), 'r+', 'MarkerSize', 30, 'LineWidth', 2);
        for n=1:height(coordinates)
           plot(coordinates(n,2),coordinates(n,1),'-r') %x is im 'column', y is im 'row'
        end
        pause
        close
    end

end
    
%Calculate time variable
tpoints=[0:T-1]*tscale;
time=tpoints(1,:);

%calculate average intensity
intensityAvg=mean(pixelIntensity,1);

%calculate tau
idx=find(intensityAvg == intensityAvg(6)/2);
tau=time(idx)/log(2);

end

cd(savename)
save([basename '_pbC'])

%plot figures
figure(1)
plot(time, intensityAvg)
title('Average Intensity vs Time')
xlabel('Time (s)')
ylabel('Average Intensity')

figure(2), hold on
for n=1:length(pixelIntensity)
    plot(time, pixelIntensity(n,:))
end
title('Location-Specific Intensity vs Time')
xlabel('Time (s)')
ylabel('Intensity')
