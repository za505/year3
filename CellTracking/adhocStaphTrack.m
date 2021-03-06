%StaphTrack.m
%Tracks S. aureus growth from fluorescent image stacks.  
%
%INSTRUCTIONS FOR USE:
%Save fluorescent image stack in a directory by itself.  If desired, also 
%save the micromanager metadata file.
%
%INPUT: 
%basename: assign a name to the run.
%dirname:the full pathname of the directory where you saved the image
%        stack.
%metaname(optional):full or relative pathname of micromanager metadata file from
    %which to extract time points.  If it is relative path name, the
    %directory in which it is saved must be on the matlab path.
%hsize:size of gaussian filter.
%sigma:width of gaussian filter.
%tscale: frame rate in seconds.
%lscale: microscope calibration in microns per pixels.
%sfac: number of time points over which to smooth length traces.
%recrunch:0 or 1.  if you've already tracked the data set and just want to
%         re-plot the data enter 1.
%
%OUTPUT:
%T: number of time points.
%nc: vector of length T containing number of cells at each time point.
%B:NxT cell array containing x and y coordinates of cell outlines vs time.
%vcell:NxT matrix containing volumes of cells at each time point.
%vtot:1xT vector containing the sum of the volumes of all cells.
%shellx,shelly,shellz: NxT cell array containing x,y,z coordinates of cell
    %surfaces at each time point.
%icell:NxT vector containing mean cell intensity vs. time.
%P:NxT cell array containing column vector of pixel values within cells.
%time: vector of length T with time points.

%Calls on the following m-files:
%norm16bit.m
%rotatexy.m
%cellmesh.m
%rotatemesh.m
%metadata.m

clear
close all

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='04292022_Exp1';
dirname=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/04292022_analysis/' basename '/' basename '_mCherry/' basename '_aligned'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/04292022_analysis/' basename '/' basename '_mCherry/' basename '_figures'];
%metaname=['/Volumes/HD-PNTU3/Matlab Ready/' basename '/' basename '.txt'];
hsize=4;
sigma=2;
sfac=1;
lscale=0.08;
tscale1=60;
tscale2=4;
tpoint1=[0:tscale1:7*60];
tpoint2=[tpoint1(end)+tscale2:tscale2:(256*4)+tpoint1(end)];
recrunch=0;
multiScale=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==1
    load([basename '_ST'])
else

%determine number of frames
cd(dirname);
directory=dir('*.tif');
T=length(directory);

nc=zeros(1,T);
allcentroids=[];
cellnum=[];
tstamp=[];

%define time
if multiScale==0
    time=[0:T-1]*tscale;
 elseif multiScale==1
    time=[tpoint1, tpoint2];
end

        
imagename=directory(1).name;
[p1, p2]=getBackground(imagename);

%pick cells to track
[rp1, rp2, count, im1]=pickCell(dirname, directory, 1);
icell_intensity = nan(count, T);
bg_intensity = nan(1, T);

%track cells
for t=1:T
    t
    
    %load image
    imagename=directory(t).name;  
    im=imread(imagename);
    imraw=im;
    [imM,imN]=size(im);
%     figure,imshow(im,[]), pause
%     message = 'loaded image'
    
    %get background
    bg_intensity(1, t)=measureBackground(im, p1, p2);
    
    %smooth image with gaussian filter
    H=fspecial('gaussian',hsize,sigma);
    im=imfilter(im,H,'replicate');
%     figure,imshow(im,[]), pause
%     message = 'filtered the image'
    
    im=imcomplement(im);
%     figure,imshow(im,[]),pause
%     message = 'complemented'
    
    %normalize image
    ppix=0.5;
    im=norm16bit(im,ppix);
    imraw=norm16bit(imraw,ppix);
%     figure,imshow(im,[]),pause
%     message = 'normalized the image'
    
    for n=1:count
        icell_intensity(n, t) = measureCell(im, rp1, rp2, n);
    end
    
end

adjintensity = icell_intensity - bg_intensity;

end

%% plot data
cd(savedir)

figure, title('Cell Intensity vs. Time')
hold on
for i=1:count
    plot(time./60,icell_intensity(i, :),'-m')
end
ylim([0 Inf])
xlabel('Time (minutes)')
ylabel('Intensity')
fig2pretty

figure, title('Adjusted Intensity vs. Time')
hold on
for i=1:count
    plot(time./60,adjintensity(i, :),'-m')
end
ylim([0 Inf])
xlabel('Time (minutes)')
ylabel('Adjusted Intensity')
fig2pretty
saveas(gcf, [basename '.fig'])
saveas(gcf, [basename '.png'])
%% save data
save([basename '_ST'])

%% Functions
function [rp1, rp2, count, im1]=pickCell(dirname, directory, t)
    
    %load image
    cd(dirname)
    imagename=directory(t).name;
    im1=imread(imagename);
    [imM,imN]=size(im1);

    ppix=0.5;
    im2=norm16bit(im1,ppix);

    %figure, imshowpair(imA, imB, 'montage')
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
 end


 function intensity = measureCell(im, rp1, rp2, n)

        values = im(rp1(n,2):rp2(n,2),rp1(n,1):rp2(n,1));
        intensity = mean(mean(values));
        
 end 
 
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

 function bglevel = measureBackground(im2, p1, p2)

         %Load last image
         %imagename=fluo_directory{i}(t).name;
         %im2=imread(imagename);

         %Determine background
         backim=im2(p1(2):p2(2),p1(1):p2(1));
         [counts,bins]=imhist(backim);
         [~,binnum]=max(counts);
         maxpos=bins(binnum);
         bglevel=mean(mean(backim));

 end 