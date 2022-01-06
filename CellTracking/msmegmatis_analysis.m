%msmegmatis_analysis.m
%This code tracks cells perfused with occlusion dye or stained with FDAA.
%The tracking is evalutating by comparing the overlap with the
%corresponding phase and GFP images. Then, the area of the plasmolysis bay
%is calculated.
%
%Instructions:
%1) Align the frame of interest for all channels w/imagealignRico.m
%2) Erase clumped cells w/Eraseimagepart.m

%Output Format:
%The data is saved as a .mat file as '<basename>_msmeg' in the 'colony'
%directory, where basename is a input string chosen by the user.  The main
%output variables are:
%
%T:      the total number of time points
%time:   1 x T vector containing the time points corresponding to each frame.
%ncells: total number of cells.
%lcell:  ncells x T matrix containing the length of the cells in the movie.  
%        each row is an individual
%        cell and each column is a different time point
%wcell:  ncells x T matrix containing the the widths of the cells in the movie.
%v:      ncells x T-1 matrix with the growth rate (strain rate) of the cells.
%vav:    1 x T-1 vector containing the average growth rate of the cells, 
%        averaged across cells as a
%        function of time.
%B:      ncells x T cell array containing the x and y coordinates of the 
%        contours for each cell.
%
%calls on:
%norm16bit.m
%polefinder.m
%cellcurvature.m
%metadata.m
%extrema.m
%fig2pretty.m
%ciplot.m

clear
close all

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='12162021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_647/' basename '_postShock'];%Directory that the image stack is saved in.
phasedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_phase/' basename '_aligned'];
cytodir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_GFP/' basename '_aligned'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_647/' basename '_figures'];%Directory to save the output .mat file to.
lscale=0.08;%%Microns per pixel.
multiScale=0;
tscale=30;
thresh=0;%For default, enter zero.
IntThresh=35000;%Threshold used to enhance contrast. Default:35000
dr=1;%Radius of dilation before watershed 
sm=2;%Parameter used in edge detection
minL=2;%Minimum cell length
minW=0.2;%Minimum cell width, default 0.2
maxW=1.5;%Maximum cell width
minA=50;%Minimum cell area. default 50
maxA=2000; %maximum cell area. default 2000
cellLink=4;%Number of frames to ignore missing cells when tracking frame to frame
recrunch=0;%Display data from previously crunched data? 0=No, 1=Yes.
vis=1;%Display cell tracking? 0=No, 1=Yes.
checkhist=0;%Display image histogram? 0=No, 1=Yes.
troubleshooting=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savedir)
    load([basename '_BTphase'])
else

%Determine number of frames
curdir=cd;
cd(dirname);
directory=dir('*.tif');
T=length(directory);

cd(curdir);
path(dirname,path)

nc=zeros(1,T);
allcentroids=[];
cellnum=[];
tstamp=[];

%Pre-allocate matrices
ewav=zeros(1,T);
ewstd=zeros(1,T);
ewste=zeros(1,T);
ewndp=zeros(1,T);
wav=zeros(1,T);
wstd=zeros(1,T);
wste=zeros(1,T);
vav=zeros(1,T-1);
vstd=zeros(1,T-1);
vste=zeros(1,T-1);
ndp=zeros(1,T-1);
avav=zeros(1,T-1);
avstd=zeros(1,T-1);
avste=zeros(1,T-1);
a=zeros(1,T);
w=zeros(1,T);
l=zeros(1,T);
DS=zeros(1,T);
boun=cell(1,T);
pole=zeros(1,T);
mline=cell(1,T);

%Load first image
imagename=directory(1).name;
im=imread(imagename);
[imM,imN]=size(im);
labels=zeros(imM,imN,T);
labels2=zeros(imM,imN,T);

%Track cells
for t=1:T
    t
    
    %Load image
    imagename=directory(t).name;

    im=imread(imagename);
    [imM,imN]=size(im);
    
    %De-speckle image
    im=medfilt2(im);
    
    %Normalize images
    ppix=0.5;
    im=norm16bit(im,ppix);
    
    %Enhance contrast
    imc=imcomplement(im);
    
    if checkhist==1;
        figure,imhist(imc),pause;
    end
    
    if thresh==0;
        [imcounts,bins]=imhist(imc);
        [imcounts,idx]=sort(imcounts);
        bins=bins(idx);
        thresh1=bins(end-1);
        if thresh1==65535
            thresh1=bins(end);
        end
    else
        thresh1=thresh;
    end
    imc=imadjust(imc,[thresh1/65535 1],[]);   

    %Find edges
    [ed2,thresh2]=edge(imc,'canny',[],sm*sqrt(2));
    %imshow(ed2),pause, close
    
    %Clean image
    cc=bwconncomp(ed2,8);
    stats=regionprops(cc,imc,'Area','MeanIntensity');
    idx=find([stats.Area]>minA&[stats.Area]<maxA&[stats.MeanIntensity]>IntThresh);
    %idx=find([stats.Area]>minA&[stats.MeanIntensity]>IntThresh);
    ed2=ismember(labelmatrix(cc),idx);
    
    %Close gaps in edges
    despurred=bwmorph(ed2,'spur');
    spurs=ed2-despurred;
    [spy,spx]=find(spurs);
    for k=1:length(spx)
        ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)=ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)+rot90(ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1),2);
        ed2(spy(k),spx(k))=1;
    end
    ed2=bwmorph(ed2,'bridge'); 
  
    se=strel('disk',dr);
    ed2=imdilate(ed2,se);
    ed2=bwmorph(ed2,'thin',2);
    
    %Identify cells based on size and intensity
    ed3=~ed2;
    ed3(1,:)=ones(1,imN);
    ed3(end,:)=ones(1,imN);
    ed3(:,1)=ones(imM,1);
    ed3(:,end)=ones(imM,1);
    
    cc=bwconncomp(ed3,4);
    stats=regionprops(cc,imc,'Area','MeanIntensity');

    %idx=find([stats.Area]>minA&[stats.Area]<maxA&[stats.MeanIntensity]>3e4);
    idx=find([stats.Area]>minA&[stats.Area]<maxA); %new edit, 12/24/2021
    %idx=find([stats.Area]>minA&[stats.MeanIntensity]>3e4);
    ed4=ismember(labelmatrix(cc),idx);
    %imshow(ed4), pause, close
        
    %Find cell areas and centroids
    bw=bwmorph(ed4,'thicken');
    [P,bw]=bwboundaries(bw,4,'noholes');
    stats=regionprops(bw,'Area','Centroid','PixelIdxList');
    
    L=bwlabel(bw);    
    labels(:,:,t)=L;
    labels2(:,:,t)=bw;
    
    nc(t)=length(P);
    areas=[stats.Area];
    cents=cat(1,stats.Centroid);
    a(nc(t),t)=0;
    a(1:nc(t),t)=[stats.Area]';
    centroids=cents;
     
    %Calculate smooth cell contours
    for n=1:nc(t)
         
        rP=[P{n}(:,2),P{n}(:,1)];
        px=[rP(1:end-1,1);rP(1:end-1,1);rP(:,1)];
        py=[rP(1:end-1,2);rP(1:end-1,2);rP(:,2)];
        sp=length(rP);
        dS=sqrt(diff(px).^2+diff(py).^2);
        S=[0 cumsum(dS)'];

        px=csaps(S,px,0.05,S);
        py=csaps(S,py,0.05,S);
        
        px=px(sp+1:2*sp);
        py=py(sp+1:2*sp);
        
        px(end)=px(1);
        py(end)=py(1);
        
        dS=sqrt(diff(px).^2+diff(py).^2);
        S=[0 cumsum(dS)];
        ls=length(S);
        DS(n,t)=S(end)/(ls-1);
        Sn=(0:DS(n,t):S(end));
        nx=spline(S,px,Sn);
        ny=spline(S,py,Sn);
        
        boun{n,t}=[nx',ny'];
        pxls{n,t}=stats(n).PixelIdxList;
        
    end
    allcentroids=[allcentroids;centroids];
    tstamp=[tstamp;ones(nc(t),1)*t];
    cellnum=[cellnum;(1:nc(t))'];
    
if vis==1 & (t >= T-10 | t <= 6)
   figure
   imshow(im)
   hold on
   for k=1:nc(t)
       plot(boun{k,t}(:,1),boun{k,t}(:,2),'-r')
   end
    
  pause
  close all
end
    toc

end

%Calculate cell length, width, etc.
for t=1:T
    t
    for n=1:nc(t)    
         X=boun{n,t}(:,1);
         Y=boun{n,t}(:,2);   
                       
         [sX,~]=size(X);
         
         %Find poles
         [X,Y,pole(n,t)]=polefinder(X,Y);
         
         %Create mesh
         npts=min(pole(n,t),sX-pole(n,t)+1);
         S=(0:DS(n,t):(sX-1)*DS(n,t));
         
         s1=(0:S(pole(n,t))/(npts-1):S(pole(n,t)));
         s2=(S(pole(n,t)):(S(end)-S(pole(n,t)))/(npts-1):S(end));
         xc1=spline(S(1:pole(n,t)),X(1:pole(n,t)),s1);
         xc2=spline(S(pole(n,t):end),X(pole(n,t):end),s2);
         yc1=spline(S(1:pole(n,t)),Y(1:pole(n,t)),s1);
         yc2=spline(S(pole(n,t):end),Y(pole(n,t):end),s2);
         xc2=fliplr(xc2);
         yc2=fliplr(yc2);
         
         %Calculate midline
         mline{n,t}=[(xc1+xc2)'/2,(yc1+yc2)'/2];
         dxy=diff(mline{n,t}).^2;
         dl=sqrt(dxy(:,1)+dxy(:,2));
         l(n,t)=sum(dl);
         
         %Calculate width
         ls=[0 cumsum(dl)'];
         [~,mpos1]=min(abs(ls/l(n,t)-0.25));
         [~,mpos2]=min(abs(ls/l(n,t)-0.75));
         
         widths=sqrt((xc1-xc2).^2+(yc1-yc2).^2);
         w(n,t)=(widths(mpos1)+widths(mpos2))/2;
          
    end
end

%Throw away cells with only one or two time points, or less time points than
%we want in general
%throw away cells that lyse pre-maturely
ncells=height(w); 
lcell=l;
wcell=w;
acell=a;
pcell=pole;
B=boun;
pixels=pxls;
mlines=mline;
[ncells,~]=size(lcell);

%remove post analysis !!!

delind=[];
for i=1:ncells
    %if length(nonzeros(lcell(i,:)))<=2|sum(cellfun(@isempty, B(i,:)))/T>0.1
    if sum(acell(i,:)<minA)>20%remove cells that area too small or big, edit 12/24/21
        delind=[delind;i];
    end
end

lcell(delind,:)=[];
wcell(delind,:)=[];
acell(delind,:)=[];
pcell(delind,:)=[];
B(delind,:)=[];
pixels(delind,:)=[];
mlines(delind,:)=[];
[ncells,~]=size(lcell);

lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;
pcell(pcell==0)=NaN;

%Dimsionalize the variables
lcell=lcell*lscale;
wcell=wcell*lscale;
acell=acell*lscale^2;

%Throw away cells that are too short or too fat or too skinny
lcell(lcell<minL|wcell>maxW|wcell<minW)=NaN;
wcell(lcell<minL|wcell>maxW|wcell<minW)=NaN;
acell(lcell<minL|wcell>maxW|wcell<minW)=NaN;

%Calculate total length of all cells
lcell(isnan(lcell))=0;
acell(isnan(acell))=0;
for t=1:T
    ltotal(t)=sum(nonzeros(l(:,t)));
    atotal(t)=sum(nonzeros(a(:,t)));
end
lcell(lcell==0)=NaN;
acell(acell==0)=NaN;

%Throw out outliers and calculate the average width,strain and strain rate 
%across cells
lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;


cd(savedir);
save([basename '_BTlab'],'labels','labels2','-v7.3')
clear labels
clear labels2
save([basename '_BT'])
end

%% Troubleshooting
if troubleshooting == 1

    %load phase image
    cd(phasedir);
    directory2=dir('*.tif');
    imagename=directory2(1).name;
    im2=imread(imagename);

    %load phase image
    cd(cytodir);
    directory3=dir('*.tif');
    imagename=directory3(1).name;
    im3=imread(imagename);

    
    for k=1:ncells
       k
       figure
       imshow(im)
       hold on

       for t=1:T
         if isempty(B{k,t})==0
            plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
         else
             continue
         end
       end
       
%        figure(2)
%        imshow(im2, [])
%        hold on
% 
%        for t=1:T
%          if isempty(B{k,t})==0
%             plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
%          else
%              continue
%          end
%        end
%        
%        figure(3)
%        imshow(im3, [])
%        hold on
% 
%        for t=1:T
%          if isempty(B{k,t})==0
%             plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
%          else
%              continue
%          end
%        end

      pause
      close all
    end
end


%% save data
cd(savedir);
save([basename '_colony1_postShock'])
