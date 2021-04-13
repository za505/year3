%PlasmolysisTrack.m
%Tracks and quantifies size of plasmolysis bays in bacteria
%To be used after imagealign.m, Eraseimagepart.m, & BacTrack.m
%Zarina Akbary, last updated 04/13/2021
clear 
close all
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='03082021_Exp3_colony3';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/PlasmolysisTrack_test/' basename '_647/'  basename '_erased'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/PlasmolysisTrack_test/' basename '_647'];%Directory to save the output .mat file to.
channel='647';
lscale=0.08;%%Microns per pixel.
tscale=10;%Frame rate.
thresh=0;%For default, enter zero.
IntThresh=10280;%Threshold used to enhance contrast. Default:35000
dr=1;%Radius of dilation before watershed %default: 1
sm=2;%Parameter used in edge detection %default: 2
minL=2;%Minimum cell length default: 2
maxL=40; %Maximum cell length
minW=0.2;%Minimum cell width default: 0.2
maxW=2.5;%Maximum cell width
minA=50;%Minimum cell area. default: 50
dotA_min=8800; %area range of the 'dots' (pillars) in the trap
dotA_max=8850;
cellLink=4;%Number of frames to ignore missing cells when tracking frame to frame
recrunch=0;%Display data from previously crunched data? 0=No, 1=Yes.
vis=1;%Display cell tracking? 0=No, 1=Yes.
checkhist=0;%Display image histogram? 0=No, 1=Yes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==1
    load([savedir '/' basename '_' channel '_BTphase.mat'])
else

%Determine number of frames
curdir=cd;
cd(dirname);
directory=dir('*.tif');
T=1; %length(directory);

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
    
    tiledlayout(5,1)
    nexttile
    imshow(im)

    %De-speckle image
    im=medfilt2(im);
     
    %Normalize images
    ppix=0.5;
    im=norm16bit(im,ppix);
    nexttile
    imshow(im)

    %Enhance contrast
    imc=imcomplement(im);
    nexttile
    imshow(imc)
    
    if checkhist==1 
        figure,imhist(imc);
    end
    
    if thresh==0;
        [imcounts,bins]=imhist(imc);
        [imcounts,idx]=sort(imcounts);
        bins=bins(idx);
        thresh1=bins(end-1);
    else
        thresh1=thresh;
    end
    thresh1
    imc=imadjust(imc,[thresh1/65535 1],[]);
    
    nexttile
    imshow(imc)
    
    [~,imcBins]=imhist(imc);
    %bSum=sum(bins);
    nBins=height(imcBins);
    binx=round(nBins*(0.03))
    IntThresh=imcBins(binx)
    
    %Find edges
    [ed2,thresh2]=edge(imc,'canny',[],sm*sqrt(2));
    thresh2
    
   
    
    %Clean image
    cc=bwconncomp(ed2,8);
    stats=regionprops(cc,imc,'Area','MeanIntensity');
    idx=find([stats.Area]>minA&[stats.Area]<1e5&[stats.MeanIntensity]>IntThresh);
    ed2=ismember(labelmatrix(cc),idx);
    
    nexttile
    imshow(ed2)
    
%     
%     %Close gaps in edges
%     despurred=bwmorph(ed2,'spur');
%     spurs=ed2-despurred;
%     [spy,spx]=find(spurs);
%     for k=1:length(spx)
%         ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)=ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)+rot90(ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1),2);
%         ed2(spy(k),spx(k))=1;
%     end
%     ed2=bwmorph(ed2,'bridge'); 
%   
%     se=strel('disk',dr);
%     ed2=imdilate(ed2,se);
%     ed2=bwmorph(ed2,'thin',2);
%     
%     %imshow(ed2)
%     %pause
%     
%     %Identify cells based on size and intensity
%     ed3=~ed2;
%     ed3(1,:)=ones(1,imN);
%     ed3(end,:)=ones(1,imN);
%     ed3(:,1)=ones(imM,1);
%     ed3(:,end)=ones(imM,1);
%     
%     cc=bwconncomp(ed3,4);
%     stats=regionprops(cc,imc,'Area','MeanIntensity');
%     idx=find([stats.Area]>minA&[stats.Area]<1e5&[stats.MeanIntensity]>3e4);
%     ed4=ismember(labelmatrix(cc),idx);
%     
%     %imshow(ed4)
%     %pause
%     
%     %Find cell areas and centroids
%     bw=bwmorph(ed4,'thicken');
%     [P,bw]=bwboundaries(bw,4,'noholes');
%     stats=regionprops(bw,'Area','Centroid','PixelIdxList');
%     
%     L=bwlabel(bw);    
%     labels(:,:,t)=L;
%     labels2(:,:,t)=bw;
%     
%     nc(t)=length(P);
%     areas=[stats.Area];
%     cents=cat(1,stats.Centroid);
%     a(nc(t),t)=0;
%     a(1:nc(t),t)=[stats.Area]';
%     centroids=cents;
%      
%     %Calculate smooth cell contours
%     for n=1:nc(t)
%          
%         rP=[P{n}(:,2),P{n}(:,1)];
%         px=[rP(1:end-1,1);rP(1:end-1,1);rP(:,1)];
%         py=[rP(1:end-1,2);rP(1:end-1,2);rP(:,2)];
%         sp=length(rP);
%         dS=sqrt(diff(px).^2+diff(py).^2);
%         S=[0 cumsum(dS)'];
% 
%         px=csaps(S,px,0.05,S);
%         py=csaps(S,py,0.05,S);
%         
%         px=px(sp+1:2*sp);
%         py=py(sp+1:2*sp);
%         
%         px(end)=px(1);
%         py(end)=py(1);
%         
%         dS=sqrt(diff(px).^2+diff(py).^2);
%         S=[0 cumsum(dS)];
%         ls=length(S);
%         DS(n,t)=S(end)/(ls-1);
%         Sn=(0:DS(n,t):S(end));
%         nx=spline(S,px,Sn);
%         ny=spline(S,py,Sn);
%         
%         boun{n,t}=[nx',ny'];
%         pxls{n,t}=stats(n).PixelIdxList;
%         
%     end
%     allcentroids=[allcentroids;centroids];
%     tstamp=[tstamp;ones(nc(t),1)*t];
%     cellnum=[cellnum;(1:nc(t))'];
% 
% if vis==1 & t >= T-10 | t <= 6
%    figure
%    imshow(im)
%    hold on
%    for k=1:nc(t)
%        plot(boun{k,t}(:,1),boun{k,t}(:,2),'-r')
%    end
%     
%   pause(0.3)
%   close all
%end
    toc

end
end


