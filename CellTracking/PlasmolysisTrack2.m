%PlasmolysisTrack2.m
%Tracks bacterial growth from phase image stacks.  
%Customized for B. subtilis.
%Bactrack2.m is customized for B. subtilis filaments.


%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save phase image stack in a directory
%by itself.  Also save the micromanager meGFPta file as 'basename.txt' in
%the matlab path.
%
%INPUT:
%basename: name of the image stack.
%dirname:the full pathname of the directory where you saved the image
%        stack.
%metaname(optional):full or relative pathname of micromanager meGFPta file from
%which to extract time points.  If it is relative path name, the
%directory in which it is saved must be on the matlab path.
%lscale: microscope calibration in microns per pixels.
%sm: width of the Gaussian filter used in edge finder equals sm*sqrt(2).
%minL: minimum length of cells;
%minW: minimum width of cells;
%maxW: maximum width of cells;
%recrunch:0 or 1.  if you've already tracked the data set and just want to
%         re-plot the data enter 1.
%
%OUTPUT:
%T: number of time points.
%time: vector of length T with time points.
%tmid: vector of length T-1 with interstitial time points.
%ncells: number of individual cells tracked.
%lcell: ncells x T matrix of cell lengths.
%wcell: ncells x T matrix of cell widths.
%acell: ncells x T matrix of cell areas
%ew: ncells x T matrix of circumferential strains.
%acell: ncells x T matrix of cell areas.
%v: ncells x T-1 matrix of cell strain rates.
%B: ncells x T cell array with cell contours.
%mlines: ncells x T cell array with cell midlines
%wav: vector of length T with average cell widths.
%wstd: vector of length T with standard deviations of cell widths.
%wste: vector of length T with standard error of cell widths.
%vav: vector of length T-1 with average cell strain rate.
%vstd: vector of length T-1 with standard deviation of strain rates.
%vste: vector of length T-1 with standard error of strain rates.
%avav: vector of length T-1 with average cell areal strain rate.
%avstd: vector of length T-1 with standard deviation of areal strain rates.
%avste: vector of length T-1 with standard error of areal strain rates.
%ndp: vecotr of lenth T-1 with number of data points averaged over.

%Calls on the following m-files:
%norm16bit.m
%polefinder.m
%cellcurvature.m
%meGFPta.m
%extrema.m
%EffectiveLength.m
%fig2pretty.m
%movingaverage.m

clear
close all

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='01212021_Exp1_colony1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/01212021_analysis/' basename '/' basename '_647/' basename '_erased'];%Directory that the image stack is saved in.
phasename=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/01212021_analysis/' basename '/' basename '_phase/' basename '_single'];
%cytoname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/01212021_analysis/' basename '/01212021_plasmolysisTrack2/01212021_GFP'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/01212021_analysis/' basename '/' basename '_figures'];%Directory to save the output .mat file to.
%metaname=['/Users/Rico/Documents/MATLAB/Matlab Ready/' basename '/meGFPta.txt'];%Name of meGFPta file.  Will only work if images were taken with micromanager.
lscale=0.08;%%Microns per pixel.
multiScale=0;
tscale=10;%Frame rate.
% tscale2=1;
% tpt1=120; %number of seconds passed by first time set
% tpt2=240; %number of seconds passed by second time set
% tpt3=480; %number of seconds passed by third time set
% tpt4=1320; %number of seconds passed by fourth time step
thresh=0;%For default, enter zero.
IntThresh=3000;%Threshold used to enhance contrast. Default:35000
dr=1;%Radius of dilation before watershed 
sm=2;%Parameter used in edge detection
minL=2;%Minimum cell length
minW=0.2;%Minimum cell width
maxW=1.5;%Maximum cell width
minA=50;%Minimum cell area. default 50
cellLink=4;%Number of frames to ignore missing cells when tracking frame to frame
recrunch=0;%Display data from previously crunched data? 0=No, 1=Yes.
vis=1;%Display cell tracking? 0=No, 1=Yes.
checkhist=0;%Display image histogram? 0=No, 1=Yes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savedir)
    load([basename '_PTcy5'])
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
    else
        thresh1=thresh;
    end
    imc=imadjust(imc,[thresh1/65535 1],[]);   
     
    %Find edges
    [ed2,thresh2]=edge(imc,'canny',[],sm*sqrt(2));
    
    %Clean image
    cc=bwconncomp(ed2,8);
    stats=regionprops(cc,imc,'Area','MeanIntensity');
    idx=find([stats.Area]>minA&[stats.Area]<1e5&[stats.MeanIntensity]>IntThresh);
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
    idx=find([stats.Area]>minA&[stats.Area]<1e5&[stats.MeanIntensity]>3e4);
    ed4=ismember(labelmatrix(cc),idx);
    
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
    
if vis==1 %& t >= T-10 | t <= 6
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


%Extract timepoints from metadata if it exists
if exist('metaname')==1
    if exist(metaname)==2
        %Extract timepoints from metadata
        tpoints=metadata(metaname);
        
        %Fix bug where micromanager screws up its timing
        dtime=diff(tpoints(1,:));
        fdt=find(dtime>2*(dtime(1)));
        if isempty(fdt)~=1
            fdt=fdt(1);
            tpoints(:,fdt+1:end)=tpoints(:,fdt+1:end)-tpoints(1,fdt+1)+tpoints(1,fdt)+(tpoints(1,fdt)-tpoints(1,fdt-1));
        end
    else
        tpoints=[0:T-1]*tscale;
    end
else
      if multiScale==0
     tpoints=[0:T-1]*tscale;
    elseif multiScale==1
     tpoint1=[0:tscale:tpt1];
     tpoint2=[tpt1+tscale2:tscale2:tpt2];
     tpoint3=[tpt2+tscale:tscale:tpt3];
     tlength=length(tpoint1)+length(tpoint2)+length(tpoint3);
     tpoint4=[tpt3+tscale2:tscale2:tpt4];
     tpoints=[tpoint1, tpoint2, tpoint3, tpoint4];
     tpoints=tpoints(1:T);
    end
end

time=tpoints(1,:);
time2=tpoints(end,:);

%Fix bug where micromanager screws up its timing
dtime=diff(time);
fdt=find(dtime>2*(dtime(1)));
if isempty(fdt)~=1
    fdt=fdt(1);
    time(fdt+1:end)=time(fdt+1:end)-dtime(fdt)+dtime(fdt-1);
    time2(fdt+1:end)=time2(fdt+1:end)-dtime(fdt)+dtime(fdt-1);
end

%Track cells frame to frame
tracks=zeros(size(im));
rcents=round(allcentroids);
linind=sub2ind(size(im),rcents(:,2),rcents(:,1));
tracks(linind)=1;

nhood=[0,1,0;1,1,1;0,1,0];
tracks=imdilate(tracks,strel('disk',cellLink));
overlay1=imoverlay(im,tracks,[.3 1 .3]);

[tracksL,ncells]=bwlabel(tracks,8);

lcell=zeros(ncells,T);
wcell=zeros(ncells,T);
acell=zeros(ncells,T);
pcell=zeros(ncells,T);
B=cell(ncells,T);
pixels=cell(ncells,T);
mlines=cell(ncells,T);
lcents=length(allcentroids);

for i=1:lcents
    cellid=tracksL(linind(i));
    lcell(cellid,tstamp(i))=l(cellnum(i),tstamp(i));
    wcell(cellid,tstamp(i))=w(cellnum(i),tstamp(i));
    acell(cellid,tstamp(i))=a(cellnum(i),tstamp(i));
    B{cellid,tstamp(i)}=boun{cellnum(i),tstamp(i)};
    pixels{cellid,tstamp(i)}=pxls{cellnum(i),tstamp(i)};
    mlines{cellid,tstamp(i)}=mline{cellnum(i),tstamp(i)};
    pcell(cellid,tstamp(i))=pole(cellnum(i),tstamp(i));
end

%Throw away cells with only one or two time points
delind=[];
for i=1:ncells
    if length(nonzeros(lcell(i,:)))<=2
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

%Calculate circumferential strain
wcell(isnan(wcell))=0;
ew=zeros(size(lcell));
for i=1:ncells
    ew(i,:)=wcell(i,:)/mean(wcell(i,:));
end

%Calculate the growth rate
deltat=time(2:end)-time(1:end-1);
v=(lcell(:,2:end)-lcell(:,1:end-1))./((lcell(:,1:end-1)+lcell(:,2:end))/2);
av=(acell(:,2:end)-acell(:,1:end-1))./((acell(:,1:end-1)+acell(:,2:end))/2);
for i=1:ncells
    v(i,:)=v(i,:)./deltat;
    av(i,:)=av(i,:)./deltat;
end

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
v(isnan(v))=0;
av(isnan(av))=0;
ew(isnan(ew))=0;

for t=1:T-1
    vav(t)=mean(nonzeros(v(:,t)));
    vstd(t)=std(nonzeros(v(:,t)));
end
vavm=ones(ncells,1)*vav;
vstdm=ones(ncells,1)*vstd;

inddel=abs(v-vavm)>2*vstdm&vstdm~=0;

v(inddel)=0;
av(inddel)=0;
lcell(inddel)=0;
acell(inddel)=0;
wcell(inddel)=0;
ew(inddel)=0;

for t=1:T
    wav(t)=mean(nonzeros(wcell(:,t)));
    wstd(t)=std(nonzeros(wcell(:,t)));
    wste(t)=wstd(t)./length(nonzeros(wcell(:,t)));
    ewav(t)=mean(nonzeros(ew(:,t)));
    ewstd(t)=std(nonzeros(ew(:,t)));
    ewste(t)=ewstd(t)./length(nonzeros(ew(:,t)));
end

for t=1:T-1
    vav(t)=mean(nonzeros(v(:,t)));
    vstd(t)=std(nonzeros(v(:,t)));
    ndp(t)=length(nonzeros(v(:,t)));
    vste(t)=vstd(t)/sqrt(ndp(t));
    avav(t)=mean(nonzeros(av(:,t)));
    avstd(t)=std(nonzeros(av(:,t)));
    avste(t)=avstd(t)/ndp(t);
end

v(v==0)=NaN;
av(av==0)=NaN;
lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;
ew(ew==0)=NaN;

tmid=(time(2:end)+time(1:end-1))/2;

cd(savedir);
save([basename '_PTlab'],'labels','labels2','-v7.3')

clear labels
clear labels2
save([basename '_PTcy5'])

end

%% let's load the phase images (note that the frames in phase should correspond to the frames in TADA)
%Determine number of frames
curdir=cd;
cd(phasename);
phasedir=dir('*.tif');

% pim=nan(imM,imN,T);
pim=nan(imM,imN);

y=1:imM; %rows are the y direction
x=1:imN; %columns are the x direction
[X,Y] = meshgrid(x,y);
total_pix=nan(ncells,T);
percent_plas=cell(ncells,T);
A_xv=cell(ncells,T);
A_yv=cell(ncells,T);
A_in=cell(ncells,T);
A_on=cell(ncells,T);
plasm=nan(ncells, T);
percent_plasmolysis=nan(ncells,T);
phase=cell(1,T);

for t=1:T
    t
    
    %Load phase image
    imagename=phasedir(t).name;
    imp=imread(imagename);
    phase{1,t}=imp;
    
    %overlay TADA boundaries
    figure
    imshow(imp, [])
    hold on
    for n=1:ncells
        plot(boun{n,t}(:,1),boun{n,t}(:,2),'-r')
    end
    pause, close
    
    %threshold so higher pixel values are set to 45000
    ptThresh=13000;
    
    for j=1:imN
        for i=1:imM
            if imp(i,j)>ptThresh
                pim(i,j)=45000;
            elseif imp(i,j)<ptThresh
                pim(i,j)=imp(i,j);
            else
                continue
            end
        end
    end
  
    pim=uint16(pim);
    
    figure
    imshow(pim)
    hold on
    for n=1:ncells
        plot(boun{n,t}(:,1),boun{n,t}(:,2),'-r')
    end
    pause, close
    
    for n=1:ncells
        A_xv{n,t}=boun{n,t}(:,1);%total boun
        A_yv{n,t}=boun{n,t}(:,2);%total boun
        [A_in{n,t},A_on{n,t}] = inpolygon(X,Y,A_xv{n,t},A_yv{n,t});%total boun
        total_pix(n,t)=sum(sum(A_in{n,t}));%total number of pixels within boundary
        percent_plas{n,t}=double(pim)-A_in{n,t};
        
        for j=1:imN %x direction
            for i=1:imM %y direction
                if A_in{n,t}(i,j)==0
                    percent_plas{n,t}(i,j)=NaN;
                end
            end
        end
        
        for j=1:imN %x direction
            for i=1:imM %y direction
                if percent_plas{n,t}(i,j)>=44999
                    percent_plas{n,t}(i,j)=NaN;
                elseif percent_plas{n,t}(i,j)>-2
                    percent_plas{n,t}(i,j)=1;   
                else
                    continue
                end
            end
        end
        
        plasm(n,t)=sum(nansum(percent_plas{n,t}(:,:)));
        percent_plasmolysis(n,t)=100-(plasm(n,t)/total_pix(n,t))*100;
    end
    
    for n=1:ncells
        for j=1:imN %x direction
            for i=1:imM %y direction
                if A_in{n,t}(i,j)==0
                    phase{1,t}(i,j)=45000;
                end
            end
        end
    end
 
    imshow(phase{1,t}), pause
end

cd(savedir)
% v = VideoWriter('pt2','MPEG-4');
% open(v);

for t=1:T
    intensity{1,t}=ones(size(boun{1,t},1),2);
    intensity{1,t}=intensity{1,t}*65530;
end

figure
hold on
for t=1:T
    plot3(boun{1,t}(:,1),boun{1,t}(:,2),intensity{1,t})
    hold on
    t=surf(X,Y,65535-phase{1,t}/1000)
    t.EdgeColor = 'interp';
    t.FaceColor = 'interp';
    view(0,95)
%     frame = getframe(gcf);
%     writeVideo(v,frame);
    pause(2)
    clf
end

% close(v)
close all

% %% let's load the GFP images (note that the frames in phase should correspond to the frames in TADA)
% %Determine number of frames
% curdir=cd;
% cd(cytoname);
% cytodir=dir('*.tif');
% 
% % pim=nan(imM,imN,T);
% cim=nan(imM,imN);
% 
% 
% 
% percent_plas_cyto=cell(ncells,T);
% plasm_cyto=nan(ncells, T);
% percent_plasmolysis_cyto=nan(ncells,T);
% cyto=cell(1,T);
% 
% for t=1:T
%     t
%     
%     %Load phase image
%     imagename=cytodir(t).name;
%     img=imread(imagename);
%     cyto{1,t}=img;
%     
%     %overlay TADA boundaries
%     figure
%     imshow(img, [])
%     hold on
%     for n=1:ncells
%         plot(boun{n,t}(:,1),boun{n,t}(:,2),'-r')
%     end
%     pause, close
%     
%     %threshold so higher pixel values are set to 45000
%     gfpThresh=1000;
%     
%     for j=1:imN
%         for i=1:imM
%             if img(i,j)<gfpThresh
%                 cim(i,j)=0;
%             elseif img(i,j)>gfpThresh
%                 cim(i,j)=img(i,j);
%             else
%                 continue
%             end
%         end
%     end
%   
%     cim=uint16(cim);
%     
%     figure
%     imshow(cim, [])
%     hold on
%     for n=1:ncells
%         plot(boun{n,t}(:,1),boun{n,t}(:,2),'-r')
%     end
%     pause, close
%     
%     for n=1:ncells
%         
%         percent_plas_cyto{n,t}=double(cim)-A_in{n,t};
%         
%         for j=1:imN %x direction
%             for i=1:imM %y direction
%                 if A_in{n,t}(i,j)==0
%                     percent_plas_cyto{n,t}(i,j)=NaN;
%                 end
%             end
%         end
%         
%         for j=1:imN %x direction
%             for i=1:imM %y direction
%                 if percent_plas_cyto{n,t}(i,j)<0
%                     percent_plas_cyto{n,t}(i,j)=NaN;
%                 elseif percent_plas_cyto{n,t}(i,j)>0
%                     percent_plas_cyto{n,t}(i,j)=1;   
%                 else
%                     continue
%                 end
%             end
%         end
%         
%         plasm_cyto(n,t)=sum(nansum(percent_plas_cyto{n,t}(:,:)));
%         percent_plasmolysis_cyto(n,t)=100-(plasm_cyto(n,t)/total_pix(n,t))*100;
%     end
% 
%     for n=1:ncells
%         for j=1:imN %x direction
%             for i=1:imM %y direction
%                 if A_in{n,t}(i,j)==0
%                     cyto{1,t}(i,j)=0;
%                 end
%             end
%         end
%     end
%     
%     imshow(cyto{1,t}, []), pause
% end
% 
% cd(savedir)
% % v = VideoWriter('pt2_cyto','MPEG-4');
% % open(v);
% 
% for t=1:T
%     intensity2{1,t}=ones(size(boun{1,t},1),2);
%     intensity2{1,t}=intensity{1,t}*2000;
% end
% 
% figure
% hold on
% for t=1:T
%     plot3(boun{1,t}(:,1),boun{1,t}(:,2),intensity2{1,t})
%     hold on
%     t=surf(X,Y,cyto{1,t})
%     rotate3d on;
%     t.EdgeColor = 'interp';
%     t.FaceColor = 'interp';
%     %view(0,95)
% %     frame = getframe(gcf);
% %     writeVideo(v,frame);
%     pause
%     clf
% end
% 
% close(v)
% close all