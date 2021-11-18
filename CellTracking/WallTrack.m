%WallTrack.m
%This code tracks cell growth from stacks of fluorescent images where
%the cell periphery has been stained. it also records the intensity of a
%tracer dye taken in a different channel simultaneously
%
%Instructions:
%1)If there is significant drift in your movie, register the images before 
%  you use this code or else it may have trouble tracking cells frame to frame
%2)Save the image sequence in a folder by itself.
%3)Modify the User Input section below accordingly.
%
%Output Format:
%The data is saved as a .mat file as 'basename_WT' in the working
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
%cd /Users/dylanfitzmaurice/Documents/MATLAB/
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input

%File Parameters
basename='05262021_Exp1';%Name to save the file to
skp=[];%frames to skip
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/PlasmolysisTrack/' basename '_TADA/' basename '_erased'];%directory where the imagestack is saved
imdir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/PlasmolysisTrack/' basename '_GFP/' basename '_aligned'];%directory where the imagestack is saved
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/PlasmolysisTrack/' basename '_TADA/' basename '_figures'];
%metaname=['/Matlab Ready/' basename '/' basename '.txt'];%full pathname of the metadata file

%Tracking Parameters
ll=0.35;%parameter used to enhance contrast of image before performing edge detection % default=0.35
ppix=0.4;%percent of pixels to saturate during normalization default=0.6, %0.26 FOR GFP
cutL=0.5;%low threshold for edge detector default=0.7
cutH=1;%high threshold for edge detector default=1
dr=1;%radius of structure element to perform morphological operations

%Analysis Parameters
lscale=0.08;%Microns per pixel
multiscale=0; 
tscale=120;%Time between frames %Default=30
minA=1;%Minimum cell area %Default=20
maxI=70000;%Maximum cell intensity
minW=0.1;%Minimum cell width
maxW=3;%Maximum cell width 
minL=0.5;%Minimum cell length %Default=2
maxD=10;%Maximum cell width

%Other Parameters
recrunch=0;%whether or not to re-plot a file using new analysis parameters
plotcells=1;%whether or not to visualize the tracking
troubleshoot=1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==1
    cd(savedir);
    load([basename '_WT'],'-regexp','^(?!lscale$|tscale$|minA$|maxI$|minW$|maxW$|minL$|maxD$).')
else

%determine number of frames
curdir=cd;
cd(dirname);
directory=dir('*.tif');
T=length(directory);

cd(curdir);
path(dirname,path)

%pre-allocate matrices
allcentroids=[];
tstamp=[];
cellnum=[];
int=zeros(1,T);
nc=zeros(1,T);
l=zeros(1,T);
w=zeros(1,T);
boun=cell(1,T);

%     % Initialize video
% myVideo = VideoWriter('GFP_tracking'); %open video file
% myVideo.FrameRate = 0.002;  %
% open(myVideo)

for t=1:T
    t
    
    %load images
    imagename=directory(t).name;
    im=imread(imagename);
    imraw=im;
    %imraw=imcomplement(im); %
    [imm,imn]=size(im);
       
    %check image histogram, if necssary 
%     message=['raw image histogram']
%     figure,imhist(im),pause, close
    
    %normalize images
    [imcounts,imbin]=imhist(im);
    csic=cumsum(fliplr(imcounts'));
    pcsic=csic/sum(imcounts);
    [~,mpos]=min(abs(pcsic-ppix/100));
    maxintensity=imbin(end-mpos);
    im=imadjust(im,[0 maxintensity/65535],[]);

    message=['normalized image']
    figure,imshow(im),pause, close
    message=['normalized histogram']
    figure,imhist(im),pause, close
    %im=norm16bit(im,ppix);%color decode edit DF? 
        
    %enhance contrast
    im=imadjust(im,[ll 1],[]);
%     message=['adjusted image']
%     figure,imshow(im),pause,close
    
    %perform edge detection
    im=medfilt2(im);
    bw=edge(im,'canny',cutL,cutH);
    bw=imcomplement(bw);
%     message=['image edge detection']
%     figure,imshow(bw),pause, close
    
    %identify cells
    se=strel('disk',dr);
    bw=imerode(bw,se);
    bwL=bwlabel(bw,8);
    bw(bwL==1)=0;
%     message=['post erosion']
%     figure,imshow(bw),pause, close
    
    %throw away cells that are too small, too bright, or too close to the
    %border
    cc=bwconncomp(bw);
    stats=regionprops(cc, 'Area','Centroid');
    idx=find([stats.Area]>minA);
    cents=cat(1,stats.Centroid);
    centsx=cents(:,1);
    centsy=cents(:,2);
    idx2=find(centsx>10);
    idx3=find(centsx<imn-10);
    idx4=find(centsy>10);
    idx5=find(centsy<imm-10);
    
    stats=regionprops(cc,im,'MeanIntensity');
    idx6=find([stats.MeanIntensity]<maxI);
    
    idx=intersect(idx,idx2);
    idx=intersect(idx,idx3);
    idx=intersect(idx,idx4);
    idx=intersect(idx,idx5);
    idx=intersect(idx,idx6);
    
    bw2=ismember(labelmatrix(cc),idx);
    
    %watershed image
    im(bw2==1)=-Inf;
%     message=['before watershed']
%     imshow(im, []), pause, close
    ws=watershed(im);
%     message=['watershed image']
%     imshow(ws, []), pause, close
    
    %throw away background watersheds
    stats=regionprops(ws,'Area');
    idx=find([stats.Area]<12000 & [stats.Area]>minA);
    ws=ismember(ws, idx);
%     message=['without background watersheds?']
%     imshow(ws, []), pause, close
    [L,n]=bwlabel(ws);
    labels(:,:,t)=L;
    
    %calculate boundaries, length and width of cells
    stats=regionprops(L,'centroid','orientation','PixelIdxList');
    centroids=cat(1,stats.Centroid);
    orientations=cat(1,stats.Orientation);
    nc(t)=n;
    nhood=[0,1,0;1,1,1;0,1,0]; 
    P=cell(nc(t),1);
    l(nc(t),t)=0;
    w(nc(t),t)=0;
    boun{nc(t),t}=0;
    for i=1:nc(t)
        %calculate boundaries
        Lsub=imdilate(L==i,nhood);
        P(i)=bwboundaries(Lsub,8,'noholes');
        rP=[P{i}(:,2),P{i}(:,1)];
        
        %smooth boundaries
        px=[rP(1:end-1,1);rP(1:end-1,1);rP(:,1)];
        py=[rP(1:end-1,2);rP(1:end-1,2);rP(:,2)];
        sp=length(rP);
        dS=sqrt(diff(px).^2+diff(py).^2);
        S=[0 cumsum(dS)'];
        ls=length(S);
        DS(i,t)=S(end)/(ls-1);
        ww=ones(1,ls)/1;
        ww(end)=10;
        ww(1)=10;
        px=csaps(S,px,0.05,S,ww);
        py=csaps(S,py,0.05,S,ww);
        rPs=[px(sp:2*sp)',py(sp:2*sp)'];
        
        px=px(sp+1:2*sp);
        py=py(sp+1:2*sp);
        
        px(end)=px(1);
        py(end)=py(1);
        
        dS=sqrt(diff(px).^2+diff(py).^2);
        S=[0 cumsum(dS)];
        ls=length(S);
        DS(i,t)=S(end)/(ls-1);
        Sn=(0:DS(i,t):S(end));
        nx=spline(S,px,Sn);
        ny=spline(S,py,Sn);
         
        boun{i,t}=[nx',ny'];
        pxls{i,t}=stats(i).PixelIdxList;
        sX=length(nx);
        
        %find poles
        [X,Y,pole(i,t)]=polefinder(nx',ny');
        
        %create mesh
        npts=min(pole(i,t),sX-pole(i,t)+1);
        S=(0:DS(i,t):(sX-1)*DS(i,t));
        
        s1=(0:S(pole(i,t))/(npts-1):S(pole(i,t)));
        s2=(S(pole(i,t)):(S(end)-S(pole(i,t)))/(npts-1):S(end));
        xc1=spline(S(1:pole(i,t)),X(1:pole(i,t)),s1);
        xc2=spline(S(pole(i,t):end),X(pole(i,t):end),s2);
        yc1=spline(S(1:pole(i,t)),Y(1:pole(i,t)),s1);
        yc2=spline(S(pole(i,t):end),Y(pole(i,t):end),s2);
        xc2=fliplr(xc2);
        yc2=fliplr(yc2);
        
        %calculate midline
        mline{i,t}=[(xc1+xc2)'/2,(yc1+yc2)'/2];
        dxy=diff(mline{i,t}).^2;
        dl=sqrt(dxy(:,1)+dxy(:,2));
        l(i,t)=sum(dl);
        
        %calculate width
        ls=[0 cumsum(dl)'];
        [~,mpos1]=min(abs(ls/l(i,t)-0.25));
        [~,mpos2]=min(abs(ls/l(i,t)-0.75));
        
        widths=sqrt((xc1-xc2).^2+(yc1-yc2).^2);
        w(i,t)=(widths(mpos1)+widths(mpos2))/2;

    end

    allcentroids=[allcentroids;centroids];
    tstamp=[tstamp;ones(nc(t),1)*t];
    cellnum=[cellnum;(1:nc(t))'];
    
    if plotcells==1
        %overlay cell boundaries with image
        imshow(imraw,[])
        hold on
        for k=1:length(P)
            cellboun=P{k};
            plot(cellboun(:,2),cellboun(:,1),'-r')
            %hold on
            plot(centsx(:,1),centsy(:,1),'b*')
            %pause(0.1)
        end
        pause
%          M(t) = getframe;
%       pause(0.01) %Pause and grab frame
%       frame = getframe(gcf); %get frame
%         open(myVideo)
%       writeVideo(myVideo, frame);
% end
% close(myVideo)
        pause(0.1)
        %pause
        clf
    end
    
    toc
end

end%of recrunch

%extract timepoints from metadata
if exist('metaname')==1
    if exist(metaname)==2
        tpoints=metadata(metaname);
       
        %fix bug where micromanager screws up its timing
        dtime=diff(tpoints(1,:));
        fdt=find(dtime>2*(dtime(1)));
        if isempty(fdt)~=1
            fdt=fdt(1);
            tpoints(:,fdt+1:end)=tpoints(:,fdt+1:end)-tpoints(1,fdt+1)+tpoints(1,fdt)+(tpoints(1,fdt)-tpoints(1,fdt-1));
        end
        time=tpoints(2,:);
    else
        tpoints=[0:T-1]*tscale;
        time=tpoints(1,:);
    end
else
    if multiscale==0
        tpoints=[0:T-1]*tscale;
        time=tpoints(1,:);
    else
        time=[490 970 1300 1320];
    end
end

%track cells frame to frame
tracks=zeros(size(im));
rcents=round(allcentroids);
linind=sub2ind(size(im),rcents(:,2),rcents(:,1));
tracks(linind)=1;

nhood=[0,1,0;1,1,1;0,1,0];
tracks=imdilate(tracks,nhood);

[tracksL,ncells]=bwlabel(tracks,8);

lcell=zeros(ncells,T);
wcell=zeros(ncells,T);
B=cell(ncells,T);
pixels=cell(ncells,T);%edit for color decode DF
lcents=length(allcentroids);
for i=1:lcents
    cellid=tracksL(linind(i));
    lcell(cellid,tstamp(i))=l(cellnum(i),tstamp(i));
    wcell(cellid,tstamp(i))=w(cellnum(i),tstamp(i));
    B{cellid,tstamp(i)}=boun{cellnum(i),tstamp(i)};
    pixels{cellid,tstamp(i)}=pxls{cellnum(i),tstamp(i)};%edit for color decode DF
end

%clean length traces for isolated pixels
for i=1:ncells
    lbw=lcell(i,:)>0;
    lbwc=bwmorph(lbw,'clean');
    lcell(i,lbwc==0)=0;
    wcell(i,lbwc==0)=0;
end

%throw away cells with only one or two time points
delind=[];
for i=1:ncells
    if length(nonzeros(lcell(i,:)))<=2 & T>10
        delind=[delind;i];
    end
end
lcell(delind,:)=[];
wcell(delind,:)=[];
B(delind,:)=[];
pixels(delind,:)=[];%edit for color decode DF
[ncells,~]=size(lcell);

lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;

%insert skipped frames
lsk=length(skp);
[lm,ln]=size(l);

for i=1:lsk
     lcell=[lcell(:,1:skp(i)-1),(lcell(:,skp(i))+lcell(:,skp(i)-1))/2,lcell(:,skp(i):end)];
     wcell=[wcell(:,1:skp(i)-1),(wcell(:,skp(i))+wcell(:,skp(i)-1))/2,wcell(:,skp(i):end)];
     B=[B(:,1:skp(i)-1),cell(ncells,1),B(:,skp(i):end)];
end
T=T+lsk;

time=time(1:T);
timeskp=time;
timeskp(skp)=[];
nskp=setdiff((1:T),skp);
lcellt=lcell;
wcellt=wcell;

warning('off','MATLAB:chckxy:IgnoreNaN');
for i=1:ncells
    lcellt(i,skp)=spline(timeskp,lcell(i,nskp),time(skp));
    wcellt(i,skp)=spline(timeskp,wcell(i,nskp),time(skp));
end
warning('on','MATLAB:chckxy:IgnoreNaN');

lcellt(isnan(lcell)==1)=NaN;
wcellt(isnan(lcell)==1)=NaN;
lcell=lcellt;
wcell=wcellt;

%dimsionalize the variables
lcell=lcell*lscale;
wcell=wcell*lscale;

%%
%throw away cells that are too short or too fat or too skinny
lcell(lcell<minL)=NaN;
lcell(wcell>maxW)=NaN;
lcell(wcell<minW)=NaN;
%%

wcell(wcell<minW)=NaN;
wcell(wcell>maxW)=NaN;
wcell(lcell<minL)=NaN;

%%
% %throw away cells that grow too fast
ld=(lcell(:,2:end)-lcell(:,1:end-1));
ld=abs([zeros(ncells,1) ld]);
lcell(ld>maxD)=NaN;
%%

%calculate the strain rate
deltat=time(2:end)-time(1:end-1);
v=(lcell(:,2:end)-lcell(:,1:end-1))./((lcell(:,1:end-1)+lcell(:,2:end))/2);

for i=1:ncells
    v(i,:)=v(i,:)./deltat;
end

%throw out outliers and calculate the average strain and strain rate across
%cells
wav=zeros(1,T);
wstd=zeros(1,T);
vav=zeros(1,T-1);
vstd=zeros(1,T-1);
vste=zeros(1,T-1);
ndp=zeros(1,T-1);
for t=1:T   
    wt=wcell(:,t);
    wt(isnan(wt))=0;
    wav(t)=mean(nonzeros(wt));
    wstd(t)=std(nonzeros(wt));
end

for t=1:T-1
    vt=v(:,t);
    vt2=vt;
    lt=lcell(:,t);
    wt=wcell(:,t);
    vt(isnan(vt))=0;
    vav(t)=mean(nonzeros(vt));
    vstd(t)=std(nonzeros(vt));
    vt(vt==0)=NaN;
    if vstd~=0 %I put this if statement in
        vt2(abs(vt-vav(t))>=3*vstd(t))=NaN;
        lt(abs(vt-vav(t))>=3*vstd(t))=NaN;
        wt(abs(vt-vav(t))>=3*vstd(t))=NaN;
    end
    v(:,t)=vt2;
    lcell(:,t)=lt;
    wcell(:,t)=wt;
    
    vt2(isnan(vt2))=0;
    vav(t)=mean(nonzeros(vt2));
    vstd(t)=std(nonzeros(vt2));
    vste(t)=std(nonzeros(vt2))/sqrt(length(nonzeros(vt2)));
    ndp(t)=length(nonzeros(vt2));
end

if troubleshoot==1
    cd(imdir);
    imfiles=dir('*.tif');
   
    for t=1:T
            t
            imagename=imfiles(t).name;
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
end

%plot data
figure(1), title('Cell Length vs. Time')
clf
hold on
for i=1:ncells 
    plot((1:size(lcell,2))*10/60,lcell(i,:))
end
xlabel('t (min)')
ylabel('l (\mum)')
fig2pretty

% figure(2), title('Cell Width vs. Time')
% hold on
% ciplot(wav-wstd,wav+wstd,time,[0.75 0.75 1])
% plot(time,wav,'-r','LineWidth',2)
% xlabel('time (s)')
% ylabel('width (\mum)')
% fig2pretty

tmid=(time(2:end)+time(1:end-1))/2;

% figure(3), title('Strain Rate vs. Time')
% hold on
% for i=1:ncells
%     plot(tmid,v(i,:))
% end
% plot(tmid,vav,'-r')
% xlabel('t (s)')
% ylabel('\mu (s^{-1})')
% fig2pretty

%%

cd(savedir)
save([basename '_WT'])
