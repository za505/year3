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
cd /Users/Rico/Documents/MATLAB/data

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='041913_3_Pos2_1';
dirname=['/Matlab Ready/' basename '/' basename '_a' ];
%metaname=['/Volumes/HD-PNTU3/Matlab Ready/' basename '/' basename '.txt'];
hsize=4;
sigma=2;
sfac=1;
lscale=0.08;
tscale=60;
recrunch=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if recrunch==1
    load([basename '_ST'])
else

%determine number of frames
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

figure
%track cells
for t=1:T
    t
    
    %load image
    imagename=directory(t).name;  
    im=imread(imagename);
    imraw=im;
    [imM,imN]=size(im);
    %figure,imshow(im,[])
    
    %smooth image with gaussian filter
    H=fspecial('gaussian',hsize,sigma);
    im=imfilter(im,H,'replicate');

    im=imcomplement(im);
    %figure,imshow(im,[]),pause
    
    %normalize image
    ppix=0.5;
    im=norm16bit(im,ppix);
    imraw=norm16bit(imraw,ppix);
    
    %identify cells
    rm=imregionalmax(im);
    se=strel('disk',1);
    rm=imdilate(rm,se); 
    
    %find septa
    iml=del2(double(im));
    H2=fspecial('gaussian',4,2);
    iml=imfilter(iml,H2,'replicate');
    iml=iml>0;
    
    %clean septa
    iml=bwmorph(iml,'thin',10);
    iml=bwmorph(iml,'spur');

    %find colony outline
    thresh=graythresh(im);
    BW=im2bw(im,thresh);
    
    %perform watershed
    wsr=iml==BW;
    wsr(BW==0)=0;
    rm(BW==0)=0;
    imc=imcomplement(im);
    imc(BW==0)=-Inf;
    imc(rm)=-Inf;
    imc(wsr)=Inf;
    L=watershed(imc);
    bw2=L>1;
    
    %calculate cell boundaries
    cc=bwconncomp(bw2,8);
    stats=regionprops(cc,imc,'Area');
    idx=find([stats.Area]>8);
    bw2=ismember(labelmatrix(cc),idx);

    P=cell(length(idx),1);
    for i=1:length(idx)
        bwt=ismember(labelmatrix(cc),idx(i));
        bwt=imdilate(bwt,[0 1 0;1 1 1;0 1 0]);
        [P(i),~]=bwboundaries(bwt,4,'noholes');
    end
    
    stats=regionprops(bw2,'Area','Centroid','Orientation');
    stats2=regionprops(bw2,imraw,'MeanIntensity','PixelValues');
    
    nc(t)=length(P);
    areas=[stats.Area];
    cents=cat(1,stats.Centroid);
    orients=[stats.Orientation];
    a(nc(t),t)=0;
    V1(nc(t),t)=0;
    V2(nc(t),t)=0;
    intensity(nc(t),t)=0;
    intensity(:,t)=[stats2.MeanIntensity];
    pixels{nc(t),t}=[];
    a(1:nc(t),t)=[stats.Area]';
    centroids=cents;
   
    for n=1:nc(t)
        
        %smooth cell contours
        bwt=ismember(labelmatrix(cc),idx(n));
        bw2=bwmorph(bwt,'thicken',1);
         
        rP=[P{n}(:,2),P{n}(:,1)];
        px=[rP(1:end-1,1);rP(1:end-1,1);rP(:,1)];
        py=[rP(1:end-1,2);rP(1:end-1,2);rP(:,2)];
        sp=length(rP);
        dS=sqrt(diff(px).^2+diff(py).^2);
        S=[0 cumsum(dS)'];
  
        px=csaps(S,px,0.1,S);
        py=csaps(S,py,0.1,S);
        
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
        pixels{n,t}=[stats2(n).PixelValues];
        
        %calculate volume
        phi=orients(n)*pi/180;
        nx=nx-cents(n,1);
        ny=ny-cents(n,2);
        
        [nx,ny]=rotatexy(nx',ny',phi);
         
        nxp=nx;
        nyp=ny;
        
        %find poles
        [~,p1]=min(abs(ny));
        nyp(sign(nx)==sign(nx(p1)))=1e5;
        [~,p2]=min(abs(nyp));
        
        [~,p3]=min(abs(nx));
        nxp(sign(ny)==sign(ny(p3)))=1e5;
        [~,p4]=min(abs(nxp));
   
       
        %shift outlines such that the first point is pole
        nx1=circshift(nx(1:end-1),-p1+1);
        ny1=circshift(ny(1:end-1),-p1+1);
        nx1=[nx1;nx1(1)];
        ny1=[ny1;ny1(1)];
        p2=mod(p2-p1,ls);
        
        nx2=circshift(nx(1:end-1),-p3+1);
        ny2=circshift(ny(1:end-1),-p3+1);
        nx2=[nx2;nx2(1)];
        ny2=[ny2;ny2(1)];
        p4=mod(p4-p3,ls);
        
        %create mesh
        s1=cellmesh(nx1,ny1,p2);
        s2=cellmesh(nx2,ny2,p2);
        V1(n,t)=s1.V;
        V2(n,t)=s2.V;
  
        %calculate 3d surface
        Res=50;
        [Sx1,Sy1,Sz1]=rotatemesh(s1,Res);
        [Sx2,Sy2,Sz2]=rotatemesh(s2,Res);
        
        [Sx1,S1]=rotatexy(Sx1,Sy1,-phi);
        [Sx2,Sy2]=rotatexy(Sx2,Sy2,-phi);
        
        shellx1{n,t}=Sx1+cents(n,1);
        shelly1{n,t}=Sy1+cents(n,2);
        shellz1{n,t}=Sz1;
        shellx{n,t}=Sx2+cents(n,1);
        shelly{n,t}=Sy2+cents(n,2);
        shellz{n,t}=Sz2;
   
    end

    allcentroids=[allcentroids;centroids];
    tstamp=[tstamp;ones(nc(t),1)*t];
    cellnum=[cellnum;(1:nc(t))'];
    
   clf
   imshow(imraw)
   hold on
   for k=1:nc(t)
       plot(boun{k,t}(:,1),boun{k,t}(:,2),'-r');      
   end
   pause(.1)
   toc
end

if exist('metaname')==1
    if exist(metaname)==2
        %extract timepoints from metadata
        tpoints=metadata(metaname);
        time=tpoints(1:2:end);
        time2=tpoints(2:2:end);
        
        %fix bug where micromanager screws up its timing
        dtime=diff(time);
        fdt=find(dtime>2*(dtime(1)));
        if isempty(fdt)~=1
            fdt=fdt(1);
            time(fdt+1:end)=time(fdt+1:end)-dtime(fdt)+dtime(fdt-1);
            time2(fdt+1:end)=time2(fdt+1:end)-dtime(fdt)+dtime(fdt-1);
            
            time=time(1:T);
            time2=time2(1:T);
        end
    else
        time=[0:T-1]*tscale;
    end
else
    time=[0:T-1]*tscale;
end

%track cells frame to frame
tracks=zeros(size(im));
rcents=round(allcentroids);
linind=sub2ind(size(im),rcents(:,2),rcents(:,1));
tracks(linind)=1;

nhood=[0,1,0;1,1,1;0,1,0];
tracks=imdilate(tracks,strel('disk',1));
overlay1=imoverlay(im,tracks,[.3 1 .3]);

[tracksL,ncells]=bwlabel(tracks,8);

acell=zeros(ncells,T);
vcell1=zeros(ncells,T);
vcell=zeros(ncells,T);
icell=zeros(ncells,T);
B=cell(ncells,T);
P=cell(ncells,T);
lcents=length(allcentroids);

for i=1:lcents
    cellid=tracksL(linind(i));
    acell(cellid,tstamp(i))=a(cellnum(i),tstamp(i));
    vcell1(cellid,tstamp(i))=V1(cellnum(i),tstamp(i));
    vcell(cellid,tstamp(i))=V2(cellnum(i),tstamp(i));
    icell(cellid,tstamp(i))=intensity(cellnum(i),tstamp(i));
    B{cellid,tstamp(i)}=boun{cellnum(i),tstamp(i)};
    P{cellid,tstamp(i)}=pixels{cellnum(i),tstamp(i)};
end

%dimensionalize the variables
acell=acell*lscale^2;
vcell1=vcell1*lscale^3;
vcell=vcell*lscale^3;

atot=sum(acell);
vtot1=sum(vcell1);
vtot=sum(vcell);

acell(acell==0)=NaN;
vcell1(vcell1==0)=NaN;
vcell(vcell==0)=NaN;
icell(icell==0)=NaN;

end

%plot data
figure, title('Cell Volume vs. Time')
hold on
for i=1:ncells
    plot(time,vcell(i,:),'-m')
end
plot(time,vtot,'-c')
xlabel('t (s)')
ylabel('V (\mum^3)')
fig2pretty

figure, title('Mean Cell Intensity vs. Time')
hold on
for i=1:ncells
    plot(time,icell(i,:),'-m')
end
yl=get(gca,'Ylim');
ylim([0 yl(2)])
xlabel('t (s)')
ylabel('Mean Intensity (A.U.)')
fig2pretty

cd /Users/Rico/Documents/MATLAB
save([basename '_ST'])
