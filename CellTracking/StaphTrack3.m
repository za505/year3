%StaphTrack2.m
%Tracks S. aureus growth from fluorescent image stacks where cell membrane is stained.  
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
%basename='041913_3_Pos2_1';
basename='c1';
dirname=['/Volumes/Lacie/Matlab Ready/022913_1_pos0/' basename];

minA=5;
sigma=3;
srad=1;
sfac=1;
lscale=0.08;
tscale=60;
thresh=0.003;
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

%figure
%track cells
for t=1:T
    t
    
    %load image
    imagename=directory(t).name;  
    im=imread(imagename);
    [imM,imN]=size(im);
    
    %normalize image
    ppix=0.5;
    im=norm16bit(im,ppix);
    imraw=im;
    
    %smooth image with gaussian filter
    H=fspecial('gaussian',sigma,1);
    %im=medfilt2(im);
    im=imfilter(im,H,'replicate');
    
    %create mask
    level=graythresh(im);
    mask=im2bw(im,level);
    mask=imcomplement(mask);
    mask=imerode(mask,strel('disk',1));
    
    srec=strel('disk',srad);
    sel1=strel('disk',1);
    sel2=strel('disk',2);
    sel4=strel('disk',5);
    imc=imcomplement(im);
    imer=imerode(imc,srec);
    imr=imreconstruct(imer,imc);
    imrm=imregionalmax(imr);
    %figure,imshow(imrm),pause
    
    [ime,thresh]=edge(imraw,'log',0.001);
    ime(mask)=0;
    imec=imcomplement(ime);
    cc=bwconncomp(imec,4);
    stats1=regionprops(cc,imraw,'MeanIntensity');
    [~,idx]=max([stats1.MeanIntensity]);
    bw2=ismember(labelmatrix(cc),idx);
    bw3=bwmorph(bw2,'skel',2);
    bw3=bwmorph(bw3,'spur');
    
    bw4=imcomplement(bw3);
    cc=bwconncomp(bw4,4);
    stats1=regionprops(cc,'Area');
    idx=find([stats1.Area]<200&[stats1.Area]>minA);
    bw4=ismember(labelmatrix(cc),idx);
    %figure,imshow(bw4),pause
    
    
%     figure,imshow(ime),pause
%     imd=imdilate(ime,sel1);
%     imdc=imcomplement(imd);
%     %figure,imshow(imdc),pause
%     cc=bwconncomp(imdc,8);
%     stats1=regionprops(cc,'Area');
%     stats2=regionprops(cc,imrm,'MaxIntensity');
%     idx=find([stats1.Area]<100&[stats1.Area]>minA&[stats2.MaxIntensity]==1);
%     %idx=find([stats1.Area]<100&[stats1.Area]>minA);
%     bw2=ismember(labelmatrix(cc),idx);
   
    clear cents
    clear orients
    P=cell(length(idx),1);
    for i=1:length(idx)
        bwt=ismember(labelmatrix(cc),idx(i));
        bwt=imdilate(bwt,sel1);
        stats=regionprops(bwt,'Centroid','Orientation');
        [P(i),~]=bwboundaries(bwt,4,'noholes');    
        orients(i)=[stats.Orientation];
        cents(i,:)=cat(1,stats.Centroid);
    end
    nc(t)=length(idx);
% %    
%     figure
%    imshow(imraw)
%    hold on
%    for k=1:length(P)
%        plot(P{k}(:,2),P{k}(:,1),'-r');      
%    end
%    pause
    
%    stats=regionprops(bw2,'Area','Centroid','Orientation');
%    stats2=regionprops(bw2,imraw,'MeanIntensity','PixelValues');
    
    nc(t)=length(P);
%    areas=[stats.Area];
  %  cents=cat(1,stats.Centroid);
  %  orients=[stats.Orientation];
    V1(nc(t),t)=0;
    V2(nc(t),t)=0;
    %intensity(nc(t),t)=0;
    %intensity(1:nc(t),t)=[stats2.MeanIntensity]';
%    pixels{nc(t),t}=[];
%    a(1:nc(t),t)=[stats.Area]';
   centroids=cents;
   
    for n=1:nc(t)
        
        %smooth cell contours
        rP=[P{n}(:,2),P{n}(:,1)];
        px=[rP(1:end-1,1);rP(1:end-1,1);rP(:,1)];
        py=[rP(1:end-1,2);rP(1:end-1,2);rP(:,2)];
        sp=length(rP);
        dS=sqrt(diff(px).^2+diff(py).^2);
        S=[0 cumsum(dS)'];
  
        px=csaps(S,px,0.5,S);
        py=csaps(S,py,0.5,S);
        
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
%        pixels{n,t}=[stats2(n).PixelValues];
        
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
%    figure,hold on,axis equal, axis off, axis ij,colormap('hot')
%    for n=1:nc(t)
%         surf(shellx{n,t},shelly{n,t},shellz{n,t})      
%    end
   
   pause(0.1)
  
   
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
time=time(1:T);

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
%    acell(cellid,tstamp(i))=a(cellnum(i),tstamp(i));
    vcell1(cellid,tstamp(i))=V1(cellnum(i),tstamp(i));
    vcell(cellid,tstamp(i))=V2(cellnum(i),tstamp(i));
%    icell(cellid,tstamp(i))=intensity(cellnum(i),tstamp(i));
    B{cellid,tstamp(i)}=boun{cellnum(i),tstamp(i)};
%    P{cellid,tstamp(i)}=pixels{cellnum(i),tstamp(i)};
end

for k=1:ncells
    if isempty(B{k,T})~=1
        text(B{k,T}(1,1),B{k,T}(1,2),num2str(k),'Color',[0 0 1])
    end
end

%dimensionalize the variables
%acell=acell*lscale^2;
vcell1=vcell1*lscale^3;
vcell=vcell*lscale^3;

%atot=sum(acell);
vtot1=sum(vcell1);
vtot=sum(vcell);

%acell(acell==0)=NaN;
vcell1(vcell1==0)=NaN;
vcell(vcell==0)=NaN;
%icell(icell==0)=NaN;

end

%plot data
figure, title('Cell Volume vs. Time')
hold on
for i=1:ncells
    plot(time,vcell(i,:),'-m')
    %text(time,vcell(i,:),num2str(i))
end
plot(time,vtot,'-c')
xlabel('t (s)')
ylabel('V (\mum^3)')
fig2pretty


cd /Users/Rico/Documents/MATLAB
save([basename '_ST'])
