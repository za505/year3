%PlasmolysisTrack.m
%Tracks and quantifies size of plasmolysis bays in bacteria
%To be used after imagealign.m, Eraseimagepart.m, & BacTrack.m. Based on
%LinTrack.m
%Zarina Akbary, last updated 04/25/2021
clear 
close all

%input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename='05262021_Exp1';
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05262021_analysis/' basename '/05262021_plasmolysisTrack'];%Directory to save the output .mat file to.
%frame=4; %frame immediately upon hyperosmotic shock
channels={'GFP', 'TADA'};
vis=0; %to visualize the boundaries of the pre-shock cells plotted in post-shock frames
preShock=2;
postShock=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load BackTrack.m data
%TADA data
%cd([dirname '/' basename '_TADA/' basename '_figures']);
cd(dirname)
load([basename '_BTtada'], 'B', 'time', 'pixels', 'directory', 'T', 'ncells', 'acell', 'im', 'mlines');

pixelsTADA=pixels;
%acellTADA=acell;
B_TADA=B;
directory_TADA=directory;
ncellsT=ncells;

%GFP data
%cd([dirname '/' basename '_GFP/' basename '_figures']);
%load([basename '_BTgfp'], 'B', 'time', 'pixels', 'directory', 'T', 'ncells', 'acell', 'im');
load([basename '_BTgfp'], 'directory');

directory_GFP=directory;
%B_GFP=B;
%pixelsGFP=pixels;

%% start with the pre-shock cells
cd(directory_TADA(1).folder);
im1=imread(directory_TADA(preShock).name);

cd(directory_GFP(1).folder);
im2=imread(directory_GFP(preShock).name);

%create a binary image based on the pre-shock TADA frame
bw1=zeros(size(im));
for n=1:ncells
    bw1(pixelsTADA{n, preShock})=1;
end

%look only at TADA-tracked cells in the GFP frame. Threshold using the Otsu method
bw2=imbinarize(im2, graythresh(im2(bw1==1)));

%what is the overlap/nonoverlap between the two?
is1=bw1-bw2;
%imshow(is1)

%the conversion is 0.08 microns per pixel, and so this gap is unlikely to
%be actual. What happens if we erode bw1 and dilate bw2 by one pixel?
%nhood=structuring element neighborhood
nhood=[0,1,0;1,1,1;0,1,0];
bw1=imerode(bw1, nhood);
bw2=imdilate(bw2, nhood);
is2=bw1-bw2;
%imshowpair(is1, is2, 'montage')

%hm...still not quite right. there are -1 value pixels, but that's not the
%main problem
bw1(bw1==1)=2;
is3=bw1-bw2;
is3(is3==-1)=1;
is3(is3==2)=0;
imshow(is3)
pause
%% let's try this strategy on the post-Shock cells
cd(directory_TADA(1).folder);
im1=imread(directory_TADA(postShock).name);

cd(directory_GFP(1).folder);
im2=imread(directory_GFP(postShock).name);

%create a binary image based on the pre-shock TADA frame
bw1=zeros(size(im));
for n=1:ncells
    bw1(pixelsTADA{n, postShock})=1;
end

%look only at TADA-tracked cells in the GFP frame. Threshold using the Otsu method
bw2=imbinarize(im2, graythresh(im2(bw1==1)));

%what is the overlap/nonoverlap between the two?
is1=bw1-bw2;
%imshow(is1)

%the conversion is 0.08 microns per pixel, and so this gap is unlikely to
%be actual. What happens if we erode bw1 and dilate bw2 by one pixel?
%nhood=structuring element neighborhood
nhood=[0,1,0;1,1,1;0,1,0];
bw1=imerode(bw1, nhood);
bw2=imdilate(bw2, nhood);
is2=bw1-bw2;
%imshowpair(is1, is2, 'montage')

%hm...still not quite right. there are -1 value pixels, but that's not the
%main problem
bw1(bw1==1)=2;
is3=bw1-bw2;
is3(is3==-1)=1;
is3(is3==2)=0;
imshow(is3)
pause
%% let's focus only on the post-Shock cells
cd(directory_TADA(1).folder);
im1=imread(directory_TADA(postShock).name);

cd(directory_GFP(1).folder);
im2=imread(directory_GFP(postShock).name);

%create a binary image based on the post-Shock TADA frame
bw1=zeros(size(im));
for n=1:ncells
    bw1(pixelsTADA{n, postShock})=1;
end

%look only at TADA-tracked cells in the GFP frame. Threshold using the Otsu method
bw2=imbinarize(im2, graythresh(im2(bw1==1)));

%nhood=structuring element neighborhood
nhood=[0,1,0;1,1,1;0,1,0];
bw1=imerode(bw1, nhood);
bw1Shell=bwmorph(bw1, 'skel');

bw1=imdilate(bw1, nhood);
%bw1=imerode(bw1, strel('disk',4));
for n=1:ncells
    midx=sub2ind(size(im), round(mlines{n, postShock}(:,1)), round(mlines{n, postShock}(:,2)));
    bw1(midx)=0;
end
%strel('disk',cellLink)
%for some reason, the septum's get merges/ignored in the binary. Since they
%are usually the brightest regions, let's see if we can find them
%sanity check
% figure(1), imshowpair(im1, im2) %title('Orignal TADA vs GFP')
% figure(2), imshowpair(bw1, bw2) %title('Binary TADA vs GFP')
% figure(3), imshowpair(im1, bw1) %title('TADA: original vs binary')
% figure(4), imshowpair(im2, bw2) %title('GFP: original vs binary')

%% identify plasmolysis bays
bw3=bw1-bw2;

%sanity check
figure(1), imshowpair(im1, im2) %title('Orignal TADA vs GFP')
figure(2), imshowpair(bw3, im2) %title('P.bays vs original GFP')
%% first, let's look at our TADA cells
cd(directory_TADA(1).folder);
im1=imread(directory_TADA(preShock).name);
im2=imread(directory_TADA(postShock).name);

% imshowpair(im1, im2, 'montage')
% pause, close

%what is the size difference b/t the pre and post-shock cells. This is the
%difference that will be negated when looking at the difference between
%pre- and post-shock GFP cells

%potential prob: what if the # of cells tracked in the pre and post
%shock TADA images don't match? have to add a line of code that takes this into
%account

%potential solution
%create a matrix of all the pixels indicies of the tracked TADA cells
allpixels=[];
allpixels2=[];
for k=1:ncells
    allpixels=[allpixels; pixels{k,preShock}];
    allpixels2=[allpixels2; pixels{k,postShock}];
end

pmem=allpixels(ismember(allpixels, allpixels2));
pdiff=setdiff(allpixels, allpixels2);

%this is the amount of area loss we expect from contraction
bw1=zeros(size(im));
bw1(pmem)=1;
bw1(pdiff)=1; %this is the full pre-shock TADA cell
bw2=bwmorph(bw1, 'remove'); %this is my work-around because 'skel' isn't working the way I expect it to
bw2(pdiff)=1;
%imshow(bw2)
%% now, let's look at the GFP cells
cd(directory_GFP(1).folder);
im3=imread(directory_GFP(preShock).name);
im4=imread(directory_GFP(postShock).name);

% imshowpair(im3, im4, 'montage')
% pause, close

%we are only looking at pixels in the pre-shock TADA-tracked
%region. Threshold using the Otsu method
bw3=imbinarize(im3, graythresh(im3(bw1==1)));
bw4=imbinarize(im4, graythresh(im4(bw1==1)));

%imshow(bw3)
%imshow(bw4)

bdiff=bw3-bw4; %this is the area difference between the pre- and post-shock GFP cells
bdiff2=bdiff-bw2; %hm, bw2 is bigger than bdiff. the overlap isn't perfect
bdiff2(bdiff2==-1)=0; %these are the putative plasmolysis bays

bdiff2=bwmorph(bdiff2, 'spur');
bdiff2=bwmorph(bdiff2, 'clean');
bdiff2=bwmorph(bdiff2, 'open');

%sanity check
%imshowpair(im4, bdiff2)

%% measure the area of the plasmolysis bays
baysCC=bwconncomp(bdiff2,8);
baysStats=regionprops(baysCC,im4,'Area','MeanIntensity', 'PixelIdxList');
idx=find([baysStats.Area]>1);
bays=ismember(labelmatrix(baysCC),idx);

%sanity check
figure(1), imshowpair(im4, im2)
figure(2), imshowpair(im4, bays)
%% measure the frequency, severity, and localization of plasmolysis
%firstly, how many cells do we have? (check pre-shock GFP binary image)
% D = bwdist(~bw3);
% D = -D;
% cells = watershed(D);
% cells(~bw3)=0;
% 
% rgb = label2rgb(cells,'jet',[.5 .5 .5]);
% imshow(rgb)

%assign each plasmolysis bay to a cell
cellsBW=bwmorph(bw3, 'open'); %does this line make much of a difference?
cellsL=bwlabel(cellsBW);

cellsCC=bwconncomp(cellsBW,8);
cellsStats=regionprops(cellsCC, im3, 'Area','MeanIntensity', 'PixelIdxList', 'Centroid');

ncells=max(max(cellsL));

pcells=[];
for b=1:height(baysStats)
    for n=1:ncells
        if sum(ismember(baysStats(b).PixelIdxList, cellsStats(n).PixelIdxList))/height(baysStats(b).PixelIdxList)>0.95
            baysStats(b).cellID=n;
            pcells=[pcells, n];
        end
    end
end

%measure the frequency of plasmolysis (# plasmolyzed cells/total cells)
freq=sum(ismember(1:ncells, baysStats.cellID))/ncells;

%measure the severity of plasmolysis (plasmolysis area/total pre-shock area of plasmolyzed cells)
severity=nan(length(pcells),1);
for j=1:length(pcells)
    n=pcells(j);
    severity(n)=sum(baysStats(baysStats.cellID==n).Area)/cellsStats(n).Area;
end


%sanity check:
%from the GFP image, set any pixels value not within the TADA region to 0
didx=find(bwe==0);
im3(didx)=0;
im4(didx)=0;
%okay, so the outline of our tada cell matches the outline of the GFP cell
%so far. Do they cancel?
im3(im3~=0)=1;
im3=imbinarize(im3);
bwf=im3-bwe;
imshow(bwf)

bw1=zeros(size(im));
for k=1:ncells
    bw1(im1(pixels{k,preShock})>1200)=1;
end
%% overlay TADA cell boundaries on GFP image stack
cd(directory_GFP(1).folder);

if vis==1
for t=1:T
    t
    
    %Load image
    imagename=directory_GFP(t).name;

    im=imread(imagename);
    [imM,imN]=size(im);
    
%     cd(directory_TADA(1).folder);
%     im1=imread(directory_TADA(postShock).name);

%     cd(directory_GFP(1).folder);
%     im2=imread(directory_GFP(postShock).name);

    figure
    %imshowpair(im1, im2)
    imshow(im, [])
    hold on
    for k=1:ncellsT
       if isempty(B_TADA{k,t})==0
%            [r1 c1]=ind2sub(size(im), pixels{k,t});
%            plot(c1, r1)
             plot(B_TADA{k,t}(:,1),B_TADA{k,t}(:,2),'-r')
       end
       
%        if isempty(B_GFP{k,t})==0
% %            [r1 c1]=ind2sub(size(im), pixels{k,t});
% %            plot(c1, r1)
%              plot(B_GFP{k,t}(:,1),B_GFP{k,t}(:,2),'--g')
%        end
    end
    
  pause
  close all
  
end

end

%from this, we know there are two "gaps" b/t the TADA and GFP images:
%around the cells and in the middle (at the septum)

%% Generate binary GFP and TADA preShock images
cd(directory_TADA(1).folder);
im1=imread(directory_TADA(preShock).name);

cd(directory_GFP(1).folder);
im2=imread(directory_GFP(preShock).name);

bw1=zeros(size(im));
for i=1:ncellsT
    if isempty(pixelsTADA{i,preShock})==0
        bw1(pixelsTADA{i,preShock})=1;
    end
end

%try watershedding
%https://www.mathworks.com/help/images/ref/watershed.html
D = bwdist(~bw1);
D = -D;
L = watershed(D);
L(~bw1) = 0;

%it didn't work

%try erosion
%https://www.mathworks.com/help/images/ref/imerode.html
%nhood=structuring element neighborhood
nhood=[0,1,0;1,1,1;0,1,0];
bwe=imerode(bw1, nhood);

bw2=zeros(size(im));
im3=imadjust(im2,[0 55000]/2^16,[]);   
bw2=imbinarize(im3, graythresh(im3));

nhood=[0,1,0;1,1,1;0,1,0];
nhood2=[0 0 1 0 0; 0 0 1 0 0; 1 1 1 1 1; 0 0 1 0 0; 0 0 1 0 0];
bw3=imdilate(bw2, nhood2);

bw4=bw1-bw2;
bw5=bwe-bw3;

%imshowpair(bw5, bw4, 'montage')

%% Generate binary GFP and TADA postShock images
imscd(directory_TADA(1).folder);
im1=imread(directory_TADA(postShock).name);

cd(directory_GFP(1).folder);
im2=imread(directory_GFP(postShock).name);

bw1=zeros(size(im));
for i=1:ncellsT
    if isempty(pixelsTADA{i,postShock})==0
        bw1(pixelsTADA{i,postShock})=1;
    end
end

bw2=zeros(size(im));
for i=1:ncells
    if isempty(pixelsGFP{i,postShock})==0
        bw2(pixelsGFP{i,postShock})=1;
    end
end

nhood=[0,1,0;1,1,1;0,1,0];
nhood2=[0 0 1 0 0; 0 0 1 0 0; 1 1 1 1 1; 0 0 1 0 0; 0 0 1 0 0];
bw3=imdilate(bw2, nhood2);

bw4=bw1-bw2;
bw5=bw1-bw3;

%% develop a method to generate in "expected" GFP cell based on a TADA frame
bw1=zeros(size(im));
for k=1:ncells
    bw1(im1(pixels{k,preShock})>1200)=1;
end

bw2=zeros(size(im));
for k=1:ncells
    bw1(pixels{k,preShock})=1;
end
%% Find the difference in cell area between TADA and GFP preShock cells
preShock_area=nan(ncellsT, preShock);
for i=1:preShock
    preShock_area(:,i)=acellTADA(:,i)-acellGFP(:,i);
end

%incredulously, there is no difference in area
%given this, we will assume the TADA area is the predicted GFP area if
%there is no plasmolysis

%load a postShock GFP image
postShock_cytoname=directory_GFP(postShock).name;
postShock_cytoim=imread(postShock_cytoname);


cd([dirname '/' basename '_phase/' basename '_erased'])

%%
for t=[1 110 120 145]
    t
    
    %Load image
    imagename=directory(t).name;

    im=imread(imagename);
    
   figure
   imshow(im)
   hold on
   for k=1:ncells
       if isempty(B{k,t})==0
        plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
       end
   end
    
  pause
  close all
  
end



%load the image of the labeled cytoplasm
cd(dirname)
directory=dir('*.tif');
cyimagename=directory(1).name;
cyimage=imread(cyimagename);

%load the image of the labeled cell wall
cwimagename=directory(2).name;
cwimage=imread(cwimagename);

%now overlap the images
merge_fc=imfuse(cyimage, cwimage, 'falsecolor');
merge_diff1=imfuse(cyimage, cwimage, 'diff');
merge_diff2=imfuse(cwimage, cyimage, 'diff');

%Load image
im=imread(imagename);
%process the phase image
[imp, impc, impb]=imageProcess([dirname '/' basename '_' channels{1}], frame, 1);
phase=cat(3, imp, impc, impb);
montage(phase)
pause

%process the TADA image
[imt, imtc, imtb]=imageProcess([dirname '/' basename '_' channels{2}], frame, 0);
tada=cat(3, imt, imtc, imtb);
montage(tada)
pause

function [im, imc, imb] = imageProcess(dirname, frame, complement)
    
    %change directory
    cd(dirname)
    directory=dir('*.tif');
    imagename=directory(frame).name;
    
    %Load image
    im=imread(imagename);
    %[imM,imN]=size(im);
    
    %De-speckle image
    im=medfilt2(im); %each pixel is replaced with the median value of its neighboring pixels (in this case 3x3)

    %Normalize images
    ppix=0.5; %I think this brings all the pixels in each image into a similar range
    im=norm16bit(im,ppix);
     
    if complement==1
        %Enhance contrast
        imc=imcomplement(im); %flips the pixel values (cells go from black to white)
    else
        imc=im;
    end
    
    [imcounts,bins]=imhist(imc);
    [imcounts,idx]=sort(imcounts);
    bins=bins(idx);
    thresh1=bins(end-1); 
    
    if thresh1==65535
        thresh1=bins(end-2);
    end
    
    imb=imadjust(imc,[thresh1/65535 1],[]); %basically the background pixels get lumped together and get set to 0 (black background) 
end
       

