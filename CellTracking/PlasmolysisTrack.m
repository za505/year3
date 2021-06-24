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
load([basename '_BTtada'], 'B', 'time', 'pixels', 'directory', 'T', 'ncells', 'acell', 'im', 'lscale');

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

%erode the TADA cell outline by two pixels
nhood=[0 0 1 0 0; 0 0 1 0 0; 1 1 1 1 1; 0 0 1 0 0; 0 0 1 0 0];
bw1=imerode(bw1, nhood);

%% look at the post-shock cells
cd(directory_TADA(1).folder);
im3=imread(directory_TADA(postShock).name);

cd(directory_GFP(1).folder);
im4=imread(directory_GFP(postShock).name);
    
%create a binary image based on the pre-shock TADA frame
bw3=zeros(size(im));
for n=1:ncells
    bw3(pixelsTADA{n, postShock})=1;
end

%look only at TADA-tracked cells in the GFP frame. Threshold using the Otsu method
bw4=imbinarize(im4, graythresh(im4(bw3==1)));

%erode the TADA cell outline by two pixels
bw3=imerode(bw3, nhood);

%% what is the difference between pre- and post-shock frames
is1=bw1-bw3; 
is1(is1==-1)=0;
is2=bw2-bw4;
is2(is2==-1)=0;

%clean and de-spur
is1=bwmorph(is1, 'clean');
is1=bwmorph(is1, 'spur');
is1=bwmorph(is1, 'bridge');

is2=bwmorph(is2, 'clean');
is2=bwmorph(is2, 'spur');
is2=bwmorph(is2, 'bridge');

%subtract the change in area from contraction (TADA) from plasmolysis (GFP)
is3=is2-is1;
is3(is3==-1)=0;

is3=bwmorph(is3, 'clean');
is3=bwmorph(is3, 'spur');

%get the stats
cc1=bwconncomp(is3,8);
stats1=regionprops(cc1,im4,'Area','MeanIntensity', 'PixelIdxList');
idx1=find([stats1.Area]>2);
by1=ismember(labelmatrix(cc1),idx1);

%which cell does the plasmolysis bay belong to?
stats=stats1(idx1); %remember these are the ones we care about

pcells=[];
for b=1:height(stats)
    for n=1:ncellsT
        if (sum(ismember(stats(b).PixelIdxList, pixelsTADA{n, postShock}))/length(stats(b).PixelIdxList))>=0.95
            stats(b).cellID=n;
            pcells=[pcells, n];
        end
    end
end

bays=nan(ncellsT,1);
for n=1:ncellsT
    idxN=find([stats.cellID]==n);
    
    if isempty(idxN)==1
        bays(n)=NaN;
    else
        totalA=sum([stats(idxN).Area]);
        bays(n)=totalA;
    end
end

%% plot the comparison
% figure
% bar([1:ncellsT], bays*lscale^2)
% xlabel('Cell ID')
% ylabel('Total Plasmolyzed Area (\mum^{2})')
cd(dirname)

figure
scatter([stats.cellID], [stats.Area]*lscale^2, 'filled')
xticks([1:ncellsT])
xlim([0,3])
xlabel('Cell ID')
ylabel('Plasmolyzed Area (\mum^{2})')
saveas(gcf, [basename,'_plasmolysis1.fig'])
saveas(gcf, [basename,'_plasmolysis1.png'])

save(['plasmolysis1'])