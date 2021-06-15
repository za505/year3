%PlasmolysisTrack.m
%Tracks and quantifies size of plasmolysis bays in bacteria
%To be used after imagealign.m, Eraseimagepart.m, & BacTrack.m. Based on
%LinTrack.m
%Zarina Akbary, last updated 04/25/2021
clear 
close all

%input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename='05262021_Exp1';
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/05262021_analysis/' basename];%Directory to save the output .mat file to.
%frame=4; %frame immediately upon hyperosmotic shock
channels={'GFP', 'TADA'};
vis=2; %to visualize the boundaries of the pre-shock cells plotted in post-shock frames
preShock=2;
postShock=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load BackTrack.m data
%TADA data
cd([dirname '/' basename '_TADA/' basename '_figures']);
load([basename '_BTtada'], 'B', 'time', 'pixels', 'directory', 'T', 'ncells', 'acell', 'im');

%pixelsTADA=pixels;
%acellTADA=acell;
%B_TADA=B;
directory_TADA=directory;
%ncellsT=ncells;

%GFP data
cd([dirname '/' basename '_GFP/' basename '_figures']);
load([basename '_BTgfp'], 'directory');

directory_GFP=directory;

%% overlay TADA cell boundaries on GFP image stack
cd(directory_GFP(1).folder);

if vis==1
for t=1:T
    t
    
    %Load image
    imagename=directory_GFP(t).name;

    im=imread(imagename);
    [imM,imN]=size(im);
    
    figure
    imshow(im, [])
    hold on
    for k=1:ncells
       if isempty(pixels{k,t})==0
%            [r1 c1]=ind2sub(size(im), pixels{k,t});
%            plot(c1, r1)
             plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
       end
    end
    
  pause
  close all
  
end

elseif vis==2
    %Load gfp post-shock image
    cd(directory_GFP(postShock).folder)
    imagename=directory_GFP(postShock).name;
    im1=imread(imagename);
    
    %Load tada post-shock image
    cd(directory_TADA(postShock).folder)
    imagename=directory_TADA(postShock).name;
    im2=imread(imagename);
    
    minVal1=[];
    maxVal1=[];
    medVal1=[];
    modVal1=[];
    for k=1:ncells
        tempCol=im1(pixels{k,postShock});
        minVal1=[minVal1 max(min(tempCol))];
        maxVal1=[maxVal1 min(max(tempCol))];
        medVal1=[medVal1 min(median(tempCol))];
        modVal1=[modVal1 min(mode(tempCol))];
    end
    
    minVal2=[];
    maxVal2=[];
    for k=1:ncells
        tempCol=im2(pixels{k,postShock});
        minVal2=[minVal2 max(min(tempCol))];
        maxVal2=[maxVal2 min(max(tempCol))];
    end
end

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
       

