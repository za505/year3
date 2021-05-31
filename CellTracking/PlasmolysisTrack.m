%PlasmolysisTrack.m
%Tracks and quantifies size of plasmolysis bays in bacteria
%To be used after imagealign.m, Eraseimagepart.m, & BacTrack.m. Based on
%LinTrack.m
%Zarina Akbary, last updated 04/25/2021
clear 
close all

%input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename='04262021_Exp1_colony1';
dirname=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/04262021_analysis/04262021_Exp1/04262021_Exp1_colony1/plasmolysis_test'];%Directory to save the output .mat file to.
frame=4; %frame immediately upon hyperosmotic shock
channels={'phase', 'TADA'};
vis=1; %to visualize the boundaries of the pre-shock cells plotted in post-shock frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        