%plasImage.m
%Loads the PT.mat data generated from the PlasTrack.m code to create an
%overlay of the post-shock image from the phase and GFP channel. Outlines
%of individual cells of interest can be seen
%Zarina Akbary

clear 
close all

%% load the PT.mat data
basenames={'12162021_Exp1', '12162021_Exp2', '12162021_Exp3', '12162021_Exp4'};

%chose the experiment of interest
basename=basenames{4};
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony2'];

cd(savedir)
load([basename '_colony2_PT.mat'])

%overlay the phase and GFP images
comp=imfuse(imGpost, imPpost); %green = imGpost, magenta = imPpost

%save the overlay image
imwrite(comp, [basename '_colony2_comp.tiff'])


%% see the outline of the cell of interest
%input cell number
n=7;

figure
imshow(comp) 
hold on
plot(round(postB{n,1}(:,1)),round(postB{n,1}(:,2)),'-w')  
           
