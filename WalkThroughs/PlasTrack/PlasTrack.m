%PlasTrack.m
%Tracks percent plasmolyis from phase/CY5 image stacks.
%Dylan Fitzmaurice
%edit: Zarina Akbary
%purpose: I want to compare the plasmolysis quantification b/t TADA/phase,
%TADA/GFP, and/or phase/TADA

clear 
close all

%% load the WallTrack.m data
basename='12162021_Exp3';
cwdir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony4/' basename '_TADA/' basename '_figures'];
phasedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony4/' basename '_phase/' basename '_aligned'];
cytodir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony4/' basename '_GFP/' basename '_aligned'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony4'];

%load cell wall images
cd(cwdir)
load([basename '_colony4_preShock'], 'im');
imCpre=im;
load([basename '_colony4_postShock'], 'im');
imCpost=im;

%load boundary coordinates
load([basename '_colony4_preShock'], 'B');
preB=B;
load([basename '_colony4_postShock'], 'B');
postB=B;

%load pixel locations
load([basename '_colony4_preShock'], 'pixels');
prePixels=pixels;
load([basename '_colony4_postShock'], 'pixels');
postPixels=pixels;

%load phase images
cd(phasedir)
directory=dir('*.tif');
imagename=directory(1).name;
imPpre=imread(imagename);
imagename=directory(2).name;
imPpost=imread(imagename);

%load cytoplasm images
cd(cytodir)
directory=dir('*.tif');
imagename=directory(1).name;
imGpre=imread(imagename);
imagename=directory(2).name;
imGpost=imread(imagename);

% %% Troubleshooting
% for n=1:height(preB)
%     n
%     figure
%     imshow(imCpre)
%     hold on
% 
%     if isempty(preB{n,1})==0
%         plot(preB{n,1}(:,1),preB{n,1}(:,2),'-r')
%     end
% 
%       pause
%       close all
% end
% 
% for n=1:height(postB)
%     n
%     figure
%     imshow(imCpost)
%     hold on
% 
%     if isempty(postB{n,1})==0
%         plot(postB{n,1}(:,1),postB{n,1}(:,2),'-r')
%     end
% 
%       pause
%       close all
% end
%% Generate a binary image of the occulded/stained cells
y=1:size(imCpre,1);
x=1:size(imCpre,2);
[X,Y] = meshgrid(x,y);

mx=max(x);
my=max(y);

pre_in = zeros(size(imCpre));
pre_on = zeros(size(imCpre));
post_in = zeros(size(imCpre));
post_on = zeros(size(imCpre));

for i=1:height(preB)
    pre_xv=preB{i,1}(:, 1);
    pre_yv=preB{i,1}(:, 2);
    
    [in, on]=inpolygon(X,Y,pre_xv, pre_yv);
    pre_in = pre_in + in;
    pre_on = pre_on + on;
    
end

for i=1:height(postB)
    post_xv=postB{i,1}(:, 1);
    post_yv=postB{i,1}(:, 2);
    
    [in, on]=inpolygon(X,Y,post_xv, post_yv);
    post_in = post_in + in;
    post_on = post_on + on;
    
end

[preL,pre_bw]=bwboundaries(pre_in,4,'noholes');
[postL,post_bw]=bwboundaries(post_in,4,'noholes');
preP=struct2cell(regionprops(pre_bw,'PixelIdxList'));
postP=struct2cell(regionprops(post_bw,'PixelIdxList'));

%calculate the area of pre- and post-shocked cells
pre_area=cellfun(@height, preP)';
post_area=cellfun(@height, postP)';

%% Generate a binary image of the GFP labelled cytoplasm
preMat=zeros(size(imCpre));
postMat=zeros(size(imCpost));

preMat_area=nan(height(pre_area),1);
postMat_area=nan(height(post_area),1);

minA=65535;
for i=1:height(pre_area)
    if min(min(imGpre(preP{i})))<minA
        minA=min(min(imGpre(preP{i})));
    end
end
for i=1:height(pre_area)
    idx=find(imGpre(preP{i})>=minA);
    preMat(preP{i}(idx,1))=1;
    preMat_area(i)=length(idx);
end

minB=65535;
for i=1:height(post_area)
    if min(min(imGpost(postP{i})))<minB
        minB=min(min(imGpost(postP{i})));
    end
end
for i=1:height(post_area)
    idx=find(imGpost(postP{i})>=minB);
    postMat(postP{i}(idx,1))=1;
    postMat_area(i)=length(idx);
end

%% Calculate the percent plasmolysis
pre_plasmolysis=pre_area-preMat_area; 
post_plasmolysis=post_area-postMat_area; 

pre_perc=(pre_plasmolysis./pre_area)*100;
post_perc=(post_plasmolysis./post_area)*100;

figure(1)
h1 = scatter(zeros(height(pre_perc), 1), pre_perc,'o');
hold on
h2 = scatter(ones(height(post_perc), 1), post_perc,'x');
title('Pre-shock Plasmolysis vs Post-Shock Plasmolysis')
ylabel('Percent Plasmolysis')
xticklabels({'pre-shock', 'post-shock'})
xticks([0 1])
xlim([-0.2 1.2])

%% Troubleshooting
%for n=1:height(preB)
%     n
figure(1)
imshow(imGpre)
hold on
    
for n=1:height(preL)
    if isempty(preL{n,1})==0
        plot(preL{n,1}(:,2),preL{n,1}(:,1),'-r')
    end
end

%       pause
%       close all
%end

% for n=1:height(postB)
%     n
figure(2)
imshow(imGpost)
hold on

for n=1:height(postL)
    if isempty(postL{n,1})==0
        plot(postL{n,1}(:,2),postL{n,1}(:,1),'-r')
    end
end

%       pause
%       close all
% end

%% save data
cd(savedir);
save([basename '_colony4_PT'])
