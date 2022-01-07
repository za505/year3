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

staining=2;

cd(cwdir)
load([basename '_colony4_preShock'], 'im');
imCpre=im;
load([basename '_colony4_postShock'], 'im');
imCpost=im;

load([basename '_colony4_preShock'], 'B');
preB=B;
load([basename '_colony4_postShock'], 'B');
postB=B;

cd(phasedir)
directory=dir('*.tif');
imagename=directory(1).name;
imPpre=imread(imagename);
imagename=directory(2).name;
imPpost=imread(imagename);

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
preP=regionprops(pre_bw,'Area','PixelIdxList');
postP=regionprops(post_bw,'Area','PixelIdxList');

if staining==1
    cw_cutoff=5400;
    mean1=0;
    for i=1:height(preP)  
        if mean(mean(imCpre(preP(i).PixelIdxList)))>mean1
            mean1=mean(mean(imCpre(preP(i).PixelIdxList)));
        end
        idx=find(imCpre(preP(i).PixelIdxList)<cw_cutoff);
        idx=setdiff(preP(i).PixelIdxList, preP(i).PixelIdxList(idx));
        pre_in(idx)=0;
    end

    mean2=0;
    for i=1:height(postP) 
        if mean(mean(imCpre(preP(i).PixelIdxList)))>mean2
            mean2=mean(mean(imCpre(preP(i).PixelIdxList)));
        end
        idx=find(imCpost(postP(i).PixelIdxList)<cw_cutoff);
        idx=setdiff(postP(i).PixelIdxList, postP(i).PixelIdxList(idx));
        post_in(idx)=0;
    end

    [preL,pre_bw]=bwboundaries(pre_in,4,'noholes');
    [postL,post_bw]=bwboundaries(post_in,4,'noholes');
    preP=regionprops(pre_bw,'Area','PixelIdxList');
    postP=regionprops(post_bw,'Area','PixelIdxList');

end

preP_area=struct2cell(preP);
preP_area=cell2mat(preP_area(1,:))';

postP_area=struct2cell(postP);
postP_area=cell2mat(postP_area(1,:))';
%% Generate a binary image of the GFP labelled cytoplasm
preMat=zeros(size(imCpre));
postMat=zeros(size(imCpost));

preMat_area=nan(height(preP),1);
postMat_area=nan(height(postP),1);

gfp_cutoff=8890;
meanA=65535;
for i=1:height(preP)
    if mean(mean(imGpre(preP(i).PixelIdxList)))<meanA
        meanA=mean(mean(imGpre(preP(i).PixelIdxList)));
    end
    idx=find(imGpre(preP(i).PixelIdxList)>gfp_cutoff);
%    preMat(preP(i).PixelIdxList(idx))=imGpre(preP(i).PixelIdxList(idx));
    preMat(preP(i).PixelIdxList(idx))=1;
    preMat_area(i)=length(idx);
end

meanB=65535;
for i=1:height(postP)
    if mean(mean(imGpre(postP(i).PixelIdxList)))<meanB
        meanB=mean(mean(imGpre(postP(i).PixelIdxList)));
    end
    idx=find(imGpost(postP(i).PixelIdxList)>gfp_cutoff);
%     postMat(postP(i).PixelIdxList(idx))=imGpost(postP(i).PixelIdxList(idx));
    postMat(postP(i).PixelIdxList(idx))=1;   
    postMat_area(i)=length(idx);
end

%% Calculate the percent plasmolysis
pre_plasmolysis=preP_area-preMat_area; 
post_plasmolysis=postP_area-postMat_area; 

pre_perc=(pre_plasmolysis./preP_area)*100;
post_perc=(post_plasmolysis./postP_area)*100;

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
for n=1:height(preL)
    n
    figure
    imshow(imGpre)
    hold on

    if isempty(preL{n,1})==0
        plot(preL{n,1}(:,2),preL{n,1}(:,1),'-r')
    end

      pause
      close all
end

for n=1:height(postL)
    n
    figure
    imshow(imGpost)
    hold on

    if isempty(postL{n,1})==0
        plot(postL{n,1}(:,2),postL{n,1}(:,1),'-r')
    end

      pause
      close all
end

%% save data
cd(savedir);
save([basename '_colony4_PT'])
