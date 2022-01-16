%PlasTrack.m
%Tracks percent plasmolyis from GFP/TADA/TADA image stacks.
%Dylan Fitzmaurice
%edit: Zarina Akbary

clear 
close all

%% load the WallTrack.m data
basename='12162021_Exp4';
cwdir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony2/' basename '_TADA/' basename '_figures'];
phasedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony2/' basename '_phase/' basename '_aligned'];
cytodir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony2/' basename '_GFP/' basename '_aligned'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony2'];

%load cell wall images, midline variable, and delind
cd(cwdir)
load([basename '_colony2_preShock'], 'im', 'mline', 'delind');
    mline(delind, :)=[];
    preMid=mline;
imCpre=im;
load([basename '_colony2_postShock'], 'im', 'mline', 'delind');
    mline(delind, :)=[];
    postMid=mline;
imCpost=im;

%load boundary coordinates
load([basename '_colony2_preShock'], 'B');
preB=B;
load([basename '_colony2_postShock'], 'B');
postB=B;

%load pixel locations
load([basename '_colony2_preShock'], 'pixels');
prePixels=pixels;
load([basename '_colony2_postShock'], 'pixels');
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

% %calculate the area of pre- and post-shocked cells
% pre_area2=cellfun(@height, prePixels);
% post_area2=cellfun(@height, postPixels);

%% Generate a binary image of the GFP labelled cytoplasm
preBay=zeros(size(pre_in));
postBay=zeros(size(post_in));

preCyto=zeros(size(pre_in));
postCyto=zeros(size(post_in));

%generate binary image of pre-shocked cells
pre_level = graythresh(imGpre);
preMat = imbinarize(imGpre, pre_level);

%generate binary image of post-shocked cells
post_level = graythresh(imGpost);
postMat = imbinarize(imGpost, post_level);

%calculate area
preMat_area=nan(height(pre_area),1);
postMat_area=nan(height(post_area),1);

prebp={};
postbp={};

for i=1:height(prePixels)
    preMat_area(i)=sum(preMat(prePixels{i})==1);
    idx=find(preMat(prePixels{i})==0);
    preBay(prePixels{i}(idx))=1;
    
    [x,y]=ind2sub(size(im), prePixels{i}(idx));
    prebp{i,1}=[x,y];
    
    idx=find(preMat(prePixels{i})==1);
    preCyto(prePixels{i}(idx))=1;
end

for i=1:height(postPixels)
    postMat_area(i)=sum(postMat(postPixels{i})==1);
    idx=find(postMat(postPixels{i})==0);
    postBay(postPixels{i}(idx))=1;
    
    [x,y]=ind2sub(size(im), postPixels{i}(idx));
    postbp{i,1}=[x,y];
    
    idx=find(postMat(postPixels{i})==1);
    postCyto(postPixels{i}(idx))=1;
end

%% Calculate the percent plasmolysis
pre_plasmolysis=pre_area-preMat_area; 
post_plasmolysis=post_area-postMat_area; 

pre_perc=(pre_plasmolysis./pre_area)*100;
post_perc=(post_plasmolysis./post_area)*100;

delind_pre=[];
delind_post=[];

%remove cells that moved between frames
for i=1:height(pre_perc)
    if pre_perc(i)==100
        delind_pre=[delind_pre, i];
    end
end

for i=1:height(post_perc)
    if post_perc(i)==100
        delind_post=[delind_post, i];
    end
end

if isempty(delind_pre)==0
    preB(delind_pre)=[];
    pre_perc(delind_pre)=[];
    prePixels(delind_pre)=[];
    preMid(delind_pre)=[];
    preP(delind_pre)=[];
    prebp(delind_pre)=[];
end

if isempty(delind_post)==0
    postB(delind_pre)=[];
    post_perc(delind_post)=[];
    postPixels(delind_post)=[];
    postMid(delind_post)=[];
    postP(delind_post)=[];
    postbp(delind_post)=[];
end

%plot percent plasmolysis
%cd(savedir);
% figure(1)
% h1 = scatter(zeros(height(pre_perc), 1), pre_perc,'k');
% hold on
% h2 = scatter(ones(height(post_perc), 1), post_perc,'k');
% title('Pre-shock Plasmolysis vs Post-shock Plasmolysis')
% ylabel('Percent Plasmolysis')
% xticklabels({'Pre-shock', 'Post-shock'})
% xticks([0 1])
% xlim([-0.5 1.5])
% ylim([0 100])
% saveas(gcf, [basename '_colony2_plasmolysis.png'])
% saveas(gcf, [basename '_colony2_plasmolysis.fig'])

%% Generate images where mid-cell point, poles, and subpolar region are all marked
% cd(cytodir)
% cd('../')
% cd(['./' basename '_figures'])

%1=polar
%2=subpolar
%3=midcell

preCount = cell(height(preMid),1);
postCount = cell(height(postMid),1);

% overlay midline, poles, and midpoint on GFP image
%imname1=[basename '_colony2_preShock.tif'];
for n=1:height(preMid)
%     figure
%     imshow(imGpre)
%     hold on
%     
%     midx=round(length(preMid{n,1})/2);
%     qidx=round(length(preMid{n,1})/8);
%     plot(round(preB{n,1}(:,1)),round(preB{n,1}(:,2)),'-w')    
%     plot(round(preMid{n,1}(midx,1)),round(preMid{n,1}(midx,2)),'b.', 'MarkerSize', 12)
%     plot(round(preMid{n,1}(1,1)),round(preMid{n,1}(1,2)),'r.', 'MarkerSize', 12)
%     plot(round(preMid{n,1}(end,1)),round(preMid{n,1}(end,2)),'r.', 'MarkerSize', 12)
%     plot(round(preMid{n,1}(end-qidx,1)),round(preMid{n,1}(end-qidx,2)),'r.', 'MarkerSize', 12)    
%     plot(round(preMid{n,1}(qidx,1)),round(preMid{n,1}(qidx,2)),'r.', 'MarkerSize', 12)
%     plot(round(preMid{n,1}(midx-qidx,1)),round(preMid{n,1}(midx-qidx,2)),'b.', 'MarkerSize', 12) 
%     plot(round(preMid{n,1}(midx+qidx,1)),round(preMid{n,1}(midx+qidx,2)),'b.', 'MarkerSize', 12) 
%     
    preCount{n,1} = 0; %input('How many plasmolysis bays did you count for this cell? What kind? ');
    
%    close
end
%saveas(gcf, imname1)
message=['end of pre-shock analysis']
n=height(postMid)
%imname2=[basename '_colony2_postShock.tif'];
for n=1:height(postMid)
    figure
    imshow(imGpost)
    hold on

    midx=round(length(postMid{n,1})/2);
    qidx=round(length(postMid{n,1})/8);
    plot(round(postB{n,1}(:,1)),round(postB{n,1}(:,2)),'-w')    
    plot(round(postMid{n,1}(midx,1)),round(postMid{n,1}(midx,2)),'b.', 'MarkerSize', 12)
    plot(round(postMid{n,1}(1,1)),round(postMid{n,1}(1,2)),'r.', 'MarkerSize', 12)
    plot(round(postMid{n,1}(end,1)),round(postMid{n,1}(end,2)),'r.', 'MarkerSize', 12)
    plot(round(postMid{n,1}(end-qidx,1)),round(postMid{n,1}(end-qidx,2)),'r.', 'MarkerSize', 12)    
    plot(round(postMid{n,1}(qidx,1)),round(postMid{n,1}(qidx,2)),'r.', 'MarkerSize', 12)
    plot(round(postMid{n,1}(midx-qidx,1)),round(postMid{n,1}(midx-qidx,2)),'b.', 'MarkerSize', 12) 
    plot(round(postMid{n,1}(midx+qidx,1)),round(postMid{n,1}(midx+qidx,2)),'b.', 'MarkerSize', 12)   
    
    pause, close
end

for n=1:height(postMid)
        postCount{n,1} = input('How many plasmolysis bays did you count for this cell? What kind? ');
end
%saveas(gcf, imname2)

% %load recently saved images
% imNpre=imread(imname1);
% imNpost=imread(imname2);
% 
% for i=1:height(preB)
%     i
%     x1=round(min(preB{n,1}(:,1)));
%     y1=round(max(preB{n,1}(:,2)));
%     x2=round(max(preB{n,1}(:,1)));
%     y2=round(min(preB{n,1}(:,2)));
%     
%     imC = imcrop(imNpre, [x1 y1 x2 y2]);
% 
%     figure
%     imshow(imC, [])
%     
%     preCount{i,1} = input('How many plasmolysis bays did you count for this cell?');
%     
%     close
% end
% 
% 
% for i=1:height(postB)
%     i
%     x1=round(min(postB{n,1}(:,1)));
%     y1=round(max(postB{n,1}(:,2)));
%     x2=round(max(postB{n,1}(:,1)));
%     y2=round(min(postB{n,1}(:,2)));
%     
%     imC = imcrop(imNpost, [x1 y1 x2 y2]);
% 
%     figure
%     imshow(imC, [])
%     
%     postCount{i,1} = input('How many plasmolysis bays did you count for this cell?');
%     
%     close
% end

%% save data
cd(savedir);
save([basename '_colony2_PT'])

cd('/Users/zarina/Documents/MATLAB/MatlabReady/Msmegmatis_analysis/01142022');
save([basename '_colony2_PT'])

% %% Functions
% function [plasm_region]=plasRegion(percDist)
%     idx=cellfun(@isempty, percDist);
%     idx=setdiff(1:height(percDist), idx);
%     
%     percDist=percDist(idx, 1);
%     plasm_region=cell(size(percDist));
%     
%     for n=1:height(percDist)
%         for p=1:height(percDist{n,1})
%             pdist=percDist{n,1}(p,1);
%             if pdist <12.5
%                plasm_region{n,1}(p,1)="polar";
%             elseif pdist >=12.5 & pdist <25
%                 plasm_region{n,1}(p,1)="subpolar";
%             elseif pdist >=25
%                 plasm_region{n,1}(p,1)="mid-cell";
%             end
%         end
%     end
%     
%     %y = discretize(x,[0 .25 .75
%     %1],'categorical',{'small','medium','large'}); maybe use this next
%     %time?
%     plasm_region=cellfun(@categorical, plasm_region, 'UniformOutput', false);
% end