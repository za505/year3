%PlasTrack.m
%Tracks percent plasmolyis from GFP/CY5/RFP image stacks.
%Dylan Fitzmaurice
%edit: Zarina Akbary

clear 
close all

%% load the WallTrack.m data
basename='12182021_Exp2';
cwdir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12182021_analysis/' basename '/' basename '_colony3/' basename '_647/' basename '_figures'];
phasedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12182021_analysis/' basename '/' basename '_colony3/' basename '_phase/' basename '_aligned'];
cytodir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12182021_analysis/' basename '/' basename '_colony3/' basename '_GFP/' basename '_aligned'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12182021_analysis/' basename '/' basename '_colony3'];

%load cell wall images, midline variable, and delind
cd(cwdir)
load([basename '_colony3_preShock'], 'im', 'mline', 'delind');
    mline(delind, :)=[];
    preMid=mline;
imCpre=im;
load([basename '_colony3_postShock'], 'im', 'mline', 'delind');
    mline(delind, :)=[];
    postMid=mline;
imCpost=im;

%load boundary coordinates
load([basename '_colony3_preShock'], 'B');
preB=B;
load([basename '_colony3_postShock'], 'B');
postB=B;

%load pixel locations
load([basename '_colony3_preShock'], 'pixels');
prePixels=pixels;
load([basename '_colony3_postShock'], 'pixels');
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

%% Troubleshooting
%how well do the boundaries and pixels overelap with the pre-shock frame?
% figure(1)
% imshow(imGpre)
% hold on
% for n=1:height(preB)
%     n
%     if isempty(preB{n,1})==0
%         plot(preB{n,1}(:,1),preB{n,1}(:,2),'-r')
%     end
% end
% 
% %how well do the boundaries and pixels overelap with the post-shock frame?
% figure(2)
% imshow(imGpost)
% hold on
% for n=1:height(postB)
%     n
%     if isempty(postB{n,1})==0
%         plot(postB{n,1}(:,1),postB{n,1}(:,2),'-r')
%     end
% end
% 
% pause, close all

% overlay midline, poles, and midpoint on GFP image
% figure
% imshow(imGpre)
% hold on
% for n=1:height(preMid)
%     n
%     midx=round(length(preMid{n,1})/2);
%     plot(round(preMid{n,1}(midx,1)),round(preMid{n,1}(midx,2)),'b*')
%     plot(round(preMid{n,1}(1,1)),round(preMid{n,1}(1,2)),'r*')
%     plot(round(preMid{n,1}(end,1)),round(preMid{n,1}(end,2)),'r*')
% end
% 
% figure
% imshow(imGpost)
% hold on
% for n=1:height(postMid)
%     n
%     midx=round(length(postMid{n,1})/2);
%     plot(round(postMid{n,1}(midx,1)),round(postMid{n,1}(midx,2)),'b*')
%     plot(round(postMid{n,1}(1,1)),round(postMid{n,1}(1,2)),'r*')
%     plot(round(postMid{n,1}(end,1)),round(postMid{n,1}(end,2)),'r*')
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

% %calculate the area of pre- and post-shocked cells
% pre_area=nan(height(prePixels),1);
% post_area=nan(height(postPixels),1);
% 
% for i=1:height(prePixels)
%     pre_area(i)=sum(pre_in(prePixels{i})==1);
% end
% 
% for i=1:height(postPixels)
%     post_area(i)=sum(post_in(postPixels{i})==1);
% end
% 
% %calculate the area of pre- and post-shocked cells
% pre_area2=cellfun(@height, prePixels);
% post_area2=cellfun(@height, postPixels);

%% how sufficient is the overlap between the binary image and the CY5 and GFP images?
% figure(3)
% imshowpair(pre_bw, imCpre)
% 
% figure(4)
% imshowpair(pre_bw, imGpre)
% 
% figure(5)
% imshowpair(post_bw, imCpost)
% 
% figure(6)
% imshowpair(post_bw, imGpost)
% 
% pause, close all
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
    pre_perc(delind_pre)=[];
    prePixels(delind_pre)=[];
    preMid(delind_pre)=[];
    preP(delind_pre)=[];
    prebp(delind_pre)=[];
end

if isempty(delind_post)==0
    post_perc(delind_post)=[];
    postPixels(delind_post)=[];
    postMid(delind_post)=[];
    postP(delind_post)=[];
    postbp(delind_post)=[];
end

%plot percent plasmolysis
cd(savedir);
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
% saveas(gcf, [basename '_colony3_plasmolysis.png'])
% saveas(gcf, [basename '_colony3_plasmolysis.fig'])

%% calculate distance from midpoint to pole and midpoint to each pixel in the plasmolysis bay
pre_distTotal=[];
post_distTotal=[];

pre_distBay={};
post_distBay={};

pre_percDist={};
post_percDist={};

% pre-shock
for n=1:height(preMid)
    n
    midx=round(length(preMid{n,1})/2);
      pre_distTotal(n,1)=midx;
      for p=1:height(prebp{n,1})
          pre_distBay{n,1}(p,1)=distCalc(round(preMid{n,1}(midx,:)), prebp{n}(p,:));
          pre_percDist{n,1}(p,1)=(pre_distBay{n,1}(p,1)/pre_distTotal(n,1))*50;
      end
end

% post-shock
for n=1:height(postMid)
    n
    midx=round(length(postMid{n,1})/2);
      post_distTotal(n,1)=midx;
      for p=1:height(postbp{n,1})
          post_distBay{n,1}(p,1)=distCalc(round(postMid{n,1}(midx,:)), postbp{n}(p,:));
          post_percDist{n,1}(p,1)=(post_distBay{n,1}(p,1)/post_distTotal(n,1))*50;
      end
end

%sort the distances into categories
preR=plasRegion(pre_percDist);
postR=plasRegion(post_percDist);

%% Troubleshooting
for n=1:height(preB)
    n
    figure(1)
    imshow(imGpre)
    hold on

    if ~isempty(pre_percDist{n,1})
        plot(preB{n,1}(:,1),preB{n,1}(:,2),'-r')

        figure(2)
        histogram(preR{n,1})
    end

      pause
      close all
end
 
% for n=12 %1:height(postB)
%     n
%     figure
%     imshow(imGpost)
%     hold on
% 
%     if isempty(postB{n,1})==0
%         plot(postB{n,1}(:,2),postB{n,1}(:,1),'-r')
%     end
% 
%       pause
%       close all
% end

% figure(2)
% imshowpair(pre_in, preBay)
% 
% figure(3)
% imshowpair(post_in, postBay)
% imshowpair(imGpre, preBay)
% imshowpair(imGpost, postBay)

% figure
% histogram(cell2mat(pre_percDist))
% xlim([0 50])
% xlabel('Percent Cell Length')
% ylabel('Number of Pixels Pre-shock')
% 
% figure
% histogram(cell2mat(post_percDist))
% xlim([0 50])
% xlabel('Percent Cell Length')
% ylabel('Number of Pixels Post-shock')
% 
% figure
% histogram(cellfun(@median, pre_percDist), height(pre_percDist))
% xlim([0 50])
% xlabel('Median Percent Cell Length')
% ylabel('Number of Cell Pre-shock')
% 
% figure
% histogram(cellfun(@median, post_percDist), height(post_percDist))
% xlim([0 50])
% xlabel('Median Percent Cell Length')
% ylabel('Number of Cell Post-shock')
%% save data
cd(savedir);
save([basename '_colony3_PT'])

cd('/Users/zarina/Documents/MATLAB/MatlabReady/Msmegmatis_analysis');
save([basename '_colony3_PT'])

%% Functions
function [midDist]=distCalc(mat1, mat2) %for mat1 (preMid), x is the first column. for mat2 (prebp), x is the second column
    x=(mat2(2)-mat1(1))^2;
    y=(mat2(1)-mat1(2))^2;
    midDist=sqrt(x+y);    
end

function [plasm_region]=plasRegion(percDist)
    idx=cellfun(@isempty, percDist);
    idx=setdiff(1:height(percDist), idx);
    
    percDist=percDist(idx, 1);
    plasm_region=cell(size(percDist));
    
    for n=1:height(percDist)
        for p=1:height(percDist{n,1})
            pdist=percDist{n,1}(p,1);
            if pdist <12.5
               plasm_region{n,1}(p,1)="polar";
            elseif pdist >=12.5 & pdist <25
                plasm_region{n,1}(p,1)="subpolar";
            elseif pdist >=25
                plasm_region{n,1}(p,1)="mid-cell";
            end
        end
    end
    
    %y = discretize(x,[0 .25 .75
    %1],'categorical',{'small','medium','large'}); maybe use this next
    %time?
    plasm_region=cellfun(@categorical, plasm_region, 'UniformOutput', false);
end