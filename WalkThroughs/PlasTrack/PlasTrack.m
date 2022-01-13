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
load([basename '_colony3_preShock'], 'im');
imCpre=im;
load([basename '_colony3_postShock'], 'im');
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

%% Generate a binary image of the GFP labelled cytoplasm
preBay=zeros(size(imGpre));
postBay=zeros(size(imGpost));

%binarize pre-shocked GFP image
pre_level = graythresh(imGpre);
preMat = imbinarize(imGpre, pre_level);

%binarize post-shocked GFP image
post_level = graythresh(imGpost);
postMat = imbinarize(imGpost, post_level);

%generate filled GFP images
preFill=imfill(preMat, 'holes');
postFill=imfill(postMat, 'holes');

%generate pixel id list
preP=struct2cell(regionprops(preFill,'PixelIdxList'));
postP=struct2cell(regionprops(postFill,'PixelIdxList'));

%calculate the area of pre- and post-shocked cells
pre_area=cellfun(@height, preP)';
post_area=cellfun(@height, postP)';

%generate binary image of pre-shocked cells
for i=1:height(preP)
    preMat_area(i)=sum(preMat(preP{i})==1);
    idx=find(preMat(preP{i})==0);
    preBay(preP{i}(idx))=1;
    
    [x,y]=ind2sub(size(im), preP{i}(idx));
    prebp{i,1}=[x,y];
end

for i=1:height(postP)
    postMat_area(i)=sum(postMat(postP{i})==1);
    idx=find(postMat(postP{i})==0);
    postBay(postP{i}(idx))=1;
    
    [x,y]=ind2sub(size(im), postP{i}(idx));
    postbp{i,1}=[x,y];
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