%PlasTrack.m
%Tracks percent plasmolyis from GFP/CY5/RFP image stacks.
%Dylan Fitzmaurice
%edit: Zarina Akbary

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
%% Generate a binary image of the occulded/stained cells
pre_bw = zeros(size(imCpre));
post_bw = zeros(size(imCpre));

%generate binary image of pre-shocked cells
for i=1:height(prePixels)
    pre_bw(prePixels{i})=1;
end

%generate binary image of post-shocked cells
for i=1:height(postPixels)
    post_bw(postPixels{i})=1;
end

%calculate the area of pre- and post-shocked cells
pre_area=cellfun(@height, prePixels);
post_area=cellfun(@height, postPixels);

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
preBay=zeros(size(pre_bw));
postBay=zeros(size(post_bw));

%generate binary image of pre-shocked cells
pre_level = graythresh(imGpre);
preMat = imbinarize(imGpre, pre_level);

%generate binary image of post-shocked cells
post_level = graythresh(imGpost);
postMat = imbinarize(imGpost, post_level);

%calculate area
preMat_area=nan(height(pre_area),1);
postMat_area=nan(height(post_area),1);

for i=1:height(prePixels)
    preMat_area(i)=sum(preMat(prePixels{i})==1);
    idx=find(preMat(prePixels{i})==0);
    preBay(prePixels{i}(idx))=1;
end

for i=1:height(postPixels)
    postMat_area(i)=sum(postMat(postPixels{i})==1);
    idx=find(postMat(postPixels{i})==0);
    postBay(postPixels{i}(idx))=1;
end

%% Calculate the percent plasmolysis
pre_plasmolysis=pre_area-preMat_area; 
post_plasmolysis=post_area-postMat_area; 

pre_perc=(pre_plasmolysis./pre_area)*100;
post_perc=(post_plasmolysis./post_area)*100;

cd(savedir);
figure(1)
h1 = scatter(zeros(height(pre_perc), 1), pre_perc,'k');
hold on
h2 = scatter(ones(height(post_perc), 1), post_perc,'k');
title('Pre-shock Plasmolysis vs Post-shock Plasmolysis')
ylabel('Percent Plasmolysis')
xticklabels({'Pre-shock', 'Post-shock'})
xticks([0 1])
xlim([-0.5 1.5])
saveas(gcf, [basename '_colony4_plasmolysis.png'])
saveas(gcf, [basename '_colony4_plasmolysis.fig'])

%% Troubleshooting
% for n=1:height(preB)
%     n
%     figure
%     imshow(imGpre)
%     hold on
% 
%     if isempty(preB{n,1})==0
%         plot(preB{n,1}(:,2),preB{n,1}(:,1),'-r')
%     end
% 
%       pause
%       close all
% end
% 
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

imshowpair(imGpre, preBay)
imshowpair(imGpost, postBay)
%% save data
cd(savedir);
save([basename '_colony4_PT'])