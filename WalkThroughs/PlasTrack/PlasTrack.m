%PlasTrack.m
%Tracks percent plasmolyis from phase/CY5 image stacks.
%Dylan Fitzmaurice
%edit: Zarina Akbary
%purpose: I want to compare the plasmolysis quantification b/t TADA/phase,
%TADA/GFP, and/or phase/647

clear 
close all

%% load the WallTrack.m data
basename='12162021_Exp1';
cwdir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_647/' basename '_figures'];
phasedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_phase/' basename '_aligned'];
cytodir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_GFP/' basename '_aligned'];

cd(cwdir)
load([basename '_colony1_preShock'], 'im');
imCpre=im;
load([basename '_colony1_postShock'], 'im');
imCpost=im;

load([basename '_colony1_preShock'], 'B');
preB=B;
load([basename '_colony1_postShock'], 'B');
postB=B;

load([basename '_colony1_preShock'], 'pixels');
prePixels=pixels;
load([basename '_colony1_postShock'], 'pixels');
postPixels=pixels;

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

% A=load([basename '_ompA_CY5_1_BT']);%load CY5 data 
% CY5_im3=A.im3;% CY5 image
% CY5_boun=A.boun;% CY5 bound
%%
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

% pre_sum=sum(sum(pre_in));
% post_sum=sum(sum(post_in));

preMat=zeros(size(imCpre));
postMat=zeros(size(imCpost));

for i=1:height(prePixels)
    preMat(prePixels{i, 1})=1;
end

for i=1:height(postPixels)
    postMat(postPixels{i, 1})=1;
end

figure(1)
imshowpair(pre_in, preMat)

figure(2)
imshowpair(post_in, postMat)

pause
% 
% frame=1;
% 
% figure(1)
% imshow(im3{1,frame})
% hold on
% plot(CY5_boun{1,frame}(:,1),CY5_boun{1,frame}(:,2),'-r')
% 
% test_im=im3{1,frame};
% ntest_im=ones(my,mx);
% for j=1:mx%x direction
%   for i=1:my%y direction
%     if test_im(i,j)>20000
%       ntest_im(i,j)=45000;
%     elseif test_im(i,j)<20000
%       ntest_im(i,j)=test_im(i,j);
%     else 
%       continue
%     end
%   end
% end
% ntest_im=uint16(ntest_im)
% ntest_im(1,1)=45000;
% ntest_im(my,1)=45000;
% ntest_im(1,mx)=45000;
% ntest_im(my,mx)=45000;
% 
% figure
% imshow(ntest_im)
% hold on
% plot(CY5_boun{1,frame}(:,1),CY5_boun{1,frame}(:,2),'-r')
% 
% % %filter
% % for sizeim3{1,40})
% 
% % for i=1:59
% % CY5_Z{1,i}=ones(size(CY5_boun{1,i},1),2);
% % CY5_Z{1,i}=CY5_Z{1,i}*65530;
% % end
% 
% for i=frame  
% A_xv=CY5_boun{1,i}(:,1);%total boun
% A_yv=CY5_boun{1,i}(:,2);%total boun
% [A_in,A_on] = inpolygon(X,Y,A_xv,A_yv);%total boun
% end
% total_pix=sum(sum(A_in));%total number of pixels within boundary
% 
% im_for_percent_plas=double(ntest_im)-A_in;
% 
% for j=1:mx%x direction
%   for i=1:my%y direction
%     if A_in(i,j)==0
%        im_for_percent_plas(i,j)=NaN;
%     end
%   end
% end
% 
% for j=1:mx%x direction
%   for i=1:my%y direction
%     if im_for_percent_plas(i,j)>=44999
%       im_for_percent_plas(i,j)=NaN;
%     else 
%       continue
%     end
%   end
% end
% 
% for j=1:mx%x direction
%   for i=1:my%y direction
%     if im_for_percent_plas(i,j)>-2
%       im_for_percent_plas(i,j)=1;
%     else 
%       continue
%     end
%   end
% end
% 
% plasm=sum(nansum(im_for_percent_plas));
% percent_plasmolysis=100-(plasm/total_pix)*100;

