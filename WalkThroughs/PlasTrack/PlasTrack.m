%PlasTrack.m
%Tracks percent plasmolyis from phase/CY5 image stacks.
%Dylan Fitzmaurice

clear 
close all

%%
basename='FC1_061221_P2';
load([basename '_ompA_BF_1' '_BT'])%load BF data 

A=load([basename '_ompA_CY5_1_BT']);%load CY5 data 
CY5_im3=A.im3;% CY5 image
CY5_boun=A.boun;% CY5 bound
%%

y=1:size(im3{1,1},1)%im3 = BF image
x=1:size(im3{1,1},2)
[X,Y] = meshgrid(x,y);

mx=max(x);
my=max(y);

% A=load([basename ‘_ompA_CY5_1_BT’])% Need to load CY5 data like 
% CY5_im3=A.im3;% CY5 image
% CY5_boun=A.boun;% CY5 bound

frame=1;

figure(1)
imshow(im3{1,frame})
hold on
plot(CY5_boun{1,frame}(:,1),CY5_boun{1,frame}(:,2),'-r')

test_im=im3{1,frame};
ntest_im=ones(my,mx);
for j=1:mx%x direction
  for i=1:my%y direction
    if test_im(i,j)>20000
      ntest_im(i,j)=45000;
    elseif test_im(i,j)<20000
      ntest_im(i,j)=test_im(i,j);
    else 
      continue
    end
  end
end
ntest_im=uint16(ntest_im)
ntest_im(1,1)=45000;
ntest_im(my,1)=45000;
ntest_im(1,mx)=45000;
ntest_im(my,mx)=45000;

figure
imshow(ntest_im)
hold on
plot(CY5_boun{1,frame}(:,1),CY5_boun{1,frame}(:,2),'-r')

% %filter
% for sizeim3{1,40})

% for i=1:59
% CY5_Z{1,i}=ones(size(CY5_boun{1,i},1),2);
% CY5_Z{1,i}=CY5_Z{1,i}*65530;
% end

for i=frame  
A_xv=CY5_boun{1,i}(:,1);%total boun
A_yv=CY5_boun{1,i}(:,2);%total boun
[A_in,A_on] = inpolygon(X,Y,A_xv,A_yv);%total boun
end
total_pix=sum(sum(A_in));%total number of pixels within boundary

im_for_percent_plas=double(ntest_im)-A_in;

for j=1:mx%x direction
  for i=1:my%y direction
    if A_in(i,j)==0
       im_for_percent_plas(i,j)=NaN;
    end
  end
end

for j=1:mx%x direction
  for i=1:my%y direction
    if im_for_percent_plas(i,j)>=44999
      im_for_percent_plas(i,j)=NaN;
    else 
      continue
    end
  end
end

for j=1:mx%x direction
  for i=1:my%y direction
    if im_for_percent_plas(i,j)>-2
      im_for_percent_plas(i,j)=1;
    else 
      continue
    end
  end
end

plasm=sum(nansum(im_for_percent_plas));
percent_plasmolysis=100-(plasm/total_pix)*100;

