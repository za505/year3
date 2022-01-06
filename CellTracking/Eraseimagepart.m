%Eraseimagepart.m
%Rico Rojas, updated 1/21/19
%This routine lets the user select portions of an image stack to erase (make background).   

%INSTRUCTIONS FOR USE:
%Save the image stack in a directory without any other .tif files in it.  When 
%you run the program, the final image in the image stack will open.  
%Select the regions you want to delete with the cursor and then press Enter.  
%The program writes over the original image stack, so if you want a backup stack, 
%save it in a separate location.

%INPUT:
%dirname:Directory in which the image stack is saved.

%Calls upon:
%norm16bit.m

clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='12162021_Exp1';
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_647/' basename '_preShock'];
%mirname=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/12162021_analysis/' basename '/' basename '_colony1/' basename '_phase/' basename '_aligned'];
split=0; %there is a break in the movie
%split1=13; %first frame for segment 2
%split2=; %last frame for segment 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curdir=cd;

cd(dirname);
directory=dir('*.tif');
T=length(directory);
path(dirname,path)
       
    
for t=1:T
    t
    
    %Pick the regions to erase     
    [rp1, rp2, count, im1]=pickRegion(dirname, directory, t);
    
    for n=1:count
        %Load image
        imagename=directory(t).name;
        
        im=imread(imagename);
        [imM,imN]=size(im);
        
        [imcounts,bins]=hist(double(nonzeros(im1)));
        [~,mpos]=max(imcounts);
        val=bins(mpos);
        %maxval=bins(end-3);
        %val=maxval;
        im(rp1(n,2):rp2(n,2),rp1(n,1):rp2(n,1))=val*ones(rp2(n,2)-rp1(n,2)+1,rp2(n,1)-rp1(n,1)+1);
        delete(imagename);
        imwrite(im,imagename);
    end
end

close all

%% Functions
function [rp1, rp2, count, im1]=pickRegion(dirname, directory, t)
    
    %load image
    cd(dirname)
    imagename=directory(t).name;
    im1=imread(imagename);
    [imM,imN]=size(im1);

    ppix=0.5;
    im2=norm16bit(im1,ppix);

    %figure, imshowpair(imA, imB, 'montage')
    figure,imshow(im1,stretchlim(im1)*65000)
    set(gcf,'Pointer','fullcross')
    hold on
    axis manual
    
    count=0;
    k=0;
    
    while k~=1
        count=count+1;
        k=waitforbuttonpress;
        point1=get(gca,'CurrentPoint');   
        finalRect=rbbox;                   
        %pause
        point2=get(gca,'CurrentPoint');    
        point1=point1(1,1:2);              
        point2=point2(1,1:2);
        point1(point1<1)=1;
        point2(point2<1)=1;
        
        if point1(2)>imM
            point1(2)=imM;
        end
        if point1(1)>imN
            point1(1)=imN;
        end
        if point2(2)>imM
            point2(2)=imM;
        end
        if point2(1)>imN
            point2(1)=imN;
        end
        p1=min(point1,point2);%Calculate locations
        p2=max(point1,point2);
        offset = abs(point1-point2);%And dimensions
        x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
        y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
        plot(x,y)

        rp1(count,:)=round(p1);
        rp2(count,:)=round(p2);
    end
 end
    
