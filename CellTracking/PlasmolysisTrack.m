%PlasmolysisTrack.m
%Tracks and quantifies size of plasmolysis bays in bacteria
%To be used after imagealign.m, Eraseimagepart.m, & BacTrack.m. Based on
%LinTrack.m
%Zarina Akbary, last updated 04/25/2021
clear 
close all

%input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basename='03082021_Exp3_colony3';
savedir=['/Users/zarina/Downloads/NYU/Year2_2021_Spring/PlasmolysisTrack_test/' basename '_phase/'  basename '_figures'];%Directory to save the output .mat file to.
frameSwitch=137; %the last frame before hyperosmotic shock
frameShock=140;
nframe=3;
vis=1; %to visualize the boundaries of the pre-shock cells plotted in post-shock frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(savedir)
load([basename '_BTphase.mat'],'T','labels','labels2','directory', 'dirname');

%Let's calculate the number of cells in our pre-shock frame
cellnumber=max(max(labels(:,:,frameSwitch)));

%Let's pull some stats from this image
bw=labels2(:,:,frameSwitch);
stats=regionprops(bw,'Area','Centroid','PixelList'); 

%pre-allocate variables
pixelIntensity=cell(cellnumber,nframe+2);
pixelMean=zeros(cellnumber, nframe+2);
pixIdx=cell(cellnumber,nframe+1);

%now let's get the row and column coordinates from the pre-shock frame
for n=1:cellnumber
    
    r1=stats(n).PixelList(:, 1);   %Find coordinates of cells in frame T
    c1=stats(n).PixelList(:, 2);
    pixelID{n,1}=[r1 c1]; % concatenate XY coords frame T 
    
    colID{n}=c1;
    rowID{n}=r1;
    
end

%Load pre-shock image
imagename=directory(frameSwitch).name;
im=imread(imagename);

%check that the coordinates are correct
figure
imshow(im)
hold on
for n=1:cellnumber
   plot(pixelID{n,1}(:,1),pixelID{n,1}(:,2),'-r')
end
pause
   
for n=1:cellnumber
   
    %find the intensity for each pixel in each identified cell
    for p=1:height(rowID{n})
        pixelIntensity{n,1}(p)=im(rowID{n}(p),colID{n}(p));
    end
    
    %calculate the mean
    pixelMean(n, 1)=mean(pixelIntensity{n,1});
    
end


%now, let's get the pixel intensities and pixel mean for post-shock cells
%in the same coordinates as the pre-shock cells
flag=2;
 
for t=frameShock:frameShock+nframe
    
    %Load image
    imagename=directory(t).name;
    im=imread(imagename);

    %check to see the coordinates are correct
    figure
    imshow(im)
    hold on
    for n=1:cellnumber
       plot(pixelID{n,1}(:,1),pixelID{n,1}(:,2),'-r')
    end
    pause

    for n=1:cellnumber
        
        for p=1:height(rowID{n})
            pixelIntensity{n,flag}(p)=im(rowID{n}(p),colID{n}(p));
        end
        
        pixelMean(n, flag)=mean(pixelIntensity{n,flag});
        
        flag2=1;
        for p=1:height(rowID{n})
            if im(rowID{n}(p),colID{n}(p)) > pixelMean(n, 1)
                pixIdx{n,flag-1}(flag2, :)=[rowID{n}(p),colID{n}(p)];
                flag2=flag2+1;
            end
        end
        
        %check to see the coordinates are correct
        figure
        imshow(im)
        hold on
        %for n=1:cellnumber
           plot(pixIdx{n,1}(:,1),pixIdx{n,1}(:,2),'-r')
        %end
        pause

    end
    
    
    flag=flag+1;
    
end

close all
