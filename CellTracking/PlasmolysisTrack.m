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
frameSwitch=136; %the last frame before hyperosmotic shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(savedir)
load([basename '_BTphase.mat'],'T','labels','im','lcents','cellnum','tstamp','l','lscale', 'directory', 'dirname', 'pixels', 'bw');

%First, let's calculate the number of cells in our preshock image
cellnumber=max(max(labels(:,:,frameSwitch)));
stats=regionprops(labels(:, :, frameSwitch),'Area','Centroid','PixelIdxList');
pixelIntensity=cell(cellnumber, 4);
r1c1=cell(cellnumber,1);

%now, let's get the pixel intensity, centroid, and area for each of those
%preshock cells
for n=1:cellnumber
    
    pixelIntensity{n,1}=stats(n).PixelIdxList;
    
    [r1,c1] = find(bwlabel(labels(:,:,frameSwtich))==n);   %Find coordinates of cells in frame T
    r1c1{n,1}=[r1 c1]; % concatenate XY coords frame T
    
end


flag=1;
%let's go through the hyperosmotic shock images
for t=frameSwitch+1:frameSwitch+3
    
    %Load image
    imagename=directory(t).name;
    im=imread(imagename);
    
    for n=1:cellnumber
        pixelIntensity{n,1}=
    end
    
    flag=flag+1;
    
end

%Load image
    imagename=directory(t).name;

    im=imread(imagename);
    [imM,imN]=size(im);
    
    %De-speckle image
    im=medfilt2(im);
    
    %Normalize images
    ppix=0.5;
    im=norm16bit(im,ppix);
    
    %Enhance contrast
    imc=imcomplement(im);
    
%% Part 2: Lineage Tracking

% %view binary movie
% for i=flip(1:T)
% imshow(labels(:,:,i))
% pause(0.1)
% clf
% end
% %

%Get coordinates of cells
for i = flip(1:T)
numofcells(i)=max(max(labels(:,:,i)));% find number of cells
mxnumofcells=max(numofcells);
    for j=1:mxnumofcells
        [r1,c1] = find(bwlabel(labels(:,:,i))==j);   %Find coordinates of cells in frame T

        r1c1{j,i}=[r1 c1]; % concatenate XY coords frame T 
        
    end
end

for i=1:T
    for j=1:mxnumofcells
        if isempty(r1c1{j,i})==0 
            cell_id_na(j,i)=j; % cell id not aligned/ not lineage tracked
        end
    end
end

for i = flip(1:T)
    for j=1:mxnumofcells
        if isempty(r1c1{j,i})==1
            r1c1{j,i} = [NaN NaN];
        end
    end
end

%Align cells from frame to frame, correct r1c1 to track lineage
for i=1:T %should it be 1:T?
    for j=1:mxnumofcells
        cell_id(j,i)=j; %makes cell id vector
        r2c2{j,i}=r1c1{j,i};%Seed r2c2
    end
end

%Tracks cells backward through time,  
for i = flip(2:T)
    for j=1:mxnumofcells
        if isnan(r2c2{j,i})==1 %if NaN, find previous nonNaN id matrix for tracking
            index=cellfun('length', r2c2(j,:));
            nindex=index>2; %I think this is supposed to be 2
            ncols=min(find(nindex));
            Lia1=ismember(r2c2{j,ncols},r1c1{j,i-1},'rows');
            if 100*sum(Lia1)/length(r1c1{j,i-1})>50==1
                cell_id(j,i-1)=j;
                r2c2{j,i-1}=r1c1{j,i-1};
            else
                cell_id(j,i-1)=NaN;
                r2c2{j,i-1}=NaN;
            end
        else
            nidx=cellfun('length', r2c2(j,:)); %Normal tracking backward 
            nnidx=nidx>1;
            ncols=min(find(nnidx));
            Lia1=ismember(r1c1{j,i-1},r2c2{j,ncols},'rows');
                if 100*sum(Lia1)/length(r2c2{j,ncols})>50==1
                            cell_id(j,i-1)=j;
                            r2c2{j,i-1}=r1c1{j,i-1};
                    else
                        r2c2{j,i-1}=NaN;
                end
        end
     end
end

% Fill NaNs tracking backward
for i = flip(1:T)
for j=1:mxnumofcells
    for k=1:mxnumofcells
        if isnan(r2c2{j,i})==1
            index=cellfun('length', r2c2(j,:));
            nindex=index==1;
            ncols=max(find(nindex));
            Lia1=ismember(r2c2{j,ncols+1},r1c1{k,i},'rows');
            if 100*sum(Lia1)/length(r1c1{k,ncols+1})>50==1
                cell_id(j,i)=k;
                r2c2{j,i}=r1c1{k,i};
            else
                %cell_id(j,i-1)=NaN;
                %r2c2{j,i-1}=NaN;
                end
        end
     end
end
end


