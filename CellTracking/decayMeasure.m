%decayMeasure.m
%Zarina Akbary, updated 02/09/22
%Calculates changes in fluor. intensity. Incorporates BTfluo.m code.

clear, close all

%INSTRUCTIONS FOR USE:
%run BacTrack.m first

%INPUT
%basename: experiments of interest
%dirname: where .mat files are stored
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%icell_intensity = cell x time matrix of fluoresc. intensities
%bg_intensity = 1 x time matrix of mean background intensities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='05062022_Exp1';%Name of the image stack, used to save file.

dirname=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/05062022_analysis/'  basename '_phase/' basename '_figures/']; %Directory that the BT.mat files is saved in
savedir=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/05062022_analysis/'  basename '_mNeonGreen/' basename '_figures']; %Directory to save the output .mat file to.
savedir2=['/Users/zarina/Downloads/NYU/Year3_2022_Spring/mNeonGreen_analysis/raw/' ]; %Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2022_Spring/05062022_analysis/'  basename '_mNeonGreen/' basename '_cell05']}; %Directory that the image stack is saved in.

recrunch=0;
troubleshoot=2;
LinTrack=0;
splitCell=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==2
    cd(savedir)
    load([basename '_cell05_dm.mat'])
    
    cd(channels{1});
        
    imagename=fluo_directory{1}(1).name;
    im=imread(imagename);
    sz = size(im);
    improxy = zeros(sz);
    
    for t=1:T
        t
        imagename=fluo_directory{1}(t).name;
        im=imread(imagename);

        for n=1:ncells
            [xv,yv] = ind2sub(sz,pixels{n,t})
            improxy(xv, yv)=1;
        end
        
        figure
        imshow(improxy, []), hold on
       for k=1:ncells
            if isempty(B{k,t})==0
               plot(B{k,t}(:,1),B{k,t}(:,2),'-g')
            else
                continue
            end
       end
          pause
          close all
    end  
        
elseif recrunch==1

    cd(savedir)
    load([basename '_cell05_dm.mat'])
    
    delind=[1,3];
    
    if ~isempty(delind)
        idx=setdiff(1:ncells, delind);
        icell_intensity=icell_intensity(idx,:);
        lcell=lcell(idx,:);
        B=B(idx,:);
        dl=dl(idx,:);
        ncells=length(idx);
        pixels=pixels(idx,:);
    end
    
else 
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    
    if LinTrack==1
        load([basename '_cell05_LT'], 'nnB', 'T', 'ncells', 'time', 'nnpixels', 'nnlcell')
        B = nnB;
        lcell = nnlcell;
        pixels = nnpixels;
    else
        load([basename '_cell05_BT'], 'B', 'T', 'ncells', 'time', 'pixels', 'lcell')
    end
        
    %pre-allocate variables
    if splitCell==1
        icell_intensity=nan(ncells*2, T);
    else
        icell_intensity=nan(ncells, T);
    end
    bg_intensity=nan(1,T);
    time=time./60;
    
    %find the final pre-lysis frame 
    dl=diff(lcell, 1, 2);
    lvg=mean(dl, 1, 'omitnan');
    [~, lidx]=min(lvg);

    for i=1:length(channels)
        
        cd(channels{i});
        
        imagename=fluo_directory{i}(lidx).name;
        [p1, p2]=getBackground(imagename);
        
     
        for t=1:T
            t
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            bg_intensity(1, t)=measureBackground(im, p1, p2);
            
            if splitCell==1
                p = 0;
                for n=1:ncells

                    pixLength = length(pixels{n,t});
                    pixHalf = round(pixLength/2);
                    c1 = pixels{n,t}(1:pixHalf);
                    c2 = pixels{n,t}(pixHalf+1:end);
                    p = p + 1;
                    icell_intensity(p,t)=mean(im(c1));      
                    p = p + 1;
                    icell_intensity(p,t)=mean(im(c2));
                end
            else
               for n=1:ncells

                %intensity_value = im(pixels{n,t});
                %pixel_location = pixels{n,t};
                %figure, hold on, plot(pixel_location, intensity_value, 'or', 'LineStyle', 'none'), yline(bg_intensity(1, t), '--k'), ylim([0, Inf]), pause, close
                icell_intensity(n,t)=mean(im(pixels{n,t}));
               end
            end
            
        end  
       
    end
end

%% Troubleshooting
if troubleshoot==1
    cd(channels{1}); 
     for t=1:T
            t
            imagename=fluo_directory{1}(t).name;
            im=imread(imagename);


           figure
           imshow(im, [])
           hold on
           for k=1:ncells
                if isempty(B{k,t})==0
                   plot(B{k,t}(:,1),B{k,t}(:,2),'-g')
                else
                    continue
                end
           end
          pause
          close all
     end
     
elseif troubleshoot==2
    cd(channels{1});
     for k=1:ncells
            k
            imagename=fluo_directory{1}(T).name;
            im=imread(imagename);


           figure
           imshow(im, [])
           hold on
           for t=1:T
                if isempty(B{k,t})==0
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-g')
                else
                    continue
                end
           end
          pause
          close all
     end
end

%% Plot data
cd(savedir)
    
%plot cellular fluorescence traces
figure(1), hold on
plot(time, icell_intensity, '-g')
title('Intensity vs Time')
xlabel('Time (min)')
ylabel('Cellular Intensity (A.U.)')
ylim([0 Inf])
saveas(gcf, [basename,'_cell05_fullIntensity_dm.fig'])
saveas(gcf, [basename,'_cell05_fullIntensity_dm.png'])

%plot background fluorescence traces
figure(2), hold on
plot(time, bg_intensity, '-b')
title('Intensity vs Time')
xlabel('Time (min)')
ylabel('Background Intensity (A.U.)')
ylim([0 Inf])
saveas(gcf, [basename,'_cell05_bgIntensity_dm.fig'])
saveas(gcf, [basename,'_cell05_bgIntensity_dm.png'])
   
save([basename '_cell05_dm.mat'])

cd(savedir2)
save([basename '_cell05_dm.mat'])

%% Functions
 function [p1, p2]=getBackground(imagename)

         %Load last image
         %imagename=fluo_directory{i}(t).name;
         im2=imread(imagename);

         %Determine Background
         figure,imshow(im2,[]), hold on, title('Select Background')
         k=waitforbuttonpress;
         set(gcf,'Pointer')
         hold on
         axis manual
         point1=get(gca,'CurrentPoint');
         finalRect=rbbox;
         point2=get(gca,'CurrentPoint');
         point1=point1(1,1:2);
         point2=point2(1,1:2);
         point1(point1<1)=1;
         point2(point2<1)=1;
         p1=min(point1,point2);%Calculate locations
         p2=max(point1,point2);
         offset = abs(point1-point2);%And dimensions
         x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
         y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
         plot(x,y);
         p1=round(p1);
         p2=round(p2);  
 end 

 function bglevel = measureBackground(im2, p1, p2)

         %Load last image
         %imagename=fluo_directory{i}(t).name;
         %im2=imread(imagename);

         %Determine background
         backim=im2(p1(2):p2(2),p1(1):p2(1));
         [counts,bins]=imhist(backim);
         [~,binnum]=max(counts);
         maxpos=bins(binnum);
         bglevel=mean(mean(backim));

 end 