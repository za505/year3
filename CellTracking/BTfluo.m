%BTfluo.m
%Rico Rojas, updated 1/21/19
%Calculates the average cytoplasmic fluorescence intensity from cell
%tracked with BacTrack.m.  

clear, close all

%INSTRUCTIONS FOR USE:
%Remove frames with poor contrast and save fluorescent image stacks
%directories by themselves. 

%INPUT
%basename: name to save the results to.
%fluordir: name of folder containing fluorescent image stacks to quantify.

%OUTPUT:
%icell: a matrix (ncellxT) with the fluorescent intensities of each
        %cell, where rows are the cells and columns are time points.
%icell_av: Each entry is a vector containing the population-
        %average of the single-cell fluorescent intensities.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='10102021_Exp6';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/10102021_analysis/' basename '/' basename '_colony4/' basename '_phase/' basename '_figures'];%Directory that the BTphase.mat file is saved in
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/10102021_analysis/'  basename '/' basename '_colony4/' basename '_FD10/' basename '_figures'];%Directory to save the output .mat file to.
fluordir=['/Users/zarina/Downloads/NYU/Year3_2021_Fall/10102021_analysis/'  basename '/' basename '_colony4/' basename '_FD10/' basename '_aligned']; 
recrunch=0;
troubleshoot=4;
fluorFrames=[12:21, 26:49]; %frames where FD10 is perfused
dotThresh=500; %this is the lightest pixel value of the dot during FD10 perfusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==0;
    
    cd(fluordir); 
    fluo_directory=dir('*.tif');

    cd(dirname)
    load([basename '_BTphase'], 'pixels', 'ncells', 'lcell', 'B', 'time', 'imM', 'imN', 'directory')

    T=length(fluo_directory);
    icell=nan(ncells, T);
    icell_av=nan(ncells, T);
    ratio=nan(ncells, T);
    ratio_av=nan(ncells, T);

    bgIntensity=nan(ncells, T);
    bgPixels=cell(ncells, T);
    Cin=nan(ncells, T);
    Cout=nan(ncells, T);
    
    %nhood=[0 0 1 0 0; 0 0 1 0 0; 1 1 1 1 1; 0 0 1 0 0; 0 0 1 0 0];
    nhood=[repmat([0 0 0 0 0 1 0 0 0 0 0],5,1); ones(1,11);repmat([0 0 0 0 0 1 0 0 0 0 0],5,1)];
    
    allPixels=[];
    pixelTemp=[];
    
    %define a matrix that has all the pixel coordinates combined
    for n=1:ncells
        for t=round([1, T/2, T])
            pixelTemp=[pixelTemp; pixels{n,t}];
        end
        allPixels=[allPixels; unique(pixelTemp)];
    end
    allPixels=unique(allPixels);
   
    cd(fluordir);
    for t=1:T
        t
        imagename=fluo_directory(t).name;
        im=imread(imagename);

        for n=1:ncells
            
            %find local bg pixel locations
            im2=zeros(imM, imN);
            im2(pixels{n,t})=1;
            im2=imdilate(im2, nhood);
            im2(allPixels)=0;
            bgPixels{n,t}=find(im2==1);
           
%             if ismember(t,fluorFrames)==1
%                 idx=find(im(bgPixels{n,t})>dotThresh);
%                 bgPixels{n,t}=bgPixels{n,t}(idx);
%             end
            
            %calculate the mean local background fluorescence
            bgIntensity(n,t)=mean(im(bgPixels{n,t}), 'omitnan');
          
            %calculate the mean cellular fluorescence
            icell(n,t)=mean(im(pixels{n,t}));    
            
        end
    end

%rename the variables so I don't have to alter the rest of the code
%(previously this is where I would subtract autofluorescence)
Cout=bgIntensity; 
Cin=icell; 

%calculate average intensity and standard deviation
avgIntensity = mean(Cin, 1, 'omitnan');
stdIntensity = std(Cin, 0, 1, 'omitnan');

%now, calculate the Cin/Cout ratio
ratio = Cin ./ Cout;
%ratio(:,[1:frameInitial-1, frameAuto])=NaN; %note: 0/0 makes weird things happen

avgRatio = mean(ratio, 1, 'omitnan');
stdRatio = std(ratio, 0, 1, 'omitnan');

%find Cin-Cout
diffC=nan(ncells,T);
for n=1:ncells
    diffC(n,:)=Cin(n,:)-Cout(n,:);
end

dCin=nan(ncells,T);
deltat=[0 time(2:end)-time(1:end-1)];
for n=1:ncells
    dCin(n,:)=[0 diff(Cin(n,:))];
    dCin(n,:)=dCin(n,:)./deltat;
end

cd(savedir)
save([basename '_BTfluo'])

elseif recrunch==1
    cd(savedir)
    load ([basename '_BTfluo'])
    troubleshoot=0;
end

%% Troubleshooting
cd(fluordir); 
 if troubleshoot==1
     for t=1:T
            t
            imagename=fluo_directory(t).name;
            im=imread(imagename);


           figure
           imshow(im, [])
           hold on
           for k=1:ncells
                if isempty(B{k,t})==0
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
                else
                    continue
                end
           end
          pause
          close all
     end
 elseif troubleshoot==2
      for k=1:ncells
            k
            imagename=fluo_directory(T).name;
            im=imread(imagename);


           figure
           imshow(im, [])
           hold on
           for t=1:T
                if isempty(B{k,t})==0
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
                else
                    continue
                end
           end
          pause
          close all
      end
 elseif troubleshoot==3
   for n=1:ncells
      
     im3=zeros(imM, imN);
     
     figure, hold on
     for t=1:T
        im3(bgPixels{n,t})=1;
     end
     imshow(im3), pause, close
   end
   
 elseif troubleshoot==4
   cd(directory(T).folder)
   
   for n=1:ncells
      
    imagename=directory(T).name;
    im3=imread(imagename);
     
     figure, hold on
     for t=1:T
        im3(bgPixels{n,t})=max(max(im3));
     end
     imshow(im3, []), pause, close
   end
   
 end
%% Plot FD10 data
cd(savedir)

%Plot data
%Let's plot cellular intensity first
figure(1), hold on, title('Cin vs Time')
for n=1:ncells
    plot(time/60,Cin(n,:))
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(5, '--', {'FITC-dextran'})
xline(14, '--', {'*FITC-dextran'})
xline(20, '--', {'FITC-dextran + 20 mM Mg^{2+}'})
ylim([-3 Inf])
saveas(gcf, [basename,'_intensity.fig'])
saveas(gcf, [basename,'_intensity.png'])

%Plot adj background fluorescence
figure(2), hold on 
title('Cout vs Time')
for n=1:ncells
    plot(time/60,Cout(n,:))
end
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
fig2pretty
xline(5, '--', {'FITC-dextran'})
xline(14, '--', {'*FITC-dextran'})
xline(20, '--', {'FITC-dextran + 20 mM Mg^{2+}'})
ylim([-3 Inf])
saveas(gcf, [basename,'_Cout.fig'])
saveas(gcf, [basename,'_Cout.png'])

%Let's plot Cin-Cout
figure, hold on, title('Cin-Cout vs Time')
for n=1:ncells
    plot(time/60,Cin(n,:)-Cout(n,:))
end
%ylim([-3 Inf])
xline(5, '--', {'FITC-dextran'})
xline(14, '--', {'*FITC-dextran'})
xline(20, '--', {'FITC-dextran + 20 mM Mg^{2+}'})
fig2pretty
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
saveas(gcf, [basename,'_diffC.fig'])
saveas(gcf, [basename,'_diffC.png'])


