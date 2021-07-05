%decayMeasure.m
%Zarina Akbary, updated 07/04/21
%Calculates changes in fluor. intensity. Incorporates BTfluo.m code.

clear, close all

%INSTRUCTIONS FOR USE:
%run BacTrack.m first

%INPUT
%basename: experiments of interest
%dirname: where .mat files are stored
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%timescale=array of time scales
%Ts=array of T (number of frames)
%switchF=array of switchFrames
%A=array of A values (A*e^(-alpha*t)+y0)
%f=cell of coeff for exponential eqxn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename='06062021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_aligned']; ['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_aligned']}; 
recrunch=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(dirname)
    load([basename '_dm.mat'])
else
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    load([basename '_BTphase'], 'B', 'T', 'ncells', 'time', 'pixels')

    %pre-allocate variables
    icell_green=[];
    icell_mcherry=[];
    nGreen=[];
    nCherry=[];
    
    %calculate time index
    tidx=1:2:T-1; %remember, fluor images were taken every 2 frames
    tpt=480; %time at which *PBS+5% NLS is perfused
    
    for i=1:length(channels)
        
        cd(channels{i});
        intensity=nan(ncells, length(tidx));
        
        for j=1:length(tidx)
            t=tidx(j);
             
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            for n=1:ncells
                intensity(n,j)=mean(im(pixels{n,t}));
            end
            
        end
        
        for n=1:ncells
            if i==1 & mean(intensity(n, 1:5)>20000)
                icell_green=[icell_green; intensity(n,:)];
                nGreen=[nGreen n];
            elseif i==2 & mean(intensity(n, 1:5)>20000)
                icell_mcherry=[icell_mcherry; intensity(n,:)];
                nCherry=[nCherry n];
            else
                continue
            end
        end
    end
    
    cd(dirname)
    save([basename '_dm.mat'])
end
%% Plot data
%plot to see single traces of mNeonGreen cells
figure(1), hold on
for i=1:height(icell_green)
    plot(time(tidx), icell_green(i,:))
end
xline(tpt, '--', {'Membrane Lysis'})
title('Cellular Intensity of mNeonGreen vs Time')
xlabel('Time (s)')
ylabel('Cellular Intensity (A.U.)')
saveas(gcf, [basename,'_fullGreen.fig'])
saveas(gcf, [basename,'_fullGreen.png'])
 
%plot to see single traces of mCherry cells
figure(2), hold on
for i=1:height(icell_mcherry)
   plot(time(tidx), icell_mcherry(i,:))
end
xline(tpt, '--', {'Membrane Lysis'})
title('Cellular Intensity of mCherry vs Time')
xlabel('Time (s)')
ylabel('Cellular Intensity (A.U.)')
saveas(gcf, [basename,'_fullCherry.fig'])
saveas(gcf, [basename,'_fullCherry.png'])

%% make movies
cg=0; %green count
cc=0; %cherry count

imgreen=cell(length(nGreen),length(tidx));
imcherry=cell(length(nCherry),length(tidx));

for n=1:ncells
    n
    
    if ismember(n, nGreen)
        cg=cg+1;
        
        [~, locB]=ismember(n, nGreen);
        cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
        
        m=min(find(cellfun('length', B(n,:))>0)); %find the initial boundaries
        dx=max(B{n,m}(:,1))-min(B{n,m}(:,1)); %columns are the x direction
        dy=max(B{n,m}(:,2))-min(B{n,m}(:,2)); %rows are the y direction
        r1=round(min(B{n,m}(:,1))-dx/2);
        r2=round(min(B{n,m}(:,2))-dy/2);
        d1=max([dx dy]);
        
        cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
        v = VideoWriter(strcat('mNeonGreen_', num2str(cg), '_intensity'),'MPEG-4');
        open(v);
        
        for j=1:length(tidx)
            
            t=tidx(j)
            cd(channels{1});
            imagename=fluo_directory{1}(t).name;
            im=imread(imagename);
            
            
            imgreen{n,j} = imcrop(im,[r1 r2 d1*2 d1*2]);
            imshow(imgreen{n,j})
            frame = getframe(gcf);
            writeVideo(v,frame);
            pause(1)
            clf

        end
        
        close(v)
        close all
        
%         [imM,imN]=size(imgreen{n,m});
%         y=1:imM; %rows are the y direction
%         x=1:imN; %columns are the x direction
%         [X,Y] = meshgrid(x,y);
%         
%         cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
%         v = VideoWriter(strcat('mNeonGreen_', num2str(cg), '_intensity'),'MPEG-4');
%         open(v);
%         
%         for j=1:length(tidx)
%             figure, hold on
%             t=surf(X,Y,imgreen{n,j})
%             %rotate3d on;
%             t.EdgeColor = 'interp';
%             t.FaceColor = 'interp';
%             view(0,95)
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%             pause(0.1)
%             clf
%         end
%         
%         close(v)
%         close all

    elseif ismember(n, nCherry)
        cc=cc+1;
        
        [~, locB]=ismember(n, nCherry);
        cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_figures'])
        
        m=min(find(cellfun('length', B(n,:))>0)); %find the initial boundaries
        dx=max(B{n,m}(:,1))-min(B{n,m}(:,1)); %columns are the x direction
        dy=max(B{n,m}(:,2))-min(B{n,m}(:,2)); %rows are the y direction
        r1=round(min(B{n,m}(:,1))-dx/2);
        r2=round(min(B{n,m}(:,2))-dy/2);
        d1=max([dx dy]);
        
        cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_figures'])
        v = VideoWriter(strcat('mCherry_', num2str(cg), '_intensity'),'MPEG-4');
        open(v);
        for j=1:length(tidx)
            
            t=tidx(j)
            cd(channels{2});
            imagename=fluo_directory{2}(t).name;
            im=imread(imagename);
            
            imcherry{n,j} = imcrop(im,[r1 r2 d1*2 d1*2]);
            imshow(imcherry{n,j})
            frame = getframe(gcf);
            writeVideo(v,frame);
            pause(1)
            clf

        end
        
        close(v)
        close all
        
%         [imM,imN]=size(imcherry{n,m});
%         y=1:imM; %rows are the y direction
%         x=1:imN; %columns are the x direction
%         [X,Y] = meshgrid(x,y);
%         
%         cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_figures'])
%         v = VideoWriter(strcat('mCherry_', num2str(cg), '_intensity'),'MPEG-4');
%         open(v);
%         
%         for j=1:length(tidx)
%             figure, hold on
%             t=surf(X,Y,imcherry{n,j})
%             %rotate3d on;
%             t.EdgeColor = 'interp';
%             t.FaceColor = 'interp';
%             view(0,95)
%             frame = getframe(gcf);
%             writeVideo(v,frame);
%             pause(0.1)
%             clf
%         end
%         
%         close(v)
%         close all
    else
        continue
    end
end
        
cd(savedir)
save([basename 'processed_dm.mat'])

%%%%%%%%%%%Functions
function [y] = exponential(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=b(1)*exp(-x./b(2));
end

function [y] = linear(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=x.*b(1)+b(2);
end