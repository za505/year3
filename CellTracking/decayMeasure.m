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

%% old code

Z=cell(ncells,length(tidx));
cg=0; %green count
cc=0; %cherry count

for n=1:ncells
    
    if ismember(n, nGreen)
        cg=cg+1;
        
        [~, locB]=ismember(n, nGreen);
        cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_figures'])
        
        dx=max(B{n,1}(:,1))-min(B{n,1}(:,1)); %columns are the x direction
        dy=max(B{n,1}(:,2))-min(B{n,1}(:,2)); %rows are the y direction
        x=min(B{n,1}(:,1))-dx/2:max(B{n,1}(:,1))+dx/2;
        y=min(B{n,1}(:,2))-dy/2:max(B{n,1}(:,2))+dy/2;
        [X,Y] = meshgrid(x,y);
        
        for j=1:length(tidx)
            t=tidx(j);
            Z{n,j}=ones(size(B{n,t},1),2);
            Z{n,j}=Z{1,t}*65530;
        end

        v = VideoWriter(strcat('mNeonGreen_', num2str(cg), '_intensity'),'MPEG-4');
        open(v);

        figure, hold on
        for j=1:length(tidx)
            t=tidx(j);
            plot3(B{n,t}(:,1),B{n,t}(:,2),Z{n,j})
            hold on
            t=surf(X,Y,cyto{1,t})
            %rotate3d on;
            t.EdgeColor = 'interp';
            t.FaceColor = 'interp';
            view(0,95)
            frame = getframe(gcf);
            writeVideo(v,frame);
            pause(1)
            clf
        end

        close(v)
        close all

    elseif ismember(n, nCherry)
        cc=cc+1;
        
        [~, locB]=ismember(n, nCherry);
        cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_figures'])
        
        dx=max(B{n,1}(:,1))-min(B{n,1}(:,1)); %columns are the x direction
        dy=max(B{n,1}(:,2))-min(B{n,1}(:,2)); %rows are the y direction
        x=min(B{n,1}(:,1))-dx/2:max(B{n,1}(:,1))+dx/2;
        y=min(B{n,1}(:,2))-dy/2:max(B{n,1}(:,2))+dy/2;
        [X,Y] = meshgrid(x,y);
        
        for j=1:length(tidx)
            t=tidx(j);
            Z{n,j}=ones(size(B{n,t},1),2);
            Z{n,j}=Z{1,t}*65530;
        end

        v = VideoWriter(strcat('mCherry_', num2str(cc), '_intensity'),'MPEG-4');
        open(v);

        figure, hold on
        for j=1:length(tidx)
            t=tidx(j);
            plot3(B{n,t}(:,1),B{n,t}(:,2),Z{n,j})
            hold on
            t=surf(X,Y,cyto{1,t})
            %rotate3d on;
            t.EdgeColor = 'interp';
            t.FaceColor = 'interp';
            view(0,95)
            frame = getframe(gcf);
            writeVideo(v,frame);
            pause(1)
            clf
        end

        close(v)
        close all
    else
        continue
    end
end
        
cd(savedir)
saveas([basename 'processed_dm.mat'])

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