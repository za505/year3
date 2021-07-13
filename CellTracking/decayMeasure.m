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
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_reanalysis'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mNeonGreen/'  basename '_aligned']; ['/Users/zarina/Downloads/NYU/Year3_2021_Summer/06062021_analysis/' basename '_colony1/' basename '_mCherry/'  basename '_aligned']}; 
recrunch=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savedir)
    load([basename '_dm.mat'])
else
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(savedir)
    load([basename '_BTphase'], 'B', 'T', 'ncells', 'time', 'pixels', 'lcell')

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
            if i==1 & mean(intensity(n, 1:5)>500)
                icell_green=[icell_green; intensity(n,:)];
                nGreen=[nGreen n];
            elseif i==2 & mean(intensity(n, 1:5)>1000)
                icell_mcherry=[icell_mcherry; intensity(n,:)];
                nCherry=[nCherry n];
            else
                continue
            end
        end
    end
    
    cd(savedir)
    save([basename '_dm.mat'])
end

%% Plot data
%plot to see single traces of mNeonGreen cells
figure(1), hold on
for i=1:height(icell_green)
    plot(time(tidx), icell_green(i,:))
    x=time(tidx(end));
    y=icell_green(i,end);
    text(x,y, num2str(i));
end
xline(tpt, '--', {'Membrane Lysis'})
title('Cellular Intensity of mNeonGreen vs Time')
xlabel('Time (s)')
ylabel('Cellular Intensity (A.U.)')
saveas(gcf, [basename,'_fullGreen.fig'])
saveas(gcf, [basename,'_fullGreen.png'])

figure(2), hold on
for i=1:length(nGreen)
    n=nGreen(i);
    plot(time(tidx), lcell(n,tidx))
    x=time(tidx(end));
    y=lcell(n,end);
    text(x,y, num2str(i));
end
xline(tpt, '--', {'Membrane Lysis'})
title('mNeonGreen Cell Length vs Time')
xlabel('Time (s)')
ylabel('Length (\mum)')
saveas(gcf, [basename,'_LTGreen.fig'])
saveas(gcf, [basename,'_LTGreen.png'])
 
%plot to see single traces of mCherry cells
figure(3), hold on
for i=1:height(icell_mcherry)
   plot(time(tidx), icell_mcherry(i,:))
   x=time(tidx(end));
    y=icell_mcherry(i,end);
    text(x,y, num2str(i));
end
xline(tpt, '--', {'Membrane Lysis'})
title('Cellular Intensity of mCherry vs Time')
xlabel('Time (s)')
ylabel('Cellular Intensity (A.U.)')
saveas(gcf, [basename,'_fullCherry.fig'])
saveas(gcf, [basename,'_fullCherry.png'])

figure(4), hold on
for i=1:length(nCherry)
    n=nCherry(i);
    plot(time(tidx), lcell(n,tidx))
    x=time(tidx(end));
    y=lcell(n,end);
    text(x,y, num2str(i));
end
xline(tpt, '--', {'Membrane Lysis'})
title('mCherry Cell Length vs Time')
xlabel('Time (s)')
ylabel('Length (\mum)')
saveas(gcf, [basename,'_LTCherry.fig'])
saveas(gcf, [basename,'_LTCherry.png'])

%%%%%%%%%%%Functions
function [y] = exponential(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=b(1)*exp(-x.*b(2));
end

function [y] = linear(b,x)
%this function calculates y=A*(e^-t/tau)
%where b(1)=A, b(2)=tau, x=t, and the cellular intensity=y;
y=x.*b(1)+b(2);
end