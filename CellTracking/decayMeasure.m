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
basename='07102021_Exp1';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/' basename '_colony2/' basename '_phase2/' basename '_figures'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/' basename '_colony2/' basename '_GFP/' basename '_figures'];%Directory to save the output .mat file to.
channels={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/' basename '_colony2/' basename '_GFP/'  basename '_aligned']}; 
recrunch=0;
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
    cd(dirname)
    load([basename '_BTphase'], 'B', 'T', 'ncells', 'time', 'pixels', 'lcell')

    %pre-allocate variables
    icell_green=nan(ncells, T);
    f={};
    
    for i=1:length(channels)
        
        cd(channels{i});
        %intensity=nan(ncells, T);
        
        for t=1:T
             
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            for n=1:ncells
                icell_green(n,t)=mean(im(pixels{n,t}));
            end
            
        end
        
        for n=1:ncells
            [xData, yData] = prepareCurveData(time, icell_green(n,:));
            f{n}=fit(xData, yData, 'exp1');
        end
    end
    
    cd(savedir)
    save([basename '_dm.mat'])
end

%% Plot data
cd(savedir)
%plot to see single traces of mNeonGreen cells
figure(1), hold on
for i=1:height(icell_green)
    plot(time, icell_green(i,:))
    x=time(end);
    y=icell_green(i,end);
    text(x,y, num2str(i));
end
%xline(tpt, '--', {'Membrane Lysis'})
title('Cellular Intensity of mNeonGreen vs Time')
xlabel('Time (s)')
ylabel('Cellular Intensity (A.U.)')
saveas(gcf, [basename,'_fullGreen.fig'])
saveas(gcf, [basename,'_fullGreen.png'])

figure(2), hold on
for n=1:length(icell_green)
    plot(time, lcell(n, :))
    x=time(end);
    y=lcell(n,end);
    text(x,y, num2str(n));
end
%xline(tpt, '--', {'Membrane Lysis'})
title('mNeonGreen Cell Length vs Time')
xlabel('Time (s)')
ylabel('Length (\mum)')
saveas(gcf, [basename,'_LTGreen.fig'])
saveas(gcf, [basename,'_LTGreen.png'])

%plot to see single traces of mNeonGreen cells
for i=1:height(icell_green)
    figure('Name', num2str(i))
    plot(f{i}, time, icell_green(i,:))
    title(['Cellular Intensity of mNeonGreen vs Time, ' '#' num2str(i)])
    xlabel('Time (s)')
    ylabel('Cellular Intensity (A.U.)')
    saveas(gcf, [basename '_' num2str(i) '_fitGreen.fig'])
    saveas(gcf, [basename '_' num2str(i) '_fitGreen.png'])
    close
end

%% Calculate time constant histograms
dir1=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/' basename '_colony1/' basename '_GFP/' basename '_figures'];
dir2=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/' basename '_colony2/' basename '_GFP/' basename '_figures'];

cd(dir1)
load([basename '_dm.mat'], 'f')
f1=f;

cd(dir2)
load([basename '_dm.mat'], 'f')
f2=f;

tconst=[];
for i=1:2
    if i==1
        for n=1:length(f1)
            tconst=[tconst; -1/f1{n}.b];
        end
    elseif n==2
        for n=1:length(f2)
            tconst=[tconst; -1/f2{n}.b];
        end
    end
end

cd(['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/'])
bins=length(f1)+length(f2);
 histogram(tconst,bins)
 title('Histogram of Time Constants, n=37')
 xlabel('Time Scale Tau (s-)')
 saveas(gcf, [basename '_histogram.fig'])
 saveas(gcf, [basename '_histogram.png'])