%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all

%% analysis of FITC-K diffusion into the cell
maindir=['/Users/zarina/Documents/MATLAB/MatlabReady/07232021_analysis'];
cd(maindir);
dataFiles=dir('*BTfluo.mat');

preRatio=[];
postRatio=[];
for f=1:height(dataFiles)
    load(dataFiles(f).name)
    preRatio=[preRatio; ratio{2}(:, 5:9)];
    postRatio=[postRatio; ratio{2}(:, 40:44)];
end

maxPre=max(preRatio,[],2);
maxPost=max(postRatio,[],2);


meanPre=mean(preRatio,2);
meanPost=mean(postRatio,2);

x=[zeros(1, height(maxPre)), ones(1, height(maxPost))];
y=[maxPre',maxPost'];

figure
scatter(x,y)
ylabel('Ratio');
xlim([-0.5, 1.5])
title('Intensity Ratio Pre- and Post-Proteinase K Treatment')
xticklabels({'', 'Pre-Proteinase Treatment', '', 'Post-Proteinase Treatment', ''});

cd(maindir);
saveas(gcf, 'maxRatios.png')
saveas(gcf, 'maxRatios.fig')
close

h1=histogram(maxPre'), hold on
h2=histogram(maxPost')
title('Intensity Ratio Pre- and Post-Proteinase K Treatment')
xlabel('Max Ratio')
cd(maindir);
saveas(gcf, 'maxHist.png')
saveas(gcf, 'maxHist.fig')

pause, close

h1=histogram(meanPre'), hold on
h2=histogram(meanPost')
pause, close

%% Separate analysis
%%User Input
basename='07222021_Exp2';%Name of the image stack, used to save file.
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony1/' basename '_phase/' basename '_aligned'];%Directory that the image stack is saved in.
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
cd(savedir)
load([basename '_BTphase'], 'T', 'pixels', 'time', 'ncells')

tidx=17;

cd(dirname);
directory=dir('*.tif');

cell_temp=nan(ncells, T-tidx+1);
icell=[];
 for t=tidx:T       
    j=t-tidx+1;        
    imagename=directory(t).name;
    im=imread(imagename);

    for n=1:ncells
        cell_temp(n,j)=mean(im(pixels{n,t}));
    end

 end
 
 icell=[icell; cell_temp];
 
dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony2/' basename '_phase/' basename '_aligned'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony2/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
cd(savedir)
load([basename '_BTphase'], 'T', 'pixels', 'time', 'ncells')
 
cd(dirname);
directory=dir('*.tif');

cell_temp=nan(ncells, T-tidx+1);
 for t=tidx:T       
    j=t-tidx+1;         
    imagename=directory(t).name;
    im=imread(imagename);

    for n=1:ncells
        cell_temp(n,j)=mean(im(pixels{n,t}));
    end

 end
 
icell=[icell; cell_temp];

dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony3/' basename '_phase/' basename '_aligned'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony3/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
cd(savedir)
load([basename '_BTphase'], 'T', 'pixels', 'time', 'ncells')
 
cd(dirname);
directory=dir('*.tif');

cell_temp=nan(ncells, T-tidx+1);
 for t=tidx:T       
    j=t-tidx+1;         
    imagename=directory(t).name;
    im=imread(imagename);

    for n=1:ncells
        cell_temp(n,j)=mean(im(pixels{n,t}));
    end

 end
 
icell=[icell; cell_temp];

dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony4/' basename '_phase/' basename '_aligned'];
savedir=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony4/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
cd(savedir)
load([basename '_BTphase'], 'T', 'pixels', 'time', 'ncells')
 
cd(dirname);
directory=dir('*.tif');

cell_temp=nan(ncells, T-tidx+1);
 for t=tidx:T       
    j=t-tidx+1;         
    imagename=directory(t).name;
    im=imread(imagename);

    for n=1:ncells
        cell_temp(n,j)=mean(im(pixels{n,t}));
    end

 end
 
icell=[icell; cell_temp];

figure, hold on
for n=1:height(icell)
    plot(time(tidx:end)-time(tidx), icell(n,:))
end
%xline(time(tidx), '--', {'post-lysis'})
ylim([0, Inf])
xlabel('Time (s)')
ylabel('Intensity (A.U.)')
title('Loss of Phase Contrast Post-Lysis')
saveas(gcf, '07222021_Exp2_phaseLoss.fig')
saveas(gcf, '07222021_Exp2_phaseLoss.png')
