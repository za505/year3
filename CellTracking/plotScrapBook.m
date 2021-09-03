%%PlotScrapBook.m
%Purpose: to save code for one-off plots I've made

clear, close all

%% analysis of FITC-K diffusion into the cell
% maindir=['/Users/zarina/Documents/MATLAB/MatlabReady/07232021_analysis'];
% cd(maindir);
% dataFiles=dir('*BTfluo.mat');
% 
% preRatio=[];
% postRatio=[];
% for f=1:height(dataFiles)
%     load(dataFiles(f).name)
%     preRatio=[preRatio; ratio{2}(:, 5:9)];
%     postRatio=[postRatio; ratio{2}(:, 40:44)];
% end
% 
% maxPre=max(preRatio,[],2);
% maxPost=max(postRatio,[],2);
% 
% 
% meanPre=mean(preRatio,2);
% meanPost=mean(postRatio,2);
% 
% x=[zeros(1, height(maxPre)), ones(1, height(maxPost))];
% y=[maxPre',maxPost'];
% 
% figure
% scatter(x,y)
% ylabel('Ratio');
% xlim([-0.5, 1.5])
% title('Intensity Ratio Pre- and Post-Proteinase K Treatment')
% xticklabels({'', 'Pre-Proteinase Treatment', '', 'Post-Proteinase Treatment', ''});
% 
% cd(maindir);
% saveas(gcf,gcf, 'maxRatios.png')
% saveas(gcf,gcf, 'maxRatios.fig')
% close
% 
% h1=histogram(maxPre'), hold on
% h2=histogram(maxPost')
% title('Intensity Ratio Pre- and Post-Proteinase K Treatment')
% xlabel('Max Ratio')
% cd(maindir);
% saveas(gcf,gcf, 'maxHist.png')
% saveas(gcf,gcf, 'maxHist.fig')
% 
% pause, close
% 
% h1=histogram(meanPre'), hold on
% h2=histogram(meanPost')
% pause, close
% 
% %% Separate analysis
% %%User Input
% basename='07222021_Exp2';%Name of the image stack, used to save file.
% dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony1/' basename '_phase/' basename '_aligned'];%Directory that the image stack is saved in.
% dirsave=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony1/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
% cd(dirsave)
% load([basename '_BTphase'], 'T', 'pixels', 'time', 'ncells')
% 
% tidx=17;
% 
% cd(dirname);
% directory=dir('*.tif');
% 
% cell_temp=nan(ncells, T-tidx+1);
% icell=[];
%  for t=tidx:T       
%     j=t-tidx+1;        
%     imagename=directory(t).name;
%     im=imread(imagename);
% 
%     for n=1:ncells
%         cell_temp(n,j)=mean(im(pixels{n,t}));
%     end
% 
%  end
%  
%  icell=[icell; cell_temp];
%  
% dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony2/' basename '_phase/' basename '_aligned'];
% dirsave=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony2/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
% cd(dirsave)
% load([basename '_BTphase'], 'T', 'pixels', 'time', 'ncells')
%  
% cd(dirname);
% directory=dir('*.tif');
% 
% cell_temp=nan(ncells, T-tidx+1);
%  for t=tidx:T       
%     j=t-tidx+1;         
%     imagename=directory(t).name;
%     im=imread(imagename);
% 
%     for n=1:ncells
%         cell_temp(n,j)=mean(im(pixels{n,t}));
%     end
% 
%  end
%  
% icell=[icell; cell_temp];
% 
% dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony3/' basename '_phase/' basename '_aligned'];
% dirsave=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony3/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
% cd(dirsave)
% load([basename '_BTphase'], 'T', 'pixels', 'time', 'ncells')
%  
% cd(dirname);
% directory=dir('*.tif');
% 
% cell_temp=nan(ncells, T-tidx+1);
%  for t=tidx:T       
%     j=t-tidx+1;         
%     imagename=directory(t).name;
%     im=imread(imagename);
% 
%     for n=1:ncells
%         cell_temp(n,j)=mean(im(pixels{n,t}));
%     end
% 
%  end
%  
% icell=[icell; cell_temp];
% 
% dirname=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony4/' basename '_phase/' basename '_aligned'];
% dirsave=['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07222021_analysis/' basename '_colony4/' basename '_phase/' basename '_figures'];%Directory to save the output .mat file to.
% cd(dirsave)
% load([basename '_BTphase'], 'T', 'pixels', 'time', 'ncells')
%  
% cd(dirname);
% directory=dir('*.tif');
% 
% cell_temp=nan(ncells, T-tidx+1);
%  for t=tidx:T       
%     j=t-tidx+1;         
%     imagename=directory(t).name;
%     im=imread(imagename);
% 
%     for n=1:ncells
%         cell_temp(n,j)=mean(im(pixels{n,t}));
%     end
% 
%  end
%  
% icell=[icell; cell_temp];
% 
% figure, hold on
% for n=1:height(icell)
%     plot(time(tidx:end)-time(tidx), icell(n,:))
% end
% %xline(time(tidx), '--', {'post-lysis'})
% ylim([0, Inf])
% xlabel('Time (s)')
% ylabel('Intensity (A.U.)')
% title('Loss of Phase Contrast Post-Lysis')
% saveas(gcf,gcf, '07222021_Exp2_phaseLoss.fig')
% saveas(gcf,gcf, '07222021_Exp2_phaseLoss.png')
datadir="/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis";
dirsave="/Users/zarina/Downloads/NYU/Year3_2021_Fall/presentations/09062021";
% %%
% cd(datadir)
% load("07162021_Exp3_colony2_dm.mat")
% colony2_LBa=dataTable(dataTable.halfie==0, 3);
% colony2_LBa=table2array(colony2_LBa);
% 
% load("07162021_Exp3_colony3_dm.mat")
% colony3_LBa=dataTable(dataTable.halfie==0, 3);
% colony3_LBa=table2array(colony3_LBa);
% 
% figure(1), hold on
% plot(time, colony2_LBa, '-r')
% plot(time, colony3_LBa, '-b')
% xlabel('Time')
% ylabel('Normalized Fluorescence')
% title('LB, fr=2 min')
% cd(dirsave)
% saveas(gcf,'LBa_positions.png');
% saveas(gcf,'LBa_positions.fig');
% % close
% % 
% % norm_green=[colony2_LBa; colony3_LBa];
% % fitScore_LBa=expFitting(norm_green, time);
% %%
% cd(datadir)
% load("08262021_Exp1_colony1_dm.mat")
% colony1_LBb=dataTable(dataTable.halfie==0, 3);
% colony1_LBb=table2array(colony1_LBb);
% 
% load("08262021_Exp1_colony2_dm.mat")
% colony2_LBb=dataTable(dataTable.halfie==0, 3);
% colony2_LBb=table2array(colony2_LBb);
% 
% load("08262021_Exp1_colony3_dm.mat")
% colony3_LBb=dataTable(dataTable.halfie==0, 3);
% colony3_LBb=table2array(colony3_LBb);
% 
% load("08262021_Exp1_colony4_dm.mat")
% colony4_LBb=dataTable(dataTable.halfie==0, 3);
% colony4_LBb=table2array(colony4_LBb);
% 
% load("08262021_Exp1_colony5_dm.mat")
% colony5_LBb=dataTable(dataTable.halfie==0, 3);
% colony5_LBb=table2array(colony5_LBb);
% 
% load("08262021_Exp1_colony6_dm.mat")
% colony6_LBb=dataTable(dataTable.halfie==0, 3);
% colony6_LBb=table2array(colony6_LBb);
% 
% figure(2), hold on
% plot(time, colony1_LBb, '-m')
% plot(time, colony2_LBb, '-c')
% plot(time, colony3_LBb, '-r')
% plot(time, colony4_LBb, '-g')
% plot(time, colony5_LBb, '-b')
% plot(time, colony6_LBb, '-k')
% xlabel('Time')
% ylabel('Normalized Fluorescence')
% title('LB, fr=1 min')
% cd(dirsave)
% saveas(gcf,'LBb_positions.png');
% saveas(gcf,'LBb_positions.fig');
% %%
% cd(datadir)
% load("07162021_Exp4_colony1_dm.mat")
% colony1_EDTA=dataTable(dataTable.halfie==0, 3);
% colony1_EDTA=table2array(colony1_EDTA);
% 
% load("07162021_Exp4_colony2_dm.mat")
% colony2_EDTA=dataTable(dataTable.halfie==0, 3);
% colony2_EDTA=table2array(colony2_EDTA);
% 
% load("07162021_Exp4_colony3_dm.mat")
% colony3_EDTA=dataTable(dataTable.halfie==0, 3);
% colony3_EDTA=table2array(colony3_EDTA);
% 
% load("07162021_Exp4_colony4_dm.mat")
% colony4_EDTA=dataTable(dataTable.halfie==0, 3);
% colony4_EDTA=table2array(colony4_EDTA);
% 
% figure(3), hold on
% plot(time, colony1_EDTA, '-r')
% plot(time, colony2_EDTA, '-g')
% plot(time, colony3_EDTA, '-b')
% plot(time, colony4_EDTA, '-c')
% xlabel('Time')
% ylabel('Normalized Fluorescence')
% title('EDTA, fr=2 min')
% cd(dirsave)
% saveas(gcf,'EDTA_positions.png');
% saveas(gcf,'EDTA_positions.fig');
% 
% cd(datadir)
% load("07152021_Exp1_colony1_dm.mat")
% colony1_Mg=dataTable(dataTable.halfie==0, 3);
% colony1_Mg=table2array(colony1_Mg);
% 
% load("07152021_Exp1_colony2_dm.mat")
% colony2_Mg=dataTable(dataTable.halfie==0, 3);
% colony2_Mg=table2array(colony2_Mg);
% 
% load("07152021_Exp1_colony3_dm.mat")
% colony3_Mg=dataTable(dataTable.halfie==0, 3);
% colony3_Mg=table2array(colony3_Mg);
% 
% figure(4), hold on
% plot(time, colony1_Mg, '-r')
% plot(time, colony2_Mg, '-c')
% plot(time, colony3_Mg, '-b')
% xlabel('Time')
% ylabel('Normalized Fluorescence')
% title('Mg, fr=2 min')
% cd(dirsave)
% saveas(gcf,'Mg_positions.png');
% saveas(gcf,'Mg_positions.fig');

%%
cd(datadir)
load("08272021_Exp1_colony1_dm.mat")
colony1_PBS=dataTable(dataTable.halfie==0, 3);
colony1_PBS=table2array(colony1_PBS);

load("08272021_Exp1_colony2_dm.mat")
colony2_PBS=dataTable(dataTable.halfie==0, 3);
colony2_PBS=table2array(colony2_PBS);

load("08272021_Exp1_colony3_dm.mat")
colony3_PBS=dataTable(dataTable.halfie==0, 3);
colony3_PBS=table2array(colony3_PBS);

figure(5), hold on
plot(time, colony1_PBS, '-r')
plot(time, colony2_PBS, '-c')
plot(time, colony3_PBS, '-b')
xlabel('Time')
ylabel('Normalized Fluorescence')
title('PBS, fr=1 min')
cd(dirsave)
saveas(gcf,'PBS_positions.png');
saveas(gcf,'PBS_positions.fig');

%%
%comparison
% figure(6), hold on
% plot(time, colony1_LBb(:,2:end), '-r')
% plot(time, colony2_LBb(:, 2:end), '-r')
% plot(time, colony3_LBb(:, 2:end), '-r')
% plot(time, colony4_LBb(:, 2:end), '-r')
% plot(time, colony5_LBb(:, 2:end), '-r')
% plot(time, colony6_LBb(:, 2:end), '-r')
% plot(time, colony1_PBS, '-b')
% plot(time, colony2_PBS, '-b')
% plot(time, colony3_PBS, '-b')
% xlabel('Time')
% ylabel('Normalized Fluorescence')
% title('Diffusion of LB (red) and PBS (blue) treated cells, fr=1 min')
% cd(dirsave)
% saveas(gcf,'LB_PBS_compare.png');
% saveas(gcf,'LB_PBS_compare.fig');

%% Fit curves
close all

norm_green = [colony1_PBS; colony2_PBS; colony3_PBS];
fitScore=expFitting(norm_green, time)

% [iM, iN]=size(im);
% canvas=zeros(iM, iN);
% 
% for n=1:ncells
%     
%     tiledlayout(1,2);
%     frame=canvas;
%     frame(pixels{n,T})=1;
%     %frame=uint16(frame);
%     
%     nexttile
%     imshow(frame)
%     
%     nexttile
%     plot(time, norm_green(n,:))
%     
%     pause, close
% end

%% functions
function [a, b, c, d] = expFitting(norm_green, time)
    fitScore=nan(height(norm_green), 1);
    for n=1:height(norm_green)
        if isnan(norm_green(n, 1))==0
            [xData, yData]=prepareCurveData(time, norm_green(n,:));
            figure
            f=fit(xData, yData, 'exp2'); %fit to a linear exponential
            plot(f, xData, yData)
            prompt = 'Q1: Is this a good fit? 0=No, 1=Yes ';
            fitScore(n) = input(prompt);
            close
        else
            fitScore(n)=NaN;
        end
    end
end