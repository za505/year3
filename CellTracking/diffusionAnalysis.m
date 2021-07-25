%diffusionAnalysis.m
%Author: Zarina Akbary
%Date: 17 July 2021
%purpose: to analyze mNeonGreen.m diffusion out of the cell

clear, close all
%% User Input
maindir=['/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis'];
photodir={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/07102021_Exp1_colony1/07102021_Exp1_GFP/07102021_Exp1_figures'];['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07102021_analysis/07102021_Exp1_colony2/07102021_Exp1_GFP/07102021_Exp1_figures']};
recrunch=1;
replot=1;
%%%%%%%%%%%%%%%%%%%

%% Part I: New Tables
if recrunch==0
    % Load data
    cd(maindir)
    dataDirectory=dir('*_dm.mat');

    fullData=table();

    for f=1:height(dataDirectory)

        load(dataDirectory(f).name)
        %tempData = dataTable(dataTable.halfie == 0, :);
        tempData = dataTable;
        
        %label experiment
        experimentName=repelem(string(dataDirectory(f).name(1:end-7)), height(tempData))';

        if f<=3
            condition=repelem(string('Mg2+'), height(tempData))';
        elseif f>3 & f<6
            condition=repelem(string('LB'), height(tempData))';
        elseif f>=6
            condition=repelem(string('EDTA'), height(tempData))';
        end
        
        %calculate cell length changes
        %from LB > LB + 5% detergent is t=6 and t=7
        dL=((lcell(:, 7)-lcell(:, 6))./lcell(:,6))*100; %lf-li/li
        
        %from LB + 5% detergent to the condition 
        tf=mean(lcell(:, 10:13), 2, 'omitnan');
        dL2=((tf-lcell(:,8))./lcell(:,8))*100; %lf-li/li
        %dL2=((lcell(:, 8)-lcell(:, 10))./lcell(:,8))*100; %lf-li/li
        
        %from LB to the condition 
        dL3=((tf-lcell(:,6))./lcell(:,6))*100; %lf-li/li
        
        %create a table from these new variables
        newVar=table(experimentName, condition, dL, dL2,dL3, 'VariableNames', {'Experiment', 'Condition', 'Lysis dL', 'Post-Lysis dL', 'Overall dL'});
        
        tempData=[newVar, tempData];
        fullData=[fullData; tempData];
    end
    
        % adjust the time variable
    time=time-time(1);
    
    %% load Cm control data
for i=1:length(photodir)
    cd(photodir{i});
    photofiles=dir('*_dm.mat');
    load(photofiles.name, 'dataTable');
    phototime=cell2mat(struct2cell(load(photofiles.name, 'time')));
    tempData=dataTable;
    
    %add label
    experimentName=repelem(string(photofiles.name(1:end-7)), height(tempData))';
    condition=repelem(string('CmControl'), height(tempData))';
    
    %create new table
    newVar=table(experimentName, condition, 'VariableNames', {'Experiment', 'Condition'});
    photoData=[newVar, tempData];
    
    normMean=mean(photoData.intensity, 1, 'omitnan');
    normStd=std(photoData.intensity, 0, 1, 'omitnan');
    [xD, yD]=prepareCurveData(phototime, normMean);
    photofit=fit(xD, yD, 'exp1');
    photoextrap=photofit.a*exp(photofit.b*time);
end

else
    % Load data
    cd(maindir)
    load('diffusionAnalysis.mat')
end

%% plot the data
% plot(photofit, xD, yD), hold on
% plot(time, photoextrap)
% title('Photobleaching Curve')
% legend({'data', 'fitted curve', 'extrapolation'})
% xlabel('Time (min)')
% ylabel('Intensity (A.U.)')
% cd(maindir)
% saveas(gcf, 'photobleachCurve.fig')
% saveas(gcf, 'photobleachCurve.png')
% 
% dI=[0 abs(photoextrap(2:1:end)-photoextrap(1:1:end-1))];
%% plot data
    
% %first, let's parse the data
% cond1 = fullData(strcmp(fullData.Condition, "LB")==1, :);
% cond2 = fullData(strcmp(fullData.Condition, "EDTA")==1, :);
% cond3 = fullData(strcmp(fullData.Condition, "Mg2+")==1, :);

% if replot==0
%     
%     cd(maindir);
%     
% %% plot a histogram of the contraction lengths
% %is there a difference in lysis dL between halved and non-halved cells
% yData1=cond1(cond1.halfie==0,:);
% yData2=cond1(cond1.halfie==1,:);
% 
% figure(1)
% h1=histogram(yData1.("Overall dL")); hold on
% h2=histogram(yData2.("Overall dL"));
% 
% h1.BinWidth = 0.5;
% h2.BinWidth = 0.5;
% title('Percent Change in Cell Wall Length in Untreated Cells')
% xlabel('Percent Change (%)')
% legend({'not halved', 'halved'})
% saveas(gcf, 'dL_LB.fig')
% saveas(gcf, 'dL_LB.png')
% 
% %in EDTA-treated cells?
% yData1=cond2(cond2.halfie==0,:);
% yData2=cond2(cond2.halfie==1,:);
% 
% figure(2)
% h1=histogram(yData1.("Overall dL")); hold on
% h2=histogram(yData2.("Overall dL"));
% 
% h1.BinWidth = 0.5;
% h2.BinWidth = 0.5;
% title('Percent Change in Cell Wall Length in EDTA-treated Cells')
% xlabel('Percent Change (%)')
% legend({'not halved', 'halved'})
% saveas(gcf, 'dL_EDTA.fig')
% saveas(gcf, 'dL_EDTA.png')
% 
% %in Mg2+-treated cells?
% yData1=cond3(cond3.halfie==0,:);
% yData2=cond3(cond3.halfie==1,:);
% 
% figure(3)
% h1=histogram(yData1.("Overall dL")); hold on
% h2=histogram(yData2.("Overall dL"));
% 
% h1.BinWidth = 0.5;
% h2.BinWidth = 0.5;
% title('Percent Change in Cell Wall Length in Mg^{2+}-treated Cells')
% xlabel('Percent Change (%)')
% legend({'not halved', 'halved'})
% saveas(gcf, 'dL_Mg.fig')
% saveas(gcf, 'dL_Mg.png')
% 
% %% plotting only non-halved cells
% %is there a difference in the lysis dL of non-halved LB vs EDTA vs Mg2+ treated
% %cells? (THERE SHOULDN'T BE)
% yData1=cond1(cond1.halfie==0,:);
% yData2=cond2(cond2.halfie==0,:);
% yData3=cond3(cond3.halfie==0,:);
% 
% figure(4)
% h1=histogram(yData1.("Lysis dL")); hold on
% h2=histogram(yData2.("Lysis dL"));
% h3=histogram(yData3.("Lysis dL"));
% 
% h1.BinWidth = 0.5;
% h2.BinWidth = 0.5;
% h3.BinWidth = 0.5;
% title('Percent Change in Cell Wall Length Upon Lysis')
% xlabel('Percent Change (%)')
% legend({'LB', 'EDTA', 'Mg^{2+}'})
% saveas(gcf, 'lysisdL.fig')
% saveas(gcf, 'lysisdL.png')
% 
% %is there a difference in the post-lysis dL of non-halved LB vs EDTA treated
% %cells? 
% figure(5)
% h1=histogram(yData1.("Post-Lysis dL")); hold on
% h2=histogram(yData2.("Post-Lysis dL"));
% h3=histogram(yData3.("Post-Lysis dL"));
% 
% h1.BinWidth = 0.5;
% h2.BinWidth = 0.5;
% h3.BinWidth = 0.5;
% title('Percent Change in Cell Wall Length Post-Lysis')
% xlabel('Percent Change (%)')
% legend({'LB', 'EDTA', 'Mg^{2+}'})
% saveas(gcf, 'postlysisdL.fig')
% saveas(gcf, 'postlysisdL.png')
% 
% %is there a difference in the overall dL of non-halved LB vs EDTA treated
% %cells? (so far, it seems not)
% figure(6)
% h1=histogram(yData1.("Overall dL")); hold on
% h2=histogram(yData2.("Overall dL"));
% h3=histogram(yData3.("Overall dL"));
% 
% h1.BinWidth = 0.5;
% h2.BinWidth = 0.5;
% h3.BinWidth = 0.5;
% title('Overall Percent Change in Cell Wall Length')
% xlabel('Percent Change (%)')
% legend({'LB', 'EDTA', 'Mg^{2+}'})
% saveas(gcf, 'overalldL.fig')
% saveas(gcf, 'overalldL.png')
% 
% figure(7)
% h1=histogram(yData1.("cell length")(:, tidx+3)); hold on
% h2=histogram(yData2.("cell length")(:, tidx+3));
% h3=histogram(yData3.("cell length")(:, tidx+3)); 
% h1.BinWidth = 0.5;
% h2.BinWidth = 0.5;
% h3.BinWidth = 0.5;
% title('Cell Length Distribution')
% xlabel('Cell Length \mum')
% legend({'LB', 'EDTA', 'Mg^{2+}'})
% saveas(gcf, 'lengthHist.fig')
% saveas(gcf, 'lengthHist.png')
% end
% %% comparing intensities
% %is there a difference in normalized intensity?
% yData1=cond1(cond1.halfie==0,:);
% yData2=cond2(cond2.halfie==0,:);
% yData3=cond3(cond3.halfie==0,:);
% 
% figure(8)
% p1=plot(time, yData1.("normalized intensity"), '-r'); hold on
% p2=plot(time, yData2.("normalized intensity"), '-b'); 
% p3=plot(time, yData3.("normalized intensity"), '-c'); 
% title('Normalized Intensity vs Time')
% txt = ['{\color{red} red = LB, \color{blue} blue = EDTA, \color{cyan} cyan = Mg2+}'];
% subtitle(txt)
% xlabel('Time (min)')
% ylabel('Normalized Intensity (A.U.)')
% saveas(gcf, 'normIntensity.fig')
% saveas(gcf, 'normIntensity.png')
% 
% figure(9)
% p1=plot(time, yData1.intensity, '-r'); hold on
% p2=plot(time, yData2.intensity, '-b'); 
% p3=plot(time, yData3.intensity, '-c'); 
% title('Cellular Intensity vs Time')
% txt = ['{\color{red} red = LB, \color{blue} blue = EDTA, \color{cyan} cyan = Mg2+}'];
% subtitle(txt)
% xlabel('Time (min)')
% ylabel('Intensity (A.U.)')
% saveas(gcf, 'Intensity.fig')
% saveas(gcf, 'Intensity.png')
% 
% %re-plot with smoothed traces
% yData1a=movingaverage(yData1.("normalized intensity"), 5);
% yData2a=movingaverage(yData2.("normalized intensity"), 5);
% yData3a=movingaverage(yData3.("normalized intensity"), 5);
% 
% yData1b=movingaverage(yData1.intensity, 5);
% yData2b=movingaverage(yData2.intensity, 5);
% yData3b=movingaverage(yData3.intensity, 5);
% 
% figure(10)
% p1=plot(time, yData1a, '-r'); hold on
% p2=plot(time, yData2a, '-b'); 
% p3=plot(time, yData3a, '-c'); 
% title('Normalized Intensity vs Time')
% txt = ['{\color{red} red = LB, \color{blue} blue = EDTA, \color{cyan} cyan = Mg2+}'];
% subtitle(txt)
% xlabel('Time (min)')
% ylabel('Normalized Intensity (A.U.)')
% saveas(gcf, 'normIntensitySmooth.fig')
% saveas(gcf, 'normIntensitySmooth.png')
% 
% figure(11)
% p1=plot(time, yData1b, '-r'); hold on
% p2=plot(time, yData2b, '-b'); 
% p3=plot(time, yData3b, '-c'); 
% title('Cellular Intensity vs Time')
% txt = ['{\color{red} red = LB, \color{blue} blue = EDTA, \color{cyan} cyan = Mg2+}'];
% subtitle(txt)
% xlabel('Time (min)')
% ylabel('Intensity (A.U.)')
% saveas(gcf, 'IntensitySmooth.fig')
% saveas(gcf, 'IntensitySmooth.png')
% 
% %% Fit the data
% %re-subset the data
% yData1=cond1(cond1.halfie==0,:);
% yData2=cond2(cond2.halfie==0,:);
% yData3=cond3(cond3.halfie==0,:);
% 
% %pre-allocate variables
% coeff1=nan(height(yData1),2);
% coeff2=nan(height(yData2),2);
% coeff3=nan(height(yData3),2);
% 
% f1=cell(height(yData1),1);
% f2=cell(height(yData2),1);
% f3=cell(height(yData3),1);
% 
% thalf1=nan(height(yData1),1);
% thalf2=nan(height(yData2),1);
% thalf3=nan(height(yData3),1);
% 
% r1=nan(height(yData1),1);
% r2=nan(height(yData2),1);
% r3=nan(height(yData3),1);
% 
% for n=1:height(yData1)
%     
%     yData=yData1.("normalized intensity")(n,:);
%     [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%     f1{n,1}=fit(xData, yData, 'exp1');
%     coeff1(n,:)=coeffvalues(f1{n});
%     
%     %yData=movingaverage(yData,5);
%     yData=yData';
%     cutoff=yData(1,1)/2;
%     
%     if yData(1,end)<=cutoff
%         csum=cumsum(yData<=cutoff);
%         idx=min(find(csum==1));
%         r1(n)=cutoff-yData(idx);
%         thalf1(n)=time(idx);
%     else
%         [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%         f=fit(xData, yData, 'poly1');
%         thalf1(n)=(cutoff-f.p2)/f.p1;
%     end
% end
% 
% for n=1:height(yData2)
%     
%     yData=yData2.("normalized intensity")(n,:);
%     [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%     f2{n,1}=fit(xData, yData, 'exp1');
%     coeff2(n,:)=coeffvalues(f2{n});
%     
%     %yData=movingaverage(yData,5);
%     yData=yData';
%     cutoff=yData(1,1)/2;
%     
%     if yData(1,end)<=cutoff
%         csum=cumsum(yData<=cutoff);
%         idx=min(find(csum==1));
%         r2(n)=cutoff-yData(idx);
%         thalf2(n)=time(idx);
%     else
%         [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%         f=fit(xData, yData, 'poly1');
%         thalf2(n)=(cutoff-f.p2)/f.p1;
%     end
% end
% 
% for n=1:height(yData3)
%     
%     yData=yData3.("normalized intensity")(n,:);
%     [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%     f3{n,1}=fit(xData, yData, 'exp1');
%     coeff3(n,:)=coeffvalues(f3{n});
%     
%     %yData=movingaverage(yData,5);
%     yData=yData';
%     cutoff=yData(1,1)/2;
%     
%     if yData(1,end)<=cutoff
%         csum=cumsum(yData<=cutoff);
%         idx=min(find(csum==1));
%         r3(n)=cutoff-yData(idx);
%         thalf3(n)=time(idx);
%     else
%         [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%         f=fit(xData, yData, 'poly1');
%         thalf3(n)=(cutoff-f.p2)/f.p1;
%     end
% end
% 
% figure(12)
% h1=histogram(-1./coeff1(:,2)); hold on
% h2=histogram(-1./coeff2(:,2));
% h3=histogram(-1./coeff3(:,2)); 
% %h1.Normalization = 'probability';
% %h1.BinWidth = 0.25;
% %h2.Normalization = 'probability';
% %h2.BinWidth = 0.25;
% h2.NumBins = 12;
% %h3.Normalization = 'probability';
% title('Time Constant Distribution')
% xlabel('Tau (-min)')
% legend({'LB', 'EDTA', 'Mg^{2+}'})
% saveas(gcf, 'tauHist.fig')
% saveas(gcf, 'tauHist.png')
% 
% figure(13)
% h1=histogram(thalf1); hold on
% h2=histogram(thalf2);
% h3=histogram(thalf3); 
% %h1.Normalization = 'probability';
% %h1.BinWidth = 0.25;
% %h2.Normalization = 'probability';
% %h2.BinWidth = 0.25;
% h2.NumBins = 12;
% %h3.Normalization = 'probability';
% %h3.BinWidth = 0.25;
% title('t_{1/2} Distribution')
% xlabel('t_{1/2} (-min)')
% legend({'LB', 'EDTA', 'Mg^{2+}'})
% saveas(gcf, 'thalfHist.fig')
% saveas(gcf, 'thalfHist.png')
% 
% %% now adjust for photobleaching
% %is there a difference in normalized intensity?
% yData1=cond1(cond1.halfie==0,:);
% yData2=cond2(cond2.halfie==0,:);
% yData3=cond3(cond3.halfie==0,:);
% 
% yData1=yData1.intensity + dI;
% yData2=yData2.intensity + dI;
% yData3=yData3.intensity + dI;
% 
% figure(14)
% p1=plot(time, yData1, '-r'); hold on
% p2=plot(time, yData2, '-b'); 
% p3=plot(time, yData3, '-c'); 
% title('Intensity vs Time, photoadjusted')
% txt = ['{\color{red} red = LB, \color{blue} blue = EDTA, \color{cyan} cyan = Mg2+}'];
% subtitle(txt)
% xlabel('Time (min)')
% ylabel('Intensity (A.U.)')
% saveas(gcf, 'IntensityAdj.fig')
% saveas(gcf, 'IntensityAdj.png')
% 
% yData1=yData1./yData1(:,1);
% yData2=yData2./yData2(:,1);
% yData3=yData3./yData3(:,1);
% 
% figure(15)
% p1=plot(time, yData1, '-r'); hold on
% p2=plot(time, yData2, '-b'); 
% p3=plot(time, yData3, '-c'); 
% title('Normalized Intensity vs Time, photoadjusted')
% txt = ['{\color{red} red = LB, \color{blue} blue = EDTA, \color{cyan} cyan = Mg2+}'];
% subtitle(txt)
% xlabel('Time (min)')
% ylabel('Normalized Intensity (A.U.)')
% saveas(gcf, 'NormIntensityAdj.fig')
% saveas(gcf, 'NormIntensityAdj.png')
% 
% %pre-allocate variables
% coeff1a=nan(height(yData1),2);
% coeff2a=nan(height(yData2),2);
% coeff3a=nan(height(yData3),2);
% 
% f1a=cell(height(yData1),1);
% f2a=cell(height(yData2),1);
% f3a=cell(height(yData3),1);
% 
% thalf1a=nan(height(yData1),1);
% thalf2a=nan(height(yData2),1);
% thalf3a=nan(height(yData3),1);
% 
% r1a=nan(height(yData1),1);
% r2a=nan(height(yData2),1);
% r3a=nan(height(yData3),1);
% 
% for n=1:height(yData1)
%     
%     yData=yData1(n,:);
%     [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%     f1a{n,1}=fit(xData, yData, 'exp1');
%     coeff1a(n,:)=coeffvalues(f1a{n});
%     
%     %yData=movingaverage(yData,5);
%     yData=yData';
%     cutoff=yData(1,1)/2;
%     
%     if yData(1,end)<=cutoff
%         csum=cumsum(yData<=cutoff);
%         idx=min(find(csum==1));
%         r1a(n)=cutoff-yData(idx);
%         thalf1a(n)=time(idx);
%     else
%         [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%         f=fit(xData, yData, 'poly1');
%         thalf1a(n)=(cutoff-f.p2)/f.p1;
%     end
% end
% 
% for n=1:height(yData2)
%     n
%     yData=yData2(n,:);
%     [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%     f2a{n,1}=fit(xData, yData, 'exp1');
%     coeff2a(n,:)=coeffvalues(f2a{n});
%     
%     %yData=movingaverage(yData,5);
%     yData=yData';
%     cutoff=yData(1,1)/2;
%     
%     if yData(1,end)<=cutoff
%         csum=cumsum(yData<=cutoff);
%         idx=min(find(csum==1));
%         r2a(n)=cutoff-yData(idx);
%         thalf2a(n)=time(idx);
%     else
%         [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%         f=fit(xData, yData, 'poly1');
%         thalf2a(n)=(cutoff-f.p2)/f.p1;
%     end
% end
% 
% for n=1:height(yData3)
%     
%     yData=yData3(n,:);
%     [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%     f3a{n,1}=fit(xData, yData, 'exp1');
%     coeff3a(n,:)=coeffvalues(f3a{n});
%     
%     %yData=movingaverage(yData,5);
%     yData=yData';
%     cutoff=yData(1,1)/2;
%     
%     if yData(1,end)<=cutoff
%         csum=cumsum(yData<=cutoff);
%         idx=min(find(csum==1));
%         r3a(n)=cutoff-yData(idx);
%         thalf3a(n)=time(idx);
%     else
%         [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
%         f=fit(xData, yData, 'poly1');
%         thalf3a(n)=(cutoff-f.p2)/f.p1;
%     end
% end
% 
% figure(16)
% h1=histogram(-1./coeff1a(:,2)); hold on
% h2=histogram(-1./coeff2a(:,2));
% h3=histogram(-1./coeff3a(:,2)); 
% %h1.Normalization = 'probability';
% %h1.BinWidth = 0.25;
% %h2.Normalization = 'probability';
% %h2.BinWidth = 0.25;
% h1.NumBins = 12;
% h2.NumBins = 12;
% h3.NumBins = 12;
% %h3.Normalization = 'probability';
% %h3.BinWidth = 0.25;
% title('Time Constant Distribution, Photobleach Adjusted')
% xlabel('Tau (-min)')
% legend({'LB', 'EDTA', 'Mg^{2+}'})
% saveas(gcf, 'tauHistAdj.fig')
% saveas(gcf, 'tauHistAdj.png')
% 
% figure(17)
% h1=histogram(thalf1a); hold on
% h2=histogram(thalf2a);
% h3=histogram(thalf3a); 
% %h1.Normalization = 'probability';
% %h1.BinWidth = 0.25;
% %h2.Normalization = 'probability';
% %h2.BinWidth = 0.25;
% h1.NumBins = 12;
% h2.NumBins = 12;
% h3.NumBins = 12;
% %h3.Normalization = 'probability';
% %h3.BinWidth = 0.25;
% title('t_{1/2} Distribution, Photobleach Adjusted')
% xlabel('t_{1/2} (-min)')
% legend({'LB', 'EDTA', 'Mg^{2+}'})
% saveas(gcf, 'thalfHistAdj.fig')
% saveas(gcf, 'thalfHistAdj.png')

%% determine the pattern of the data
pattern1=[];
pattern2=[];
pattern3=[];

yData1=cond1(cond1.halfie==0,:);
yData2=cond2(cond2.halfie==0,:);
yData3=cond3(cond3.halfie==0,:);

yData1c=yData1.intensity;
yData2c=yData2.intensity;
yData3c=yData3.intensity;

figure
for n=1:height(yData1)
    plot(time, yData1c(n,:))
    prompt = 'What pattern fits the data? 0=none, 1=linear, 2=sigmoidal, 3=exponential decay, 4=hyperbolic: ';
    answer = input(prompt)
    pattern1=[pattern1; answer];
    close
end

%create a table from these new variables
newVar1=table(pattern1, 'VariableNames', {'Pattern'});
yData1=[yData1, newVar1];
%% save the data
%cd(maindir)
%save('diffusionAnalysis.mat')

%close all
%% Functions
% function [y]=sigmoidal(b,x)
%     xData=x.*b(1);
%     c=b(1)*b(2);
%     denom=1+exp(xData-c);
%     y=1./denom;
% end