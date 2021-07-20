%diffusionAnalysis.m
%Author: Zarina Akbary
%Date: 17 July 2021
%purpose: to analyze mNeonGreen.m diffusion out of the cell

clear, close all
%% User Input
dirname=['/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis'];
recrunch=0;
replot=1;
%%%%%%%%%%%%%%%%%%%

%% Part I: New Tables
if recrunch==0
    % Load data
    cd(dirname)
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
else
    % Load data
    cd(dirname)
    load('diffusionAnalysis')
end

%% adjust the time variable
time=time-time(1);

%% plot data
    
%first, let's parse the data
cond1 = fullData(strcmp(fullData.Condition, "LB")==1, :);
cond2 = fullData(strcmp(fullData.Condition, "EDTA")==1, :);
cond3 = fullData(strcmp(fullData.Condition, "Mg2+")==1, :);

if replot==0
    
    cd(savedir);
    
%% plot a histogram of the contraction lengths
    %is there a difference in lysis dL between halved and non-halved cells
    yData1=cond1(cond1.halfie==0,:);
    yData2=cond1(cond1.halfie==1,:);

    figure(1)
    h1=histogram(yData1.("Overall dL")); hold on
    h2=histogram(yData2.("Overall dL"));

    h1.BinWidth = 0.5;
    h2.BinWidth = 0.5;
    title('Percent Change in Cell Wall Length in Untreated Cells')
    xlabel('Percent Change (%)')
    legend({'not halved', 'halved'})
    saveas(gcf, 'dL_LB.fig')
    saveas(gcf, 'dL_LB.png')

    %in EDTA-treated cells?
    yData1=cond2(cond2.halfie==0,:);
    yData2=cond2(cond2.halfie==1,:);

    figure(2)
    h1=histogram(yData1.("Overall dL")); hold on
    h2=histogram(yData2.("Overall dL"));

    h1.BinWidth = 0.5;
    h2.BinWidth = 0.5;
    title('Percent Change in Cell Wall Length in EDTA-treated Cells')
    xlabel('Percent Change (%)')
    legend({'not halved', 'halved'})
    saveas(gcf, 'dL_EDTA.fig')
    saveas(gcf, 'dL_EDTA.png')

    %in Mg2+-treated cells?
    yData1=cond3(cond3.halfie==0,:);
    yData2=cond3(cond3.halfie==1,:);

    figure(3)
    h1=histogram(yData1.("Overall dL")); hold on
    h2=histogram(yData2.("Overall dL"));

    h1.BinWidth = 0.5;
    h2.BinWidth = 0.5;
    title('Percent Change in Cell Wall Length in Mg^{2+}-treated Cells')
    xlabel('Percent Change (%)')
    legend({'not halved', 'halved'})
    saveas(gcf, 'dL_Mg.fig')
    saveas(gcf, 'dL_Mg.png')
    
%% plotting only non-halved cells
    %is there a difference in the lysis dL of non-halved LB vs EDTA vs Mg2+ treated
    %cells? (THERE SHOULDN'T BE)
    yData1=cond1(cond1.halfie==0,:);
    yData2=cond2(cond2.halfie==0,:);
    yData3=cond3(cond3.halfie==0,:);

    figure(4)
    h1=histogram(yData1.("Lysis dL")); hold on
    h2=histogram(yData2.("Lysis dL"));
    h3=histogram(yData3.("Lysis dL"));

    h1.BinWidth = 0.5;
    h2.BinWidth = 0.5;
    h3.BinWidth = 0.5;
    title('Percent Change in Cell Wall Length Upon Lysis')
    xlabel('Percent Change (%)')
    legend({'LB', 'EDTA', 'Mg^{2+}'})
    saveas(gcf, 'lysisdL.fig')
    saveas(gcf, 'lysisdL.png')

    %is there a difference in the post-lysis dL of non-halved LB vs EDTA treated
    %cells? 
    figure(5)
    h1=histogram(yData1.("Post-Lysis dL")); hold on
    h2=histogram(yData2.("Post-Lysis dL"));
    h3=histogram(yData3.("Post-Lysis dL"));

    h1.BinWidth = 0.5;
    h2.BinWidth = 0.5;
    h3.BinWidth = 0.5;
    title('Percent Change in Cell Wall Length Post-Lysis')
    xlabel('Percent Change (%)')
    legend({'LB', 'EDTA', 'Mg^{2+}'})
    saveas(gcf, 'postlysisdL.fig')
    saveas(gcf, 'postlysisdL.png')

    %is there a difference in the overall dL of non-halved LB vs EDTA treated
    %cells? (so far, it seems not)
    figure(6)
    h1=histogram(yData1.("Overall dL")); hold on
    h2=histogram(yData2.("Overall dL"));
    h3=histogram(yData3.("Overall dL"));

    h1.BinWidth = 0.5;
    h2.BinWidth = 0.5;
    h3.BinWidth = 0.5;
    title('Overall Percent Change in Cell Wall Length')
    xlabel('Percent Change (%)')
    legend({'LB', 'EDTA', 'Mg^{2+}'})
    saveas(gcf, 'overalldL.fig')
    saveas(gcf, 'overalldL.png')
    
    figure(7)
    h1=histogram(yData1.("cell length")(:, tidx+3)); hold on
    h2=histogram(yData2.("cell length")(:, tidx+3));
    h3=histogram(yData3.("cell length")(:, tidx+3)); 
    h1.BinWidth = 0.5;
    h2.BinWidth = 0.5;
    h3.BinWidth = 0.5;
    title('Cell Length Distribution')
    xlabel('Cell Length \mum')
    legend({'LB', 'EDTA', 'Mg^{2+}'})
    saveas(gcf, 'lengthHist.fig')
    saveas(gcf, 'lengthHist.png')
    
%% comparing intensities
    %is there a difference in normalized intensity?
    yData1=cond1(cond1.halfie==0,:);
    yData2=cond2(cond2.halfie==0,:);
    yData3=cond3(cond3.halfie==0,:);

    figure(8)
    p1=plot(time, yData1.("normalized intensity"), '-r'); hold on
    p2=plot(time, yData2.("normalized intensity"), '-b'); 
    p3=plot(time, yData3.("normalized intensity"), '-c'); 
    title('Normalized Intensity vs Time')
    txt = ['{\color{red} red = LB, \color{blue} blue = EDTA, \color{cyan} cyan = Mg2+}'];
    subtitle(txt)
    xlabel('Time (min)')
    ylabel('Normalized Intensity (A.U.)')
    saveas(gcf, 'normIntensity.fig')
    saveas(gcf, 'normIntensity.png')

    figure(9)
    p1=plot(time, yData1.intensity, '-r'); hold on
    p2=plot(time, yData2.intensity, '-b'); 
    p3=plot(time, yData3.intensity, '-c'); 
    title('Cellular Intensity vs Time')
    txt = ['{\color{red} red = LB, \color{blue} blue = EDTA, \color{cyan} cyan = Mg2+}'];
    subtitle(txt)
    xlabel('Time (min)')
    ylabel('Intensity (A.U.)')
    saveas(gcf, 'normIntensity.fig')
    saveas(gcf, 'normIntensity.png')

%% Fit the data
cutoff=0.4;

%for LB
yData1a=cond1(cond1.halfie==0 & cond1.("normalized intensity")(:,end)>cutoff,:);
yData1b=cond1(cond1.halfie==0 & cond1.("normalized intensity")(:, end)<cutoff,:);

figure, hold on
for n=1:height(yData1)
    plot(time, yData1(n,:))
end

%for EDTA
yData2=cond2(cond2.halfie==0,:);

%for Mg2+
yData3=cond3(cond3.halfie==0,:);

%for time
xData=time;

%pre-allocate variables
coeff0=[0,0];
% coeff1=cell(height(yData1), 2);
coeff2=nan(height(yData2), 2);
coeff3=nan(height(yData3), 2);

f1=cell(height(yData1),1);

yhat1=nan(height(yData1), length(time));
yhat2=nan(height(yData2), length(time));
yhat3=nan(height(yData3), length(time));

%create custom fit object
sig=fittype('1/(1+a*exp((x-b)/c))', 'independent', {'x'}, 'coefficients', {'a', 'b', 'c'});
%for LB
for n=1:height(yData1)
    n
    [xData, yData]=prepareCurveData(time, yData1.("normalized intensity")(n,:));
    if yData(end,1)>cutoff
        f1{n,1}=fit(xData, yData, 'poly1');
    else
        f1{n,1}=fit(xData, yData, sig);
    end
    
    plot(f1{n,1}, xData, yData), pause, close
end

%for EDTA
for n=1:height(yData2)
    coeff2(n,:)=nlinfit(xData, yData2.("normalized intensity")(n,:), @sigmoidal, coeff0);
    yhat2(n,:)=sigmoidal(coeff2(n,:), xData);
end

%for Mg2+
for n=1:height(yData3)
    coeff3(n,:)=nlinfit(xData, yData3.("normalized intensity")(n,:),@sigmoidal, coeff0);
    yhat3(n,:)=sigmoidal(coeff3(n,:), xData);
end

for n=1:height(yData1)
    n
    figure
    plot(xData, yhat1(n,:), '-', 'Color', '#D95319'), hold on
    scatter(xData, yData1.("normalized intensity")(n,:), '.', 'b')
    title('Normalized Intensity vs Time Fit')
    xlabel('Time (min)')
    ylabel('Normalized Intensity (A.U.)')
    pause, close
end

figure(9), hold on
for n=1:height(yData2)
    plot(xData, yhat2(n,:), '-', 'Color', '#D95319')
    scatter(xData, yData2.("normalized intensity")(n,:), '.', 'b')
    title('Normalized Intensity vs Time Fit')
    xlabel('Time (min)')
    ylabel('Normalized Intensity (A.U.)')
end

for n=1:height(yData3)
    figure
    plot(xData, yhat3(n,:), '-', 'Color', '#D95319'), hold on
    scatter(xData, yData3.("normalized intensity")(n,:), '.', 'b')
    title('Normalized Intensity vs Time Fit')
    xlabel('Time (min)')
    ylabel('Normalized Intensity (A.U.)')
    pause, close
end
else
%% find t_{1/2}
yData1=cond1(cond1.halfie==0,:);
yData1=yData1.("normalized intensity");
yData2=cond2(cond2.halfie==0,:);
yData2=yData2.("normalized intensity");
yData3=cond3(cond3.halfie==0,:);
yData3=yData3.("normalized intensity");

%first, smooth the data by calculating a movingaverage
yData1=movingaverage(yData1, 5);
yData2=movingaverage(yData2, 5);
yData3=movingaverage(yData3, 5);

%pre-allocate variables
tconst=[];
r1=[];
expon=0;
%now, go through each trace and determine t1/2
for n=1:height(yData1)
    if yData1(n,end)<=yData1(n,1)/2
%         idx=find(yData1<=(yData1(n,1)/2), 'last');
%         r1=[r1; (yData1(n,1)/2-yData1(idx))];
%         tconst=[tconst time(idx)];
        f=fit(time', yData1(n,:)', 'exp1')
        plot(f, time, yData1(n,:))
        title('Where intensity ends with less than half the original')
        pause, close
        expon=expon+1;
        
    else
%         f=fit(time, yData1(n,:), 'exp2')
%         plot(time, yData1(n,:))
%         title('Where intensity ends with more than half the original')
%         pause, close
    end
end
end
%now, calculate the dy/dx
% d1=diff(yData1, 1,2);
% d2=diff(yData2, 1,2);
% d3=diff(yData3, 1,2);
% 
% cd(dirname)
% save('diffusionAnalysis.mat')

%% Functions
function [y]=sigmoidal(b,x)
    xData=x.*b(1);
    c=b(1)*b(2);
    denom=1+exp(xData-c);
    y=1./denom;
end