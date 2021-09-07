%diffusionAnalysis.m
%Author: Zarina Akbary
%Date: 17 July 2021
%purpose: to analyze mNeonGreen.m diffusion out of the cell

clear, close all
%% User Input
maindir=['/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis'];
photodir={['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07242021_analysis/07242021_Exp1_colony1/07242021_Exp1_mNeonGreen/07242021_Exp1_figures'];['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07242021_analysis/07242021_Exp1_colony3/07242021_Exp1_mNeonGreen/07242021_Exp1_figures'];['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07242021_analysis/07242021_Exp1_colony4/07242021_Exp1_mNeonGreen/07242021_Exp1_figures'];['/Users/zarina/Downloads/NYU/Year3_2021_Summer/07242021_analysis/07242021_Exp1_colony6/07242021_Exp1_mNeonGreen/07242021_Exp1_figures']};
recrunch=0;
replot=1;
%%%%%%%%%%%%%%%%%%%

%% Part I: New Tables
if recrunch==0
    % Load data
    cd(maindir)
    dataDirectory=dir('*_dm.mat');

    fullData=table();
    pbsData=table();
    
    for f=1:height(dataDirectory)

        load(dataDirectory(f).name)
        tempData = dataTable;
        
        %label experiment
        experimentName=repelem(string(dataDirectory(f).name(1:end-7)), height(tempData))';

        if f<=3
            condition=repelem(string('Mg2+'), height(tempData))';
        elseif f>3 & f<6
            condition=repelem(string('LB'), height(tempData))';
        elseif f>=6 & f<=9
            condition=repelem(string('EDTA'), height(tempData))';
        elseif f>9
            condition=repelem(string('PBS'), height(tempData))';
        end
        
        if f>9
            pbsTime=time;
        else
            fullTime=time;
        end
        
        if f<=9
            tG=6; %frame in which cells are still profused with LB
            tL=7; %frame in which cells lyse
            tL2=8; %frame with detergent
            tC=10:13; %frame in condition
        elseif f>9
            tG=13; %frame in which cells are still profused with LB
            tL=16; %frame in which cells lyse
            tL2=16; %frame with detergent
            tC=20:23; %frame in condition
        end
        
        %calculate cell length changes
        %from LB > LB + 5% detergent is t=6 and t=7
        dL=((lcell(:, tL)-lcell(:, tG))./lcell(:,tG))*100; %lf-li/li
        
        %from LB + 5% detergent to the condition 
        tf=mean(lcell(:, tC), 2, 'omitnan');
        dL2=((tf-lcell(:,tL2))./lcell(:,tL2))*100; %lf-li/li
        %dL2=((lcell(:, 8)-lcell(:, 10))./lcell(:,8))*100; %lf-li/li
        
        %from LB to the condition 
        dL3=((tf-lcell(:,tG))./lcell(:,tG))*100; %lf-li/li
        
        if f<=9
            
            %create a table from these new variables
            newVar=table(experimentName, condition, dL, dL2,dL3, 'VariableNames', {'Experiment', 'Condition', 'Lysis dL', 'Post-Lysis dL', 'Overall dL'});
            tempData=[newVar, tempData];
            fullData=[fullData; tempData];
            
        elseif f>9
            %create a table from these new variables
            newVar=table(experimentName, condition, dL, dL2,dL3, 'VariableNames', {'Experiment', 'Condition', 'Lysis dL', 'Post-Lysis dL', 'Overall dL'});
            tempData=[newVar, tempData];
            pbsData=[pbsData; tempData];
        end
        
    end
    
 % adjust the time variable
 time=fullTime-fullTime(1);
    
%% Part II: load control data
photoData=table();
for i=1:length(photodir)
    
    cd(photodir{i});
    photofiles=dir('*_dm.mat');
    load(photofiles.name, 'dataTable');
    phototime=cell2mat(struct2cell(load(photofiles.name, 'time')));
    tempData=dataTable;
    
    %add label
    experimentName=repelem(string(photofiles.name(1:end-7)), height(tempData))';
    condition=repelem(string('AzideControl'), height(tempData))';
    
    %create new table
    newVar=table(experimentName, condition, 'VariableNames', {'Experiment', 'Condition'});
    tempData=[newVar, tempData];
    
   photoData=[photoData; tempData];
    
end

normMean=mean(photoData.intensity, 1, 'omitnan');
normStd=std(photoData.intensity, 0, 1, 'omitnan');
[xD, yD]=prepareCurveData(phototime, normMean);
photofit=fit(xD, yD, 'poly1');

else
    % Load data
    cd(maindir)
    load('diffusionAnalysis.mat')
end

%% Part III: plot the photobleaching data
figure
plot(photofit, xD, yD)
title('Photobleaching Fit')
xlabel('Time (min)')
ylabel('Intensity (A.U.)')

cd(maindir)
saveas(gcf,'photobleachFit.png')
saveas(gcf,'photobleachFit.fig')

close

figure
plot(phototime, photoData.intensity)
title('Photobleaching Traces')
xlabel('Time (min)')
ylabel('Intensity (A.U.)')

cd(maindir)
saveas(gcf,'photobleachTraces.png')
saveas(gcf,'photobleachTraces.fig')

close
%% Part IV: adjust for photobleaching
    
% %first, let's parse the data
cond1 = fullData(strcmp(fullData.Condition, "LB")==1, :);
cond2 = fullData(strcmp(fullData.Condition, "EDTA")==1, :);
cond3 = fullData(strcmp(fullData.Condition, "Mg2+")==1, :);

yData1=cond1(cond1.halfie==0,:);
yData2=cond2(cond2.halfie==0,:);
yData3=cond3(cond3.halfie==0,:);
yData4=pbsData(pbsData.halfie==0, :);

yadj=[];
for n=1:height(yData1)
    yData=yData1.intensity(n,:);
    xData=time;
    yexp=(photofit.p1.*xData)+max(yData);
    dy=[0, diff(yexp)];
    yadj=[yadj; yData+abs(dy)];
end
newVar=table(yadj, yadj./yadj(:,1), 'VariableNames', {'intensityAdj', 'normalAdj'});
yData1=[yData1, newVar];

yadj=[];
for n=1:height(yData2)
    yData=yData2.intensity(n,:);
    xData=time;
    yexp=(photofit.p1.*xData)+max(yData);
    dy=[0, diff(yexp)];
    yadj=[yadj; yData+abs(dy)];
end
newVar=table(yadj, yadj./yadj(:,1), 'VariableNames', {'intensityAdj', 'normalAdj'});
yData2=[yData2, newVar];

yadj=[];
for n=1:height(yData3)
    yData=yData3.intensity(n,:);
    xData=time;
    yexp=(photofit.p1.*xData)+max(yData);
    dy=[0, diff(yexp)];
    yadj=[yadj; yData+abs(dy)];
    
    if n==1
        plot(xData, yexp), hold on
        plot(xData, yData),
        plot(xData, yData+abs(dy))
        xlabel('Time (min)')
        ylabel('Intensity (A.U.)')
        legend({'Fit', 'Data', 'Adjusted Data'})
        title('Intensity vs Time, Mg^{2+}')
        cd(maindir)
        saveas(gcf, 'mgPhotoFit.png')
        saveas(gcf, 'mgPhotoFit.fig')
        pause, close
    end
end
newVar=table(yadj, yadj./yadj(:,1), 'VariableNames', {'intensityAdj', 'normalAdj'});
yData3=[yData3, newVar];

yadj=[];
for n=1:height(yData4)
    yData=yData4.intensity(n,:);
    xData=pbsTime;
    yexp=(photofit.p1.*xData)+max(yData);
    dy=[0, diff(yexp)];
    yadj=[yadj; yData+abs(dy)];
    
    if n==1
        plot(xData, yexp), hold on
        plot(xData, yData),
        plot(xData, yData+abs(dy))
        legend({'Fit', 'Data', 'Adjusted Data'})
        xlabel('Time (min)')
        ylabel('Intensity (A.U.)')
        title('Intensity vs Time, *PBS')
        cd(maindir)
        saveas(gcf, 'pbsPhotoFit.png')
        saveas(gcf, 'pbsPhotoFit.fig')
        pause, close
    end
end
newVar=table(yadj, yadj./yadj(:,1), 'VariableNames', {'intensityAdj', 'normalAdj'});
yData4=[yData4, newVar];

%% Part V: extract time constants
% %pre-allocate variables
coeff1=nan(height(yData1),2);
coeff2=nan(height(yData2),2);
coeff3=nan(height(yData3),2);
coeff4=nan(height(yData4),2);

f1=cell(height(yData1),1);
f2=cell(height(yData2),1);
f3=cell(height(yData3),1);
f4=cell(height(yData4),1);

thalf1=nan(height(yData1),1);
thalf2=nan(height(yData2),1);
thalf3=nan(height(yData3),1);
thalf4=nan(height(yData4),1);

r1=nan(height(yData1),1);
r2=nan(height(yData2),1);
r3=nan(height(yData3),1);
r4=nan(height(yData4),1);

for n=1:height(yData1)
    
    yData=yData1.normalAdj(n,:);
    [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
    f1{n,1}=fit(xData, yData, 'exp1');
    coeff1(n,:)=coeffvalues(f1{n});

    yData=yData';
    cutoff=yData(1,1)/2;
    
    if yData(1,end)<=cutoff
        csum=cumsum(yData<=cutoff);
        idx=min(find(csum==1));
        r1(n)=cutoff-yData(idx);
        thalf1(n)=time(idx);
    else
        [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
        f=fit(xData, yData, 'poly1');
        thalf1(n)=(cutoff-f.p2)/f.p1;
    end
end

for n=1:height(yData2)
    
    yData=yData2.normalAdj(n,:);
    [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
    f2{n,1}=fit(xData, yData, 'exp1');
    coeff2(n,:)=coeffvalues(f2{n});

    yData=yData';
    cutoff=yData(1,1)/2;
    
    if yData(1,end)<=cutoff
        csum=cumsum(yData<=cutoff);
        idx=min(find(csum==1));
        r2(n)=cutoff-yData(idx);
        thalf2(n)=time(idx);
    else
        [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
        f=fit(xData, yData, 'poly1');
        thalf2(n)=(cutoff-f.p2)/f.p1;
    end
end

for n=1:height(yData3)
    
    yData=yData3.normalAdj(n,:);
    [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
    f3{n,1}=fit(xData, yData, 'exp1');
    coeff3(n,:)=coeffvalues(f3{n});

    yData=yData';
    cutoff=yData(1,1)/2;
    
    if yData(1,end)<=cutoff
        csum=cumsum(yData<=cutoff);
        idx=min(find(csum==1));
        r3(n)=cutoff-yData(idx);
        thalf3(n)=time(idx);
    else
        [xData, yData]=prepareCurveData(time(1:length(yData)), yData);
        f=fit(xData, yData, 'poly1');
        thalf3(n)=(cutoff-f.p2)/f.p1;
    end
end

for n=1:height(yData4)
    
    yData=yData4.normalAdj(n,:);
    [xData, yData]=prepareCurveData(pbsTime(1:length(yData)), yData);
    f4{n,1}=fit(xData, yData, 'exp1');
    coeff4(n,:)=coeffvalues(f4{n});

    yData=yData';
    cutoff=yData(1,1)/2;
    
    if yData(1,end)<=cutoff
        csum=cumsum(yData<=cutoff);
        idx=min(find(csum==1));
        r4(n)=cutoff-yData(idx);
        thalf4(n)=time(idx);
    else
        [xData, yData]=prepareCurveData(pbsTime(1:length(yData)), yData);
        f=fit(xData, yData, 'poly1');
        thalf4(n)=(cutoff-f.p2)/f.p1;
    end
end

%% Part VI: plot the data
figure
p1=plot(time, yData1.normalAdj, 'Color', [0, 0.4470, 0.7410]); hold on
p2=plot(time, yData2.normalAdj, 'Color', [0.4940, 0.1840, 0.5560]); 
p3=plot(time, yData3.normalAdj, 'Color', [0.4660, 0.6740, 0.1880]); 
p4=plot(pbsTime, yData4.normalAdj, 'Color', [0.6350, 0.0780, 0.1840]); 
title('Normalized Intensity vs Time, photobleach adjusted')
txt = ['{\color[rgb]{0, 0.4470, 0.7410} blue = LB, \color[rgb]{0.4940, 0.1840, 0.5560} purple = EDTA, \color[rgb]{0.4660, 0.6740, 0.1880} green = Mg2+, \color[rgb]{0.6350, 0.0780, 0.1840} red = *PBS,}'];
subtitle(txt)
xlabel('Time (min)')
ylabel('Normalized Intensity (A.U.)')
close
% saveas(gcf, 'normIntensityAdj.fig')
% saveas(gcf, 'normIntensityAdj.png')

figure
p1=plot(time, yData1.intensityAdj, 'Color', [0, 0.4470, 0.7410]); hold on
p2=plot(time, yData2.intensityAdj, 'Color', [0.4940, 0.1840, 0.5560]); 
p3=plot(time, yData3.intensityAdj, 'Color', [0.4660, 0.6740, 0.1880]); 
p4=plot(pbsTime, yData4.intensityAdj, 'Color', [0.6350, 0.0780, 0.1840]); 
title('Intensity vs Time, photobleach adjusted')
txt = ['{\color[rgb]{0, 0.4470, 0.7410} blue = LB, \color[rgb]{0.4940, 0.1840, 0.5560} purple = EDTA, \color[rgb]{0.4660, 0.6740, 0.1880} green = Mg2+, \color[rgb]{0.6350, 0.0780, 0.1840} red = *PBS,}'];
subtitle(txt)
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
saveas(gcf, 'IntensityAdj.fig')
saveas(gcf, 'IntensityAdj.png')
close

figure
p1=plot(time, yData1.intensity, 'Color', [0, 0.4470, 0.7410]); hold on
p2=plot(time, yData2.intensity, 'Color', [0.4940, 0.1840, 0.5560]); 
p3=plot(time, yData3.intensity, 'Color', [0.4660, 0.6740, 0.1880]); 
p4=plot(pbsTime, yData4.intensity, 'Color', [0.6350, 0.0780, 0.1840]); 
title('Intensity vs Time')
txt = ['{\color[rgb]{0, 0.4470, 0.7410} blue = LB, \color[rgb]{0.4940, 0.1840, 0.5560} purple = EDTA, \color[rgb]{0.4660, 0.6740, 0.1880} green = Mg2+, \color[rgb]{0.6350, 0.0780, 0.1840} red = *PBS,}'];
subtitle(txt)
xlabel('Time (min)')
ylabel('Intensity (A.U.)')
saveas(gcf, 'intensity.fig')
saveas(gcf, 'intensity.png')
close

figure
h1=histogram(-1./coeff1(:,2), 'FaceColor', [0, 0.4470, 0.7410]); hold on
h2=histogram(-1./coeff2(:,2), 'FaceColor', [0.4940, 0.1840, 0.5560]);
h3=histogram(-1./coeff3(:,2), 'FaceColor', [0.4660, 0.6740, 0.1880]); 
%h4=histogram(abs(coeff4(:,2))); 
h1.NumBins=12;
h2.NumBins=34;
h3.NumBins=12;
%h4.BinWidth=12;
title('Time Coefficient Distribution')
xlabel('\tau (min)')
ylabel('Counts')
legend({'LB', 'EDTA', 'Mg^{2+}'})
saveas(gcf, 'tauHist.fig')
saveas(gcf, 'tauHist.png')
close

figure
h1=histogram(thalf1(:,1), 'FaceColor', [0, 0.4470, 0.7410]); hold on
h2=histogram(thalf2(:,1), 'FaceColor', [0.4940, 0.1840, 0.5560]);
h3=histogram(thalf3(:,1), 'FaceColor', [0.4660, 0.6740, 0.1880]); 
%h4=histogram(abs(coeff4(:,2))); 
h1.NumBins=12;
h2.NumBins=34;
h3.NumBins=12;
%h4.BinWidth=12;
title('t_{1/2} Distribution')
xlabel('t_{1/2} (min^{-})')
legend({'LB', 'EDTA', 'Mg^{2+}'})
saveas(gcf, 'thalfHist.fig')
saveas(gcf, 'thalfHist.png')
close

%is there a difference in the overall dL of non-halved LB vs EDTA treated
%cells? (so far, it seems not)
figure
h1=histogram(yData1.("Overall dL"), 'FaceColor', [0, 0.4470, 0.7410]); hold on
h2=histogram(yData2.("Overall dL"), 'FaceColor', [0.4940, 0.1840, 0.5560]);
h3=histogram(yData3.("Overall dL"), 'FaceColor', [0.4660, 0.6740, 0.1880]);

h1.BinWidth = 0.5;
h2.BinWidth = 0.5;
h3.BinWidth = 0.5;
title('Overall Percent Change in Cell Wall Length')
xlabel('Percent Change (%)')
legend({'LB', 'EDTA', 'Mg^{2+}'})
saveas(gcf, 'overalldL.fig')
saveas(gcf, 'overalldL.png')

figure
h1=histogram(yData1.("cell length")(:, tidx+3), 'FaceColor', [0, 0.4470, 0.7410]); hold on
h2=histogram(yData2.("cell length")(:, tidx+3), 'FaceColor', [0.4940, 0.1840, 0.5560]);
h3=histogram(yData3.("cell length")(:, tidx+3), 'FaceColor', [0.4660, 0.6740, 0.1880]); 
h1.BinWidth = 0.5;
h2.BinWidth = 0.5;
h3.BinWidth = 0.5;
title('Cell Length Distribution')
xlabel('Cell Length \mum')
legend({'LB', 'EDTA', 'Mg^{2+}'})
saveas(gcf, 'lengthHist.fig')
saveas(gcf, 'lengthHist.png')
%% save the data
cd(maindir)
save('diffusionAnalysis.mat')

%close all
%% Functions
% function [y]=sigmoidal(b,x)
%     xData=x.*b(1);
%     c=b(1)*b(2);
%     denom=1+exp(xData-c);
%     y=1./denom;
% end