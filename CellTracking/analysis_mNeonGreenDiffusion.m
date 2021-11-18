%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

%inputs
lengthTraces=1;
intensityTraces=1;
adjTraces=1;
normTraces=1;
modelFit=1;
modelCheck=0;
modelInterp=0;
bleachCorrect=1;

colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
colorcode2={'#A1D5F7', '#FFBB9E', '#FFE3A3', '#EAA6F7', '#CBF09C', '#C7EEFF', '#FFB3C1'};

% colorcode={[0 0.45 0.74], [0.85 0.33 0.1], [0.49 0.18 0.56], [0.47 0.67 0.19]};
% colorcode2={[0.63 0.84 0.97], [1 0.73 0.62], [0.92 0.65 0.97], [0.8 0.94 0.61]};
%% load the data
dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/11152021_analysis';
cd([dirsave '/MatFiles'])

datadir1=dir(['10302021_Exp2' '*dm.mat']); %LB, frame rate=15 s
datadir2=dir(['11152021_Exp1' '*dm.mat']); %sodium azide, frame rate=15 s

datadir3=dir(['10232021_Exp1' '*dm.mat']); %LB, rep 1, frame rate=1 min
datadir4=dir(['10262021_Exp1' '*dm.mat']); %LB, rep 2, frame rate=1 min

datadir5=dir(['10232021_Exp2' '*dm.mat']); %PBS, rep 1, frame rate=1 min
datadir6=dir(['10262021_Exp2' '*dm.mat']); %PBS, rep 2, frame rate=1 min

datadir7=dir(['11152021_Exp2' '*dm.mat']); %lysozyme, frame rate=1 min

[intensity1, adjintensity1, normintensity1, positions1, lCell1, time1, t1, lidx1]=dataInput(datadir1, 0);
[intensity2, adjintensity2, normintensity2, positions2, lCell2, time2, t2, lidx2]=dataInput(datadir2, 17);

[intensity3, adjintensity3, normintensity3, positions3, lCell3, time3, t3, lidx3]=dataInput(datadir3, 0);
[intensity4, adjintensity4, normintensity4, positions4, lCell4, time4, t4, lidx4]=dataInput(datadir4, 0);

[intensity5, adjintensity5, normintensity5, positions5, lCell5, time5, t5, lidx5]=dataInput(datadir5, 0);
[intensity6, adjintensity6, normintensity6, positions6, lCell6, time6, t6, lidx6]=dataInput(datadir6, 0);

[intensity7, adjintensity7, normintensity7, positions7, lCell7, time7, t7, lidx7]=dataInput(datadir7, 20);

%combine LB, 1 min replicates
intensityLB = [intensity3(:, 1:98);intensity4];
adjintensityLB = [adjintensity3(:, 1:98);adjintensity4];
normintensityLB = [normintensity3(:, 1:93);normintensity4];
lCellLB = [lCell3(:, 1:98); lCell4];
timeLB = time4;
tLB=t4;

%combine PBS replicates
intensityPBS = [intensity5;intensity6(:, 1:51)];
adjintensityPBS = [adjintensity5;adjintensity6(:, 1:51)];
normintensityPBS = [normintensity5;normintensity6(:, 1:48)];
lCellPBS = [lCell5;lCell5(:, 1:51)];
timePBS = time5;
tPBS=t5;

%% sanity check: length traces
if lengthTraces==1
    
    %make sure lysed cells lyse and sodium azide treated cells don't
    cd([dirsave '/LengthTraces'])
    [l_avg1, l_std1]=lengthPlot(time1, lCell1,  'LB perfusion, frame rate = 15 s', 1, 1, '10302021_Exp2'); 
    [l_avg2, l_std2]=lengthPlot(time2, lCell2, 'sodium azide treatment, frame rate = 15 s', 1, 1, '11152021_Exp1');
    [l_avgLB, l_stdLB]=lengthPlot(timeLB, lCellLB,  'LB perfusion, frame rate = 1 min', 1, 1, 'LB'); 
    [l_avgPBS, l_stdPBS]=lengthPlot(timePBS, lCellPBS, 'PBS incubation, frame rate = 1 min', 1, 1, 'PBS');
    [l_avg7, l_std7]=lengthPlot(time7, lCell7, 'lysozyme treatment, frame rate = 1 min', 1, 1, '11152021_Exp2');
    
end

%% sanity check: fluorescence traces
if intensityTraces==1
    
    %check to see that the initial intensity values are comparable 
    cd([dirsave '/IntensityTraces'])
    [intensity_avg1, intensity_std1]=intensityPlot(time1, intensity1, 'LB perfusion, frame rate = 15 s', 1, 1, '10302021_Exp2');
    [intensity_avg2, intensity_std2]=intensityPlot(time2, intensity2, 'sodium azide treatment, frame rate = 15 s', 1, 1, '11152021_Exp1');
    [intensity_avgLB, intensity_stdLB]=intensityPlot(timeLB, intensityLB,  'LB perfusion, frame rate = 1 min', 1, 1, 'LB'); 
    [intensity_avgPBS, intensity_stdPBS]=intensityPlot(timePBS, intensityPBS, 'PBS incubation, frame rate = 1 min', 1, 1, 'PBS');
    [intensity_avg7, intensity_std7]=intensityPlot(time7, intensity7, 'lysozyme treatment, frame rate = 1 min', 1, 1, '11152021_Exp2');
    
%     figure, hold on
%     for n=1:height(intensity1)
%         plot(time1, intensity1(n,:), 'Color', colorcode{1})
%     end
%     for n=1:height(intensity2)
%         plot(time2, intensity2(n,:), 'Color', colorcode{2})
%     end
%     ylim([0, Inf])
%     text(25,3000,'LB','Color', colorcode{1},'FontSize',10)
%     text(25,2500,'Sodium Azide','Color', colorcode{2},'FontSize',10)
%     saveas(gcf, 'intensityTraces.fig')
%     saveas(gcf, 'intensityTraces.png')
end

%% sanity check: adjusted fluorescence traces
if adjTraces==1
    
    %do the adjusted traces both go to zero?
    cd([dirsave '/AdjTraces'])
    [adjintensity_avg1, adjintensity_std1]=adjintensityPlot(time1, adjintensity1, 'LB perfusion, frame rate = 15 s', 1, 1, '10302021_Exp2');
    [adjintensity_avg2, adjintensity_std2]=adjintensityPlot(time2, adjintensity2, 'sodium azide treatment, frame rate = 15 s', 1, 1, '11152021_Exp1');
    [adjintensity_avgLB, adjintensity_stdLB]=adjintensityPlot(timeLB, adjintensityLB, 'LB perfusion, frame rate = 1 min', 1, 1, 'LB');
    [adjintensity_avgPBS, adjintensity_stdPBS]=adjintensityPlot(timePBS, adjintensityPBS, 'PBS incubation, frame rate = 1 min', 1, 1, 'PBS');
    [adjintensity_avg7, adjintensity_std7]=adjintensityPlot(time7, adjintensity7, 'lysozyme treatment, frame rate = 1 min', 1, 1, '11152021_Exp2');

end

%% normalized fluorescence traces
if normTraces==1
    
    %do the traces start at 1 and end at zero?
    cd([dirsave '/NormTraces'])
    [normintensity_avg1, normintensity_std1]=normintensityPlot(t1, normintensity1, 'LB perfusion, frame rate = 15 s', 1, 1, '10302021_Exp2');
    [normintensity_avg2, normintensity_std2]=normintensityPlot(t2, normintensity2, 'sodium azide treatment, frame rate = 15 s', 1, 1, '11152021_Exp1');
    [normintensity_avgLB, normintensity_stdLB]=normintensityPlot(tLB, normintensityLB, 'LB perfusion, frame rate = 1 min', 1, 1, 'LB');
    [normintensity_avgPBS, normintensity_stdPBS]=normintensityPlot(tPBS, normintensityPBS, 'PBS incubation, frame rate = 1 min', 1, 1, 'PBS');
    [normintensity_avg7, normintensity_std7]=normintensityPlot(t7, normintensity7, 'lysozyme treatment, frame rate = 1 min', 1, 1, '11152021_Exp2');

end

%% fit the data to an exponential model
if modelFit==1

    %remember, none of the tau values should be zero
    cd([dirsave '/ModelFit'])
    [fit1, tau1, yhat1]=expModel(t1, normintensity1, 'LB perfusion, frame rate = 15 s', 1, 1, '10302021_Exp2');
    [fit2, tau2, yhat2]=expModel(t2, normintensity2,  'sodium azide treatment, frame rate = 15 s', 1, 1, '11152021_Exp1');
    [fitLB, tauLB, yhatLB]=expModel(tLB, normintensityLB, 'LB perfusion, frame rate = 1 min', 1, 1, 'LB');
    [fitPBS, tauPBS, yhatPBS]=expModel(tPBS, normintensityPBS,  'PBS incubation, frame rate = 1 min', 1, 1, 'PBS');
    [fit7, tau7, yhat7]=expModel(t7, normintensity7,  'lysozyme treatment, frame rate = 1 min', 1, 1, '11152021_Exp2');

    
    % compare the time constants
    % find the mean and standard deviation of the time constants
    tau_means = [mean(tau1, 'omitnan'), mean(tau2, 'omitnan'), mean(tauLB, 'omitnan'), mean(tauPBS, 'omitnan'), mean(tau7, 'omitnan')];
    tau_std = [std(tau1, 0, 'omitnan'), std(tau2, 0, 'omitnan'), std(tauLB, 0, 'omitnan'), std(tauPBS, 0, 'omitnan'), std(tau7, 0, 'omitnan')];

    figure, hold on
    scatter([1:5], tau_means, 'MarkerFaceColor', 'black')
    errorbar(1, tau_means(1), tau_std(1))
    errorbar(2, tau_means(2), tau_std(2))
    errorbar(3, tau_means(3), tau_std(3))
    errorbar(4, tau_means(4), tau_std(4))
    errorbar(5, tau_means(5), tau_std(5))
    xlim([0, 5.5])
    ylabel('\tau (min^{-1})')
    xticks([0 1 2 3 4 5])
    xticklabels({'', 'LB, 15 s', 'sodium azide', 'LB, 1 min', 'PBS', 'lysozyme'})
    saveas(gcf, 'meanTau.fig')
    saveas(gcf, 'meanTau.png')
    pause, close all
    
    %Perform an ANOVA to determine statistical significance
    anovaTau = [tau1(1:42)', tau2', tauLB(1:42)', tauPBS(1:42)', tau7(1:42)'];
    [p, tbl, stats]=anova1(anovaTau, categorical({'LB, 15 s', 'sodium azide', 'LB, 1 min', 'PBS', 'lysozyme'}));
    compTau=multcompare(stats);
    
    %what is tau as a function of initial fluor.?
    %initFluor = [intensity1(:, lidx1-1)', intensity2(:, 17)', intensity3(:, lidx3-1)', intensity4(:, lidx4-1)', intensity5(:, lidx5-1)', intensity6(:, lidx6-1)', intensity7(:, 20)'];
    figure, hold on
    scatter(intensity1(:, lidx1-1), tau1, 'MarkerFaceColor', colorcode{1}, 'MarkerEdgeColor', colorcode{1})
    scatter(intensity2(:, 17), tau2, 'MarkerFaceColor', colorcode{2}, 'MarkerEdgeColor', colorcode{2})
    scatter([intensity3(:, lidx3-1); intensity4(:, lidx4-1)], tauLB, 'MarkerFaceColor', colorcode{3}, 'MarkerEdgeColor', colorcode{3})
    %scatter(intensity4(:, lidx4-1), tauLB(49:end), 'MarkerFaceColor', colorcode{3}, 'MarkerEdgeColor', colorcode{3})
    scatter([intensity5(:, lidx5-1); intensity6(:, lidx6-1)], tauPBS, 'MarkerFaceColor', colorcode{4}, 'MarkerEdgeColor', colorcode{4})
    %scatter(intensity6(:, lidx6-1), tauPBS(68:end), 'MarkerFaceColor', colorcode{4}, 'MarkerEdgeColor', colorcode{4})
    scatter(intensity7(:, 20), tau7, 'MarkerFaceColor', colorcode{5}, 'MarkerEdgeColor', colorcode{5})
    xlabel('Initial Fluorescence (A.U.)')
    ylabel('\tau (min^{-1})')
    ylim([0 35])
    legend('LB, 15 s', 'sodium azide', 'LB, 1 min', 'PBS', 'lysozyme')
    
end
%% determine the quality of the fit
if modelCheck==1
    
   cd([dirsave '/ModelCheck'])
   [residuals1, est1]=residualPlot(normintensity1, yhat1, t1, fit1, 'LB perfusion, frame rate = 15 s', 1, 1, '10302021_Exp2');
   [residuals2, est2]=residualPlot(normintensity2, yhat2, t2, fit2, 'sodium azide treatment, frame rate = 15 s', 1, 1, '11152021_Exp1');

end

%% interpolate t_{1/2}
if modelInterp==1
    
    cd([dirsave '/ModelInterp'])
    [xq1, vq1, thalf1]=interpPlot(normintensity1, t1, 'LB perfusion, frame rate = 15 s', 1, 1, '10302021_Exp2');
    [xq2, vq2, thalf2]=interpPlot(normintensity2, t2, 'sodium azide treatment, frame rate = 15 s', 1, 1, '11152021_Exp1');

end

%% Photobleach Correction Model
if bleachCorrect==1
    
    cd(dirsave)
    [Cnew1, dCP1, unb_frac1, Cbl_exp1, gamma1] = photoCorrect(normintensity1, t1, tau1, 3, 1, 1, '10302021_Exp2');
    [Cnew2, dCP2, unb_frac2, Cbl_exp2, gamma2] = photoCorrect(normintensity2, t2, tau2, 3, 1, 1, '11152021_Exp1');
    [CnewLB, dCPLB, unb_fracLB, Cbl_expLB, gammaLB] = photoCorrect(normintensityLB, tLB, tauLB, 3, 1, 1, 'LB');
    [CnewPBS, dCPPBS, unb_fracPBS, Cbl_expPBS, gammaPBS] = photoCorrect(normintensityPBS, tPBS, tauPBS, 3, 1, 1, 'PBS');
    [Cnew7, dCP7, unb_frac7, Cbl_exp7, gamma7] = photoCorrect(normintensity7, t7, tau7, 3, 1, 1, '11152021_Exp2');
    
    figure, hold on
    %ciplot(mean(Cnew1,1, 'omitnan')-std(Cnew1,0, 1, 'omitnan'), mean(Cnew1,1, 'omitnan')+std(Cnew1,0, 1, 'omitnan'), t1, colorcode2{1})
    plot(t1, mean(Cnew1,1, 'omitnan'), 'Color', colorcode{1})
    %ciplot(mean(Cnew2,1, 'omitnan')-std(Cnew2,0, 1, 'omitnan'), mean(Cnew2,1, 'omitnan')+std(Cnew2,0, 1, 'omitnan'), t2, colorcode2{2})
    plot(t2, mean(Cnew2,1, 'omitnan'), 'Color', colorcode{2})
    %ciplot(mean(Cnew3,1, 'omitnan')-std(Cnew3,0, 1, 'omitnan'), mean(Cnew3,1, 'omitnan')+std(Cnew3,0, 1, 'omitnan'), t3, colorcode2{3})
    plot(tLB, mean(CnewLB,1, 'omitnan'), 'Color', colorcode{3})
    %ciplot(mean(Cnew4,1, 'omitnan')-std(Cnew4,0, 1, 'omitnan'), mean(Cnew4,1, 'omitnan')+std(Cnew4,0, 1, 'omitnan'), t4, colorcode2{4})
    plot(tPBS, mean(CnewPBS,1, 'omitnan'), 'Color', colorcode{4})
    %ciplot
    plot(t7, mean(Cnew7,1, 'omitnan'), 'Color', colorcode{5})
    ylim([0 Inf])
    xlabel('Time (minutes)')
    ylabel('Normalized Fluorescence (A.U.)')
    %legend('LB std, 15 s', 'LB avg, 15 s', 'SA std', 'SA avg', 'LB std, 1 min', 'LB avg, 1 min', 'lyso std', 'lyso avg')
    legend('LB, 15 s', 'SA', 'LB, 1 min', 'PBS', 'lysozyme')
    saveas(gcf, 'correctedTraces.fig')
    saveas(gcf, 'correctedTraces.png')
    
end

%% Functions
function [intensity, adjintensity, normintensity, positions, lCell, time, t, lidx]=dataInput(datadir, imstart)
    
    if imstart==0
        %to determine the initial lysis frame, go to the 647 image stack of
        %colony 1, which is usually to the right-most part of the chip
        cd(datadir(1).folder)
        load(datadir(1).name, 'channels')

        for t=5:9 %it doesn't start at 1 because the imaging sometimes starts during PBS perfusion, and that has 647 in it as well

            %to determine the lysis frame, go to the 647 image stack
            new_channel = strrep(channels,'mNeonGreen','647');
            cd(new_channel{1}); 
            fluo_directory=dir('*.tif');

            %binarize the image
            imagename=fluo_directory(t).name;
            im=imread(imagename);
            nim=imbinarize(im);

            if mean(mean(nim))<1 %frames without 647 will have a mean of 1
                lidx=t;
                break
            end
        end

    %pre-allocate variables
    intensity=[];
    lCell=[];
    positions=[0]; %each element is the cumulative sum of the cells in each movie
    
    for i=1:length(datadir)
        cd(datadir(i).folder)
        load(datadir(i).name, 'icell_intensity', 'lcell', 'time')
        count=0;
        
        for n=1:height(icell_intensity)
            if ~isnan(icell_intensity(n, lidx-1))&~isnan(icell_intensity(n, end)) %want to make sure that the pre-lysis and end time points have readings, otherwise the normalization and adjustment will be NaN
                intensity=[intensity; icell_intensity(n, :)];
                lCell=[lCell; lcell(n, :)];
                count=count+1;
            end
        end
        
        positions(i+1)=count+positions(i); %this is to distinguish which rows are from which colony
     
    end
    
        %I don't know why there is a spike in fluor. upon detergent
        %perfusion, but maybe it's messing with the downstream analysis, so
        %I'll just set that point to NaN
        adjintensity = intensity-intensity(:, end);

        %normalize to pre-lysis frame instead of initial frame. The
        %earliest lysis frame is lidx
        normintensity=adjintensity(:, lidx-1:end)./adjintensity(:,lidx-1);
        normintensity(normintensity<0)=0;
        
        %new time vector (0 = just prior to detergent perfusion)
        t=time(lidx-1:end)-time(lidx-1);
    
    else 
        
        %pre-allocate variables
        intensity=[];
        lCell=[];
        positions=[0];
    
        for i=1:length(datadir)
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'lcell', 'time')
            count=0;

            for n=1:height(icell_intensity)
                if ~isnan(icell_intensity(n, imstart))&~isnan(icell_intensity(n, end)) %make sure there are fluor readings during the initial frame and the final frame, otherwise the adjust. and norm. will look off
                    intensity=[intensity; icell_intensity(n, :)];
                    lCell=[lCell; lcell(n, :)];
                    count=count+1;
                end
            end

            positions(i+1)=count+positions(i); %this is to distinguish which rows are from which colony

        end

            %adjust for background
            adjintensity = intensity-intensity(:, end);

            %normalize to initial frame
            normintensity=adjintensity(:, imstart:end)./adjintensity(:,imstart);
            normintensity(normintensity<0)=0;

            %new time vector
            t=time(imstart:end)-time(imstart);
            lidx=NaN;

        end
        
end

function [l_avg, l_std]=lengthPlot(time, lCell, text, graph, save, saveAs)

    l_avg = mean(lCell, 1, 'omitnan');
    l_std = std(lCell, 0, 1, 'omitnan');
    
    if graph==1
        figure(1), hold on
        for i=1:height(lCell)
            plot(time, lCell(i,:), 'Color', '#0072BD')    
        end
        xlabel('Time (minutes)')
        ylabel('Length (\mum)')
        title(text)

        if save==1
            saveas(gcf, [saveAs '_lTraces.png'])
            saveas(gcf, [saveAs '_lTraces.fig'])
        end

        figure(2), hold on
        ciplot((l_avg-l_std),(l_avg+l_std),time,[0.75 0.75 1])
        plot(time,l_avg,'-r')
        ylim([0 Inf])
        xlabel('Time (minutes)')
        ylabel('Length Average (\mum)')
        title(text)

        if save==1
            saveas(gcf, [saveAs '_lavgTraces.png'])
            saveas(gcf, [saveAs '_lavgTraces.fig'])
        end

        pause, close all
    end
end

function [intensity_avg, intensity_std]=intensityPlot(time, intensity, text, graph, save, saveAs)
    
    intensity_avg=mean(intensity, 1, 'omitnan');
    intensity_std=std(intensity, 0, 1, 'omitnan');
    
    if graph==1
        figure(1), hold on
        for i=1:height(intensity)
            plot(time, intensity(i, :), 'Color', '#77AC30')
        end
        ylim([0 Inf])
        xlabel('Time (minutes)')
        ylabel('Fluorescence (A.U.)')
        title(text)
        if save==1
            saveas(gcf, [saveAs '_intensityTraces.png'])
            saveas(gcf, [saveAs '_intensityTraces.fig'])
        end
        
        figure(2), hold on
        ciplot((intensity_avg-intensity_std),(intensity_avg+intensity_std),time,[0.75 0.75 1])
        plot(time,intensity_avg,'-r')
        ylim([0 Inf])
        xlabel('Time (minutes)')
        ylabel('Average Fluorescence (A.U.)')
        title(text)
        
        if save==1
            saveas(gcf, [saveAs '_intensityAvg.png'])
            saveas(gcf, [saveAs '_intensityAvg.fig'])
        end

        pause, close all
    end
end

function [adjintensity_avg, adjintensity_std]=adjintensityPlot(time, adjintensity, text, graph, save, saveAs)
    
    adjintensity_avg=mean(adjintensity, 1, 'omitnan');
    adjintensity_std=std(adjintensity, 0, 1, 'omitnan');
    
    if graph==1
        figure, hold on
        for i=1:height(adjintensity)
            plot(time, adjintensity(i, :), 'Color', '#77AC30')
        end
        xlabel('Time (minutes)')
        ylabel('Adjusted Fluorescence (A.U.)')
        ylim([0 Inf])
        title(text)
        if save==1
            saveas(gcf, [saveAs '_adjintensityTraces.png'])
            saveas(gcf, [saveAs '_adjintensityTraces.fig'])
        end

        figure(2), hold on
        ciplot((adjintensity_avg-adjintensity_std),(adjintensity_avg+adjintensity_std),time,[0.75 0.75 1])
        plot(time,adjintensity_avg,'-r')
        ylim([0 Inf])
        xlabel('Time (minutes)')
        ylabel('Adjusted Fluorescence (A.U.)')
        title(text)

        if save==1
            saveas(gcf, [saveAs '_adjintensityAvg.png'])
            saveas(gcf, [saveAs '_adjintensityAvg.fig'])
        end

        pause, close all
    end
end

function [normintensity_avg, normintensity_std]=normintensityPlot(time, normintensity, text,  graph, save, saveAs)
    
    normintensity_avg=mean(normintensity, 1, 'omitnan');
    normintensity_std=std(normintensity, 0, 1, 'omitnan');
    
    if graph==1
        figure, hold on
       for i=1:height(normintensity)
            plot(time, normintensity(i, :), 'Color', '#77AC30')
        end
        ylim([0 Inf])
        xlabel('Time (minutes)')
        ylabel('Normalized Fluorescence (A.U.)')
        title(text)
        if save==1
            saveas(gcf, [saveAs '_normintensityTraces.png'])
            saveas(gcf, [saveAs '_normintensityTraces.fig'])
        end

        figure(2), hold on
        ciplot((normintensity_avg-normintensity_std),(normintensity_avg+normintensity_std),time,[0.75 0.75 1])
        plot(time,normintensity_avg,'-r')
        ylim([0 Inf])
        xlabel('Time (minutes)')
        ylabel('Normalized Fluorescence (A.U.)')
        title(text)

        if save==1
            saveas(gcf, [saveAs '_normintensityAvg.png'])
            saveas(gcf, [saveAs '_normintensityAvg.fig'])
        end

        pause, close all
    end
end

function [fit, tau, yhat]=expModel(time, normintensity, text,  graph, save, saveAs)
  
    modelfun=@(tau,x)exp(-x./tau);
    tau0=1;

    %pre-allocate variables
    fit=[];
    tau=[];
    yhat=[];
    
    for i=1:height(normintensity)
            idx=find(normintensity(i,:)<=1); %this will skip the detergent perfusion reading
            tau_temp=nlinfit(time(idx), normintensity(i,idx), modelfun, tau0);
            y_hat=modelfun(tau_temp, time);   
            
%             residuals=abs(normintensity(i,:)-y_hat);
%             est=residuals./normintensity(i,:);
%             est(est==Inf)=0;
            
            fit=[fit, i];
            tau=[tau, tau_temp];
            yhat=[yhat; y_hat]; 

    end

    if graph==1
        figure, hold on
        for i=1:length(fit)
            plot(time, normintensity(fit(i),:),'Color', '#77AC30')
            plot(time, yhat(i,:), '--k')
        end
        ylim([0 Inf])
        legend({'data', 'model'})
        title(text)
        xlabel('Time (min)')
        ylabel('Normalized Intensity (A.U.)')
        if save==1
            saveas(gcf, [saveAs '_expFit.png'])
            saveas(gcf, [saveAs '_expFit.fig'])
        end

        pause, close all
    end
end

function [residuals, est]=residualPlot(normintensity, yhat, time, fit, text, graph, save, saveAs)

    residuals=abs(normintensity(fit,:)-yhat);
    est=residuals./normintensity(fit, :);
    est(est==Inf)=0;

    if graph==1
        figure, hold on
        for i=1:height(est)
            plot(time, est(i,:), 'Color', '#7E2F8E')
        end
        xlabel('Time (minutes)')
        ylabel('$\frac{|y-\hat{y}|}{y}$','Interpreter','latex', 'FontSize', 20)
        title(text)


        if save==1
            saveas(gcf, [saveAs '_est.png'])
            saveas(gcf, [saveAs '_est.fig'])
        end

        pause, close all
    end
end

function [xq, vq, thalf]=interpPlot(normintensity,time, text, graph, save, saveAs)

    xq=[0:0.0167:time(end)];
    vq=nan(height(normintensity), length(xq));
    
    for i=1:height(normintensity)
        vq(i,:)=interp1(time,normintensity(i,:),xq);
    end
    
    if graph==1
        figure, hold on
        for i=1:height(normintensity)
            plot(xq,vq(i,:),':.', 'Color', '#D95319');
            plot(time,normintensity(i,:),'o', 'Color','#0072BD');
        end
        ylim([0 Inf])
        xlabel('Time (minutes)')
        ylabel('Normalized Fluoresence (A.U.)')
        legend({'interpolation', 'data'})
        title([text ', Linear Interpolation']);

        [~, tidx]=min(abs(vq-0.5), [], 2);
        thalf=xq(tidx);

        if save==1
            saveas(gcf, [saveAs '_interp.png'])
            saveas(gcf, [saveAs '_interp.fig'])
        end

        close all
    end 
end

function [Cnew, dCP, unb_frac, Cbl_exp, gamma] = photoCorrect(normintensity, t, tau, parameter, graph, save, saveAs)
    
    %calculate dt (there's variation of dt during initial values, so it's
    %easier to use end values)
    dt=abs(t(end-1)-t(end));
    beta = 1/((16.0806)*dt + 0.4234); %1/tau
    alpha=1/beta/dt; %tau*dt
    
    %assume that the initial 'measured' fluorescence values and corrected
    %fluor. values will be equal. I prefer to pre-allocate with nan in case 
    %some values are missing in the raw data
    Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
    Cnew(:, 1)=normintensity(:, 1);

    %this is the dCP, or loss attributable to permeability
    dCP=nan(height(normintensity), length(t)-1);

    unb_frac=nan(size(normintensity)); %fraction of unvbleached fluor. 
    unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

    Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
    Cbl_exp(:, 1)=0;
    
    for n=1:height(normintensity)
        
        if parameter==1
            for i=1:length(t)-1

                dCB = normintensity(n,i)/(alpha); %this is the amount of photobleaching that occured in our measured value
    %             if dCB<0 %no negative dCB for normintensity1
    %                 disp(['negative dCB for cell ' num2str(n) ', t = ' num2str(i)])
    %             end

                dCT = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
    %             if dCT>0 %many positive dCT for normintensity1
    %                 disp(['positive dCT for cell ' num2str(n) ', t = ' num2str(i)])
    %                 %dCT=-dCT;
    %                 %continue 
    %             end

                dCP(n, i) = dCT + dCB; %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores

            end
        elseif parameter==2
            
              for i=1:length(t)-1

                dCB = ((normintensity(n,i)+normintensity(n,i+1))/2)/(alpha); %this is the amount of photobleaching that occured in our measured value
    %             if dCB<0 %no negative dCB for normintensity1
    %                 disp(['negative dCB for cell ' num2str(n) ', t = ' num2str(i)])
    %             end

                dCT = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
    %             if dCT>0 %many positive dCT for normintensity1
    %                 disp(['positive dCT for cell ' num2str(n) ', t = ' num2str(i)])
    %                 %dCT=-dCT;
    %                 %continue 
    %             end

                dCP(n, i) = dCT + dCB; %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores

              end
              
        elseif parameter==3
            
              for i=1:length(t)-1

                dCB = normintensity(n,i+1)/(alpha); %this is the amount of photobleaching that occured in our measured value
    %             if dCB<0 %no negative dCB for normintensity1
    %                 disp(['negative dCB for cell ' num2str(n) ', t = ' num2str(i)])
    %             end

                dCT = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
    %             if dCT>0 %many positive dCT for normintensity1
    %                 disp(['positive dCT for cell ' num2str(n) ', t = ' num2str(i)])
    %                 %dCT=-dCT;
    %                 %continue 
    %             end

                dCP(n, i) = dCT + dCB; %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores

              end    
        end
        
    end
    
    gamma = nan(height(normintensity), 1);
    for n=1:height(normintensity)
        p = dCP(n,end)/dt;
        gamma(n,1) = p/Cnew(n, end-1);
    end
    
    if graph==1
       figure, hold on
       for i=1:height(Cnew)
            plot(t, Cnew(i, :), 'Color', '#77AC30')
        end
        ylim([0 Inf])
        xlabel('Time (minutes)')
        ylabel('Normalized Fluorescence (A.U.)')
        if save==1
            saveas(gcf, [saveAs '_normintensityTraces.png'])
            saveas(gcf, [saveAs '_normintensityTraces.fig'])
        end

        pause, close all
    end
end

function ciplot(lower,upper,x,colour,Trans)
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Raymond Reynolds 24/11/06

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if nargin<4
    colour='b';
end

if nargin<3
    x=1:length(lower);
end

if nargin<5
    Trans=1;
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

fill([x fliplr(x)],[upper fliplr(lower)],colour,'LineStyle','none','FaceAlpha',Trans)
end