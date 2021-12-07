%Author: Zarina Akbary
%Date: 05/04/2021
%Purpose: To combine and analyze dm data from mNeonGreen/mCherry diffusion
%experiments. 

clear, close all

% colorcode={'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};
% colorcode2={'#A1D5F7', '#FFBB9E', '#FFE3A3', '#EAA6F7', '#CBF09C', '#C7EEFF', '#FFB3C1'};

colorcode={[0 0.45 0.74], [0.85 0.33 0.1], [0.49 0.18 0.56], [0.47 0.67 0.19]};
colorcode2={[0.63 0.84 0.97], [1 0.73 0.62], [0.92 0.65 0.97], [0.8 0.94 0.61]};

%inputs
check1=1; %view norm. fluor vs time for LB experiments 
check2=1; 
check3=1;

%% determining whether the cell wall is rate-limiting
%for the following experiments, the goal is to determine whether the cell
%wall is the rate-limiting factor in mNeonGreen diffusion. The frame rate
%is 1 min for each experiment.

dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/12062021_analysis';
cd([dirsave '/MatFiles'])

datadir1=dir(['10232021_Exp1' '*dm.mat']); %LB, rep1
datadir2=dir(['10262021_Exp1' '*dm.mat']); %LB, rep2

datadir3=dir(['10232021_Exp2' '*dm.mat']); %PBS, 1 hour, rep 1
datadir4=dir(['10262021_Exp2' '*dm.mat']); %PBS, 1 hour, rep 2

datadir5=dir(['11302021_Exp1' '*dm.mat']); %PBS, 20 min
datadir6=dir(['12052021_Exp1' '*dm.mat']); %PBS, 2 min


[intensity1, adjintensity1, normintensity1, lCell1, time1, t1]=dataInput(datadir1, 6);
[intensity2, adjintensity2, normintensity2, lCell2, time2, t2]=dataInput(datadir2, 6);

[intensity3, adjintensity3, normintensity3, lCell3, time3, t3]=dataInput(datadir3, 3);
[intensity4, adjintensity4, normintensity4, lCell4, time4, t4]=dataInput(datadir4, 5);

[intensity5, adjintensity5, normintensity5, lCell5, time5, t5]=dataInput(datadir5, 5);
[intensity6, adjintensity6, normintensity6, lCell6, time6, t6]=dataInput(datadir6, 11);

%combine the data of the 'LB' replicates
intensity_LB=[intensity1(:, 1:98); intensity2];
adjintensity_LB=[adjintensity1(:, 1:98); adjintensity2];
normintensity_LB=[normintensity1(:, 1:93); normintensity2];
lCell_LB=[lCell1(:, 1:98); lCell2];

%combine the data of the 'PBS, 1 hour' replicates
intensity_PBS60=[intensity3; intensity4(:, 1:51)];
adjintensity_PBS60=[adjintensity3; adjintensity4(:, 1:51)];
normintensity_PBS60=[normintensity3; normintensity4(:, 1:49)];
lCell_PBS60=[lCell3; lCell4(:, 1:51)];

%% sanity check: normalization
if check1==1
    cd(dirsave)
    figure, hold on
    for i=1:height(normintensity_LB)
        plot(t2, normintensity_LB(i,:), 'Color', colorcode{1})
    end
    for i=1:height(normintensity_PBS60)
        plot(t3, normintensity_PBS60(i,:), 'Color', colorcode{2})
    end
    for i=1:height(normintensity5)
        plot(t5, normintensity5(i,:), 'Color', colorcode{3})
    end
    for i=1:height(normintensity6)
        plot(t6, normintensity6(i,:), 'Color', colorcode{4})
    end
    ylabel('Normalized Fluorescence (A.U.)')
    xlabel('Time (min)')
    title('Normalized Fluorescence vs Time')
    txt = ['{\color[rgb]{0 0.45 0.74}0 min \color[rgb]{0.85 0.33 0.1} 1 hour \color[rgb]{0.49 0.18 0.56} 20 min \color[rgb]{0.47 0.67 0.19} 2 min}'];
    subtitle(txt)
    saveas(gcf, 'normintensity.png')
    saveas(gcf, 'normintensity.fig')
end

%% fit the plots to an exponential model to find tau
graph=0;

[fit_LB, tau_LB, yhat_LB]=expModel(t2, normintensity_LB, 'LB',  graph);
[fit3_PBS60, tau_PBS60, yhat_PBS60]=expModel(t3, normintensity_PBS60, 'PBS, 1 hour',  graph);
[fit5, tau5, yhat5]=expModel(t5, normintensity5, 'PBS, 20 min',  graph);
[fit6, tau6, yhat6]=expModel(t6, normintensity6, 'PBS, 2 min',  graph);

%plot tau as a function of frame rate and determine the slope when 1 min is
%included vs when 1 min data is not included

% compare the time constants
% find the mean and standard deviation of the time constants
tau_means = [mean(tau_LB, 'omitnan'), mean(tau6, 'omitnan'), mean(tau5, 'omitnan'), mean(tau_PBS60, 'omitnan')];
tau_std = [std(tau_LB, 0, 'omitnan'), std(tau6, 0, 'omitnan'), std(tau5, 0, 'omitnan'), std(tau_PBS60, 0, 'omitnan')];

ts_LB=repelem(0, length(tau_LB));
ts_PBS60=repelem(60, length(tau_PBS60));
ts5=repelem(20, length(tau5));
ts6=repelem(2, length(tau6));

linearCoef1 = polyfit([ts_LB, ts6, ts5, ts_PBS60],[tau_LB, tau6, tau5, tau_PBS60],1);
linearFit1= polyval(linearCoef1,[ts_LB, ts6, ts5, ts_PBS60]);
    
if check2==1
    figure, hold on
    scatter(ts_LB, tau_LB, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
    scatter(ts6, tau6, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
    scatter(ts5, tau5, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
    scatter(ts_PBS60, tau_PBS60, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')   
    scatter([0, 2, 20, 60], tau_means, 'MarkerFaceColor', 'black')
    errorbar(0, tau_means(1), tau_std(1), 'Color', 'black')
    errorbar(2, tau_means(2), tau_std(2), 'Color', 'black')
    errorbar(20, tau_means(3), tau_std(3), 'Color', 'black')
    errorbar(60, tau_means(4), tau_std(4), 'Color', 'black')
    plot([ts_LB, ts6, ts5, ts_PBS60], linearFit1, '--b')
    xlim([-2, 61])
    ylabel('\tau (min^{-1})')
    %xticks([0 2 0.25 0.5 1])
    %xticklabels({'', 'LB, 1.2 s', 'LB, 15 s', 'LB, 30 s', 'LB, 1 min',})
    title('Tau vs Time in PBS')
    saveas(gcf, 'tau.png')
    saveas(gcf, 'tau.fig')
    %txt = ['{\color {blue} slope calculated with 1 min frame rate \color {red} slope calculated without 1 min frame rate}'];
    %subtitle(txt)
    %pause, close all
end


%% Correct LB and PBS traces for photobleaching
%correct for photobleaching
parameter=16.9441;
[Cnew_LB, dCP_LB, dCB_LB, dCT_LB, unb_frac_LB, Cbl_exp_LB, gamma_LB] = photoCorrect(normintensity_LB, t2, parameter);
[Cnew_PBS60, dCP_PBS60, dCB_PBS60, dCT_PBS60, unb_frac_PBS60, Cbl_exp_PBS60, gamma_PBS60] = photoCorrect(normintensity_PBS60, t3, parameter);
[Cnew5, dCP5, dCB5, dCT5, unb_frac5, Cbl_exp5, gamma5] = photoCorrect(normintensity5, t5, parameter);
[Cnew6, dCP6, dCB6, dCT6, unb_frac6, Cbl_exp6, gamma6] = photoCorrect(normintensity6, t6, parameter);

if check3==1
    figure, hold on
    for i=1:height(normintensity_LB)
        plot(t2, Cnew_LB(i,:), 'Color', colorcode{1})
    end
    for i=1:height(normintensity_PBS60)
        plot(t3, Cnew_PBS60(i,:), 'Color', colorcode{2})
    end
    for i=1:height(normintensity5)
        plot(t5, Cnew5(i,:), 'Color', colorcode{3})
    end
    for i=1:height(normintensity6)
        plot(t6, Cnew6(i,:), 'Color', colorcode{4})
    end
    ylabel('Normalized Fluorescence (A.U.)')
    xlabel('Time (min)')
    title('Normalized Fluorescence vs Time')
    txt = ['{\color[rgb]{0 0.45 0.74} Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour  \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
    subtitle(txt)
    saveas(gcf, ['correction_' num2str(parameter) '.png'])
    saveas(gcf, ['correction_' num2str(parameter) '.fig'])
end

% Cavg_LB = mean(Cnew_LB, 1, 'omitnan');
% Cavg_PBS60 = mean(Cnew_PBS60, 1, 'omitnan');
% Cavg8 = mean(Cnew8, 1, 'omitnan');
% Cavg9 = mean(Cnew9, 1, 'omitnan');

% if check5==1
%     [Cnew_LB, dCP_LB, dCB_LB, dCT_LB, unb_frac_LB, Cbl_exp_LB, gamma_LB] = photoCorrect(normintensity_LB, t2, 1);
%     [Cnew_PBS60, dCP_PBS60, dCB_PBS60, dCT_PBS60, unb_frac_PBS60, Cbl_exp_PBS60, gamma_PBS60] = photoCorrect(normintensity_PBS60, t6, 1);
%     [Cnew8, dCP8, dCB8, dCT8, unb_frac8, Cbl_exp8, gamma8] = photoCorrect(normintensity8, t8, 1);
%     [Cnew9, dCP9, dCB9, dCT9, unb_frac9, Cbl_exp9, gamma9] = photoCorrect(normintensity9, t9, 1);
% 
%     figure(1), hold on
%     plot(t2, Cavg_LB, 'Color', colorcode{1})
%     plot(t6, Cavg_PBS60, 'Color', colorcode{2})
%     plot(t8, Cavg8, 'Color', colorcode{3})
%     plot(t9, Cavg9, 'Color', colorcode{4})
%     title('Normalized Fluorescence vs Time')
%     txt = ['{\color[rgb]{0 0.45 0.74}Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
%     subtitle(txt)
%     xlabel('Time (min)')
%     ylabel('Normalized Fluorescence (A.U.)')
%     
%     [Cnew_LB, dCP_LB, dCB_LB, dCT_LB, unb_frac_LB, Cbl_exp_LB, gamma_LB] = photoCorrect(normintensity_LB, t2, 2);
%     [Cnew_PBS60, dCP_PBS60, dCB_PBS60, dCT_PBS60, unb_frac_PBS60, Cbl_exp_PBS60, gamma_PBS60] = photoCorrect(normintensity_PBS60, t6, 2);
%     [Cnew8, dCP8, dCB8, dCT8, unb_frac8, Cbl_exp8, gamma8] = photoCorrect(normintensity8, t8, 2);
%     [Cnew9, dCP9, dCB9, dCT9, unb_frac9, Cbl_exp9, gamma9] = photoCorrect(normintensity9, t9, 2);
%     
%     figure(2), hold on
%     plot(t2, Cavg_LB, 'Color', colorcode{1})
%     plot(t6, Cavg_PBS60, 'Color', colorcode{2})
%     plot(t8, Cavg8, 'Color', colorcode{3})
%     plot(t9, Cavg9, 'Color', colorcode{4})
%     title('Normalized Fluorescence vs Time, parameter 2')
%     txt = ['{\color[rgb]{0 0.45 0.74}Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
%     subtitle(txt)
%     xlabel('Time (min)')
%     ylabel('Normalized Fluorescence (A.U.)')
%         
%     [Cnew_LB, dCP_LB, dCB_LB, dCT_LB, unb_frac_LB, Cbl_exp_LB, gamma_LB] = photoCorrect(normintensity_LB, t2, 3);
%     [Cnew_PBS60, dCP_PBS60, dCB_PBS60, dCT_PBS60, unb_frac_PBS60, Cbl_exp_PBS60, gamma_PBS60] = photoCorrect(normintensity_PBS60, t6, 3);
%     [Cnew8, dCP8, dCB8, dCT8, unb_frac8, Cbl_exp8, gamma8] = photoCorrect(normintensity8, t8, 3);
%     [Cnew9, dCP9, dCB9, dCT9, unb_frac9, Cbl_exp9, gamma9] = photoCorrect(normintensity9, t9, 3);
%     
%     figure(3), hold on
%     plot(t2, Cavg_LB, 'Color', colorcode{1})
%     plot(t6, Cavg_PBS60, 'Color', colorcode{2})
%     plot(t8, Cavg8, 'Color', colorcode{3})
%     plot(t9, Cavg9, 'Color', colorcode{4})
%     title('Normalized Fluorescence vs Time, parameter 3')
%     txt = ['{\color[rgb]{0 0.45 0.74}Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
%     subtitle(txt)
%     xlabel('Time (min)')
%     ylabel('Normalized Fluorescence (A.U.)')
% end
%% Functions
function [intensity, adjintensity, normintensity, lCell, tme, t]=dataInput(datadir, imstart)
        
        %pre-allocate variables
        intensity=[];
        lCell=[];
    
        for i=1:length(datadir)
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'lcell', 'time')

            for n=1:height(icell_intensity)
                if ~isnan(icell_intensity(n, imstart))&~isnan(icell_intensity(n, end)) %make sure there are fluor readings during the initial frame and the final frame, otherwise the adjust. and norm. will look off
                    intensity=[intensity; icell_intensity(n, :)];
                    lCell=[lCell; lcell(n, :)];
                end
            end

            if i==1
                tme=time;
                t=time(imstart:end)-time(imstart); %new time vector
            end
        end

            %adjust for background
            adjintensity = intensity-intensity(:, end);

            %normalize to initial frame
            normintensity=adjintensity(:, imstart:end)./adjintensity(:,imstart);
            normintensity(normintensity<0)=0;
            
end

function [fit, tau, yhat]=expModel(time, normintensity, text,  graph)
  
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

        pause, close all
    end
end

function [Cnew, dCP, dCB, dCT, unb_frac, Cbl_exp, gamma] = photoCorrect(normintensity, t, parameter)
    
    %calculate dt (there's variation of dt during initial values, so it's
    %easier to use end values)
    dt=abs(t(end-1)-t(end));
    %beta = 1/((16.9441)*dt + 0.1823); %1/tau
    beta = 1/((parameter)*dt + 0.1823); %1/tau
    alpha=1/beta/dt; %tau*dt
    
    %assume that the initial 'measured' fluorescence values and corrected
    %fluor. values will be equal. I prefer to pre-allocate with nan in case 
    %some values are missing in the raw data
    Cnew=nan(size(normintensity));%Corrected concentration of fluorophores
    Cnew(:, 1)=normintensity(:, 1);

    %pre-allocate variables
    dCB=nan(height(normintensity), length(t)-1);
    dCT=nan(height(normintensity), length(t)-1);
    dCP=nan(height(normintensity), length(t)-1); %this is the dCP, or loss attributable to permeability

    unb_frac=nan(size(normintensity)); %fraction of unbleached fluor. 
    unb_frac(:, 1)=1;%all fluorophores are unbleached at the initial time point

    Cbl_exp=nan(size(normintensity));%Calculated (from experiment and photobleaching constant) concentration of bleached flurophores
    Cbl_exp(:, 1)=0;
    
    for n=1:height(normintensity)
        
        %if parameter==1            
            for i=1:length(t)-1

                dCB(n,i) = normintensity(n,i)/(alpha); %this is the amount of photobleaching that occured in our measured value

%                 if dCB(n,i)<10e-3
%                     break
%                 end
                
                dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
                
                dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability

                dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached

                Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);

                Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
                unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores

            end       
%         elseif parameter==2
%               for i=1:length(t)-1
% 
%                 dCB(n,i) = ((normintensity(n,i)+normintensity(n,i+1))/2)/(alpha); %this is the amount of photobleaching that occured in our measured value
% 
% %                 if dCB(n,i)<10e-3
% %                     break
% %                 end
%                 
%                 dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
%                 
%                 dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability
% 
%                 dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached
% 
%                 Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);
% 
%                 Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
%                 unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores
% 
%               end  
%         elseif parameter==3
%                for i=1:length(t)-1
% 
%                 dCB(n,i) = normintensity(n,i+1)/(alpha); %this is the amount of photobleaching that occured in our measured value
% 
% %                 if dCB(n,i)<10e-3
% %                     break
% %                 end
%                 
%                 dCT(n,i) = normintensity(n, i+1) - normintensity(n, i); %this is the total fluor. loss for the measured value
%                 
%                 dCP(n, i) = dCT(n,i) + dCB(n,i); %this is the amount of loss attributable to permeability
% 
%                 dCP(n,i)=dCP(n,i)/unb_frac(n,i);%Correcting for the fact that a fraction of fluorophores are unbleached
% 
%                 Cnew(n,i+1)=Cnew(n,i)+dCP(n,i);
% 
%                 Cbl_exp(n,i+1)=Cbl_exp(n,i)+dCB(n,i)+dCP(n,i)*(1-unb_frac(n,i));%Accounting fo the change in concentration fo bleached fluorophores
%                 unb_frac(n,i+1)=(normintensity(n,i+1))/(normintensity(n,i+1)+Cbl_exp(i+1));%Calculate the new fraction of unbleached fluorophores
% 
%               end       
%         end
    end
    
    gamma = nan(height(normintensity), 1);
    
    for n=1:height(normintensity)
        p = dCP(n,end)/dt;
        gamma(n,1) = p/Cnew(n, end-1);
    end
end