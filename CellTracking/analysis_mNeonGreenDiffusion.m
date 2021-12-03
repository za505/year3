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
check1=0; %view norm. fluor vs time for LB experiments 
check2=0;
check3=1;
check4=0;
check5=0;

%% determining the photobleaching constant, beta
%in the following experiments, cells were perfused with LB, lysed with
%detergent, and then perfused with LB once more. Given that mNeonGreen is
%greater than the putative size threshold for permeability, I expected that
%the fluorophore does not leak out of the cell and that the decrease in
%fluorescence is attributable to photobleaching. However, even if the
%fluorophore was cell wall permeable, the decrease in fluorscence would
%still be a function of frame rate, and from that linear relationship, the
%rate of photobleaching can be deduced. 

dirsave='/Users/zarina/Documents/MATLAB/MatlabReady/mNeonGreenDiffusion_analysis/12032021_analysis';
cd([dirsave '/MatFiles'])

datadir1=dir(['10232021_Exp1' '*dm.mat']); %LB, frame rate=1 min, rep1
datadir2=dir(['10262021_Exp1' '*dm.mat']); %LB, frame rate=1 min, rep2

datadir3=dir(['10302021_Exp1' '*dm.mat']); %LB, frame rate=30 s
datadir4=dir(['10302021_Exp2' '*dm.mat']); %LB, frame rate=15 s

datadir5=dir(['11202021_Exp1' '*dm.mat']); %LB, frame rate=1.2 s

[intensity1, adjintensity1, normintensity1, lCell1, t1]=dataInput(datadir1, 5);
[intensity2, adjintensity2, normintensity2, lCell2, t2]=dataInput(datadir2, 5);

[intensity3, adjintensity3, normintensity3, lCell3, t3]=dataInput(datadir3, 5);
[intensity4, adjintensity4, normintensity4, lCell4, t4]=dataInput(datadir4, 5);

[intensity5, adjintensity5, normintensity5, lCell5, t5]=dataInput(datadir5, 6);

%% Check the length traces to make sure cells aren't lysed
if check1==1
    figure(1), hold on
    for i=1:height(lCell1)
        plot(t1, lCell1(i,:), 'Color', colorcode{1})
    end
    ylabel('Length (\mum)')
    xlabel('Time (min)')
    title('Length vs Time, LB frame rate = 1 min')

    figure(2), hold on
    for i=1:height(lCell2)
        plot(t2, lCell2(i,:), 'Color', colorcode{1})
    end
    ylabel('Length (\mum)')
    xlabel('Time (min)')
    title('Length vs Time, LB frame rate = 1 min')

    figure(3), hold on
    for i=1:height(lCell3)
        plot(t3, lCell3(i,:), 'Color', colorcode{1})
    end
    ylabel('Length (\mum)')
    xlabel('Time (min)')
    title('Length vs Time, LB frame rate = 30 s')

    figure(4), hold on
    for i=1:height(lCell4)
        plot(t4, lCell4(i,:), 'Color', colorcode{1})
    end
    ylabel('Length (\mum)')
    xlabel('Time (min)')
    title('Length vs Time, LB frame rate = 15 s')

    figure(5), hold on
    for i=1:height(lCell5)
        plot(t5, lCell5(i,:), 'Color', colorcode{1})
    end
    ylabel('Length (\mum)')
    xlabel('Time (min)')
    title('Length vs Time, LB frame rate = 1.2 s')
end    

%% combine the data of the 'LB, frame rate=1 min' replicates
intensity_LB=[intensity1; intensity2];
adjintensity_LB=[adjintensity1; adjintensity2];
normintensity_LB=[normintensity1; normintensity2];

%% sanity check: normalization
if check2==1
    figure(1), hold on
    for i=1:height(normintensity_LB)
        plot(t2, normintensity_LB(i,:), 'Color', colorcode{1})
    end
    figure(2), hold on
    for i=1:height(normintensity3)
        plot(t3, normintensity3(i,:), 'Color', colorcode{2})
    end
    figure(3), hold on
    for i=1:height(normintensity4)
        plot(t4, normintensity4(i,:), 'Color', colorcode{3})
    end
    figure(4), hold on
    for i=1:height(normintensity5)
        plot(t5, normintensity5(i,:), 'Color', colorcode{4})
    end
    ylabel('Normalized Fluorescence (A.U.)')
    xlabel('Time (min)')
    title('Normalized Fluorescence vs Time')
    txt = ['{\color[rgb]{0 0.45 0.74}1 min \color[rgb]{0.85 0.33 0.1} 30 s \color[rgb]{0.49 0.18 0.56} 15 s \color[rgb]{0.47 0.67 0.19} 1.2 s}'];
    subtitle(txt)
end

%Observations: there seems to be a spike in fluorescence for some traces,
%but I am not sure why. Is it just noise?

%fit the plots to an exponential model to find tau
graph=0;

[fit_LB, tau_LB, yhat_LB]=expModel(t2, normintensity_LB, 'LB, frame rate=1 min',  graph)
[fit3, tau3, yhat3]=expModel(t3, normintensity3, 'LB, frame rate=30 s',  graph);
[fit4, tau4, yhat4]=expModel(t4, normintensity4, 'LB, frame rate=15 s',  graph);
[fit5, tau5, yhat5]=expModel(t5, normintensity5, 'LB, frame rate=1.2 s',  graph);

%% plot tau as a function of frame rate and determine the slope
% compare the time constants
% find the mean and standard deviation of the time constants
tau_means = [mean(tau5, 'omitnan'), mean(tau4, 'omitnan'), mean(tau3, 'omitnan'), mean(tau_LB, 'omitnan')];
tau_std = [std(tau5, 0, 'omitnan'), std(tau4, 0, 'omitnan'), std(tau3, 0, 'omitnan'), std(tau_LB, 0, 'omitnan')];

linearCoef1 = polyfit([repelem(1.2/60, length(tau5)), repelem(0.25, length(tau4)), repelem(0.5, length(tau3)), ones(1,length(tau_LB))],[tau5, tau4, tau3, tau_LB],1);
linearFit1= polyval(linearCoef1,[0, repelem(1.2/60, length(tau5)), repelem(0.25, length(tau4)), repelem(0.5, length(tau3)), ones(1,length(tau_LB))]);

linearCoef2 = polyfit([repelem(1.2/60, length(tau5)), repelem(0.25, length(tau4)), repelem(0.5, length(tau3))],[tau5, tau4, tau3],1);
linearFit2= polyval(linearCoef2,[0, repelem(1.2/60, length(tau5)), repelem(0.25, length(tau4)), repelem(0.5, length(tau3)), ones(1,length(tau_LB))]);
    
if check3==1
    figure, hold on
    scatter(repelem(1.2/60, length(tau5)), tau5, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
    scatter(repelem(0.25, length(tau4)), tau4, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
    scatter(repelem(0.5, length(tau3)), tau3, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')
    scatter(ones(1,length(tau_LB)), tau_LB, 'MarkerFaceColor', '#A1D5F7', 'MarkerEdgeColor', '#A1D5F7')   
    scatter([1.2/60, 0.25, 0.5, 1], tau_means, 'MarkerFaceColor', 'black')
    errorbar(1.2/60, tau_means(1), tau_std(1), 'Color', 'black')
    errorbar(0.25, tau_means(2), tau_std(2), 'Color', 'black')
    errorbar(0.5, tau_means(3), tau_std(3), 'Color', 'black')
    errorbar(1, tau_means(4), tau_std(4), 'Color', 'black')
    plot([0, repelem(1.2/60, length(tau5)), repelem(0.25, length(tau4)), repelem(0.5, length(tau3)), ones(1,length(tau_LB))], linearFit1, '--b')
    plot([0, repelem(1.2/60, length(tau5)), repelem(0.25, length(tau4)), repelem(0.5, length(tau3)), ones(1,length(tau_LB))], linearFit2, '--r')
    xlim([0, 1.2])
    ylabel('\tau (min^{-1})')
    xticks([0 1.2/60 0.25 0.5 1])
    xticklabels({'', 'LB, 1.2 s', 'LB, 15 s', 'LB, 30 s', 'LB, 1 min',})
    title('Tau vs Frame Rate')
    txt = ['{\color {blue} slope calculated with 1 min frame rate \color {red} slope calculated without 1 min frame rate}'];
    subtitle(txt)
    pause, close all
end

%without the data from 'LB, frame rate = 1 min', the slope is slightly
%bigger, but they are very close otherwise

%% Correct LB and PBS traces for photobleaching

%load the PBS data
datadir6=dir(['10232021_Exp2' '*dm.mat']); %PBS, incubation = 1 hr, rep1
datadir7=dir(['10262021_Exp2' '*dm.mat']); %PBS, incubation = 1 hr, rep2
datadir8=dir(['11192021_Exp1' '*dm.mat']); %PBS, incubation = 20 min
datadir9=dir(['11192021_Exp2' '*dm.mat']); %PBS, incubation = 2 min

[intensity6, adjintensity6, normintensity6, lCell6, time6, t6]=dataInput(datadir6, 3);
[intensity7, adjintensity7, normintensity7, lCell7, time7, t7]=dataInput(datadir7, 5);
[intensity8, adjintensity8, normintensity8, lCell8, time8, t8]=dataInput(datadir8, 4);
[intensity9, adjintensity9, normintensity9, lCell9, time9, t9]=dataInput(datadir9, 10);

%combine 'PBS, incubation = 1 hr' data
intensity_PBS60 = [intensity6; intensity7(:, 1:51)];
adjintensity_PBS60 = [adjintensity6; adjintensity7(:, 1:51)];
normintensity_PBS60 = [normintensity6; normintensity7(:, 1:49)];

%how does the data look without correcting for photobleaching?
if check3==1
    figure, hold on
    plot(t2, normintensity_LB, 'Color', colorcode{1})
    plot(t6, normintensity_PBS60, 'Color', colorcode{2})
    plot(t8, normintensity8, 'Color', colorcode{3})
    plot(t9, normintensity9, 'Color', colorcode{4})
    title('Normalized Fluorescence vs Time')
    txt = ['{\color[rgb]{0 0.45 0.74}Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
    subtitle(txt)
    xlabel('Time (min)')
    ylabel('Normalized Fluorescence (A.U.)')
end

%correct for photobleaching
[Cnew_LB, dCP_LB, dCB_LB, dCT_LB, unb_frac_LB, Cbl_exp_LB, gamma_LB] = photoCorrect(normintensity_LB, t2, 3);
[Cnew_PBS60, dCP_PBS60, dCB_PBS60, dCT_PBS60, unb_frac_PBS60, Cbl_exp_PBS60, gamma_PBS60] = photoCorrect(normintensity_PBS60, t6, 3);
[Cnew8, dCP8, dCB8, dCT8, unb_frac8, Cbl_exp8, gamma8] = photoCorrect(normintensity8, t8, 3);
[Cnew9, dCP9, dCB9, dCT9, unb_frac9, Cbl_exp9, gamma9] = photoCorrect(normintensity9, t9, 3);

if check4==1
    figure, hold on
    for i=1:height(normintensity_LB)
        plot(t2, Cnew_LB(i,:), 'Color', colorcode{1})
    end
    for i=1:height(normintensity_PBS60)
        plot(t6, Cnew_PBS60(i,:), 'Color', colorcode{2})
    end
    for i=1:height(normintensity8)
        plot(t8, Cnew8(i,:), 'Color', colorcode{3})
    end
    for i=1:height(normintensity9)
        plot(t9, Cnew9(i,:), 'Color', colorcode{4})
    end
    ylabel('Normalized Fluorescence (A.U.)')
    xlabel('Time (min)')
    title('Normalized Fluorescence vs Time')
    txt = ['{\color[rgb]{0 0.45 0.74}Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour  \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
    subtitle(txt)
end

Cavg_LB = mean(Cnew_LB, 1, 'omitnan');
Cavg_PBS60 = mean(Cnew_PBS60, 1, 'omitnan');
Cavg8 = mean(Cnew8, 1, 'omitnan');
Cavg9 = mean(Cnew9, 1, 'omitnan');

if check5==1
    [Cnew_LB, dCP_LB, dCB_LB, dCT_LB, unb_frac_LB, Cbl_exp_LB, gamma_LB] = photoCorrect(normintensity_LB, t2, 1);
    [Cnew_PBS60, dCP_PBS60, dCB_PBS60, dCT_PBS60, unb_frac_PBS60, Cbl_exp_PBS60, gamma_PBS60] = photoCorrect(normintensity_PBS60, t6, 1);
    [Cnew8, dCP8, dCB8, dCT8, unb_frac8, Cbl_exp8, gamma8] = photoCorrect(normintensity8, t8, 1);
    [Cnew9, dCP9, dCB9, dCT9, unb_frac9, Cbl_exp9, gamma9] = photoCorrect(normintensity9, t9, 1);

    figure(1), hold on
    plot(t2, Cavg_LB, 'Color', colorcode{1})
    plot(t6, Cavg_PBS60, 'Color', colorcode{2})
    plot(t8, Cavg8, 'Color', colorcode{3})
    plot(t9, Cavg9, 'Color', colorcode{4})
    title('Normalized Fluorescence vs Time')
    txt = ['{\color[rgb]{0 0.45 0.74}Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
    subtitle(txt)
    xlabel('Time (min)')
    ylabel('Normalized Fluorescence (A.U.)')
    
    [Cnew_LB, dCP_LB, dCB_LB, dCT_LB, unb_frac_LB, Cbl_exp_LB, gamma_LB] = photoCorrect(normintensity_LB, t2, 2);
    [Cnew_PBS60, dCP_PBS60, dCB_PBS60, dCT_PBS60, unb_frac_PBS60, Cbl_exp_PBS60, gamma_PBS60] = photoCorrect(normintensity_PBS60, t6, 2);
    [Cnew8, dCP8, dCB8, dCT8, unb_frac8, Cbl_exp8, gamma8] = photoCorrect(normintensity8, t8, 2);
    [Cnew9, dCP9, dCB9, dCT9, unb_frac9, Cbl_exp9, gamma9] = photoCorrect(normintensity9, t9, 2);
    
    figure(2), hold on
    plot(t2, Cavg_LB, 'Color', colorcode{1})
    plot(t6, Cavg_PBS60, 'Color', colorcode{2})
    plot(t8, Cavg8, 'Color', colorcode{3})
    plot(t9, Cavg9, 'Color', colorcode{4})
    title('Normalized Fluorescence vs Time, parameter 2')
    txt = ['{\color[rgb]{0 0.45 0.74}Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
    subtitle(txt)
    xlabel('Time (min)')
    ylabel('Normalized Fluorescence (A.U.)')
        
    [Cnew_LB, dCP_LB, dCB_LB, dCT_LB, unb_frac_LB, Cbl_exp_LB, gamma_LB] = photoCorrect(normintensity_LB, t2, 3);
    [Cnew_PBS60, dCP_PBS60, dCB_PBS60, dCT_PBS60, unb_frac_PBS60, Cbl_exp_PBS60, gamma_PBS60] = photoCorrect(normintensity_PBS60, t6, 3);
    [Cnew8, dCP8, dCB8, dCT8, unb_frac8, Cbl_exp8, gamma8] = photoCorrect(normintensity8, t8, 3);
    [Cnew9, dCP9, dCB9, dCT9, unb_frac9, Cbl_exp9, gamma9] = photoCorrect(normintensity9, t9, 3);
    
    figure(3), hold on
    plot(t2, Cavg_LB, 'Color', colorcode{1})
    plot(t6, Cavg_PBS60, 'Color', colorcode{2})
    plot(t8, Cavg8, 'Color', colorcode{3})
    plot(t9, Cavg9, 'Color', colorcode{4})
    title('Normalized Fluorescence vs Time, parameter 3')
    txt = ['{\color[rgb]{0 0.45 0.74}Untreated \color[rgb]{0.85 0.33 0.1} PBS, 1 hour \color[rgb]{0.49 0.18 0.56} PBS, 20 min \color[rgb]{0.47 0.67 0.19} PBS, 2 min}'];
    subtitle(txt)
    xlabel('Time (min)')
    ylabel('Normalized Fluorescence (A.U.)')
end
%% Functions
function [intensity, adjintensity, normintensity, lCell, t]=dataInput(datadir, imend)
        
        %pre-allocate variables
        intensity=[];
        lCell=[];
    
        for i=1:length(datadir)
            cd(datadir(i).folder)
            load(datadir(i).name, 'icell_intensity', 'lcell', 'time')

            for n=1:height(icell_intensity)
                if ~isnan(icell_intensity(n, imend))&~isnan(icell_intensity(n, end)) %make sure there are fluor readings during the initial frame and the final frame, otherwise the adjust. and norm. will look off
                    intensity=[intensity; icell_intensity(n, 1:imend)];
                    lCell=[lCell; lcell(n, 1:imend)];
                end
            end

            if i==1
                t=time(1:imend); %new time vector
            end
        end

            %adjust for background
            adjintensity = intensity-intensity(:, end);

            %normalize to initial frame
            normintensity=adjintensity./adjintensity(:,1);
            normintensity(normintensity<0)=0;
            
end


function [fit, tau, yhat]=expModel(time, normintensity, text,  graph)
    fit=[];
    tau=[];
    yhat=[];
    
    for i=1:height(normintensity)

            dni=diff(normintensity(i,:));
            idx=find(dni<0); %this will skip the detergent perfusion reading
            tau_temp = polyfit(time(idx), normintensity(i,idx),1);
            y_hat= polyval(tau_temp, time);  
            
            if tau_temp(1,1)~=0
                fit=[fit, i];
                tau=[tau, abs(tau_temp(1,1))];
                yhat=[yhat; y_hat]; 
            end

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
    beta = 1/((16.9441)*dt + 0.1823); %1/tau
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
        
        if parameter==1            
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
        elseif parameter==2
              for i=1:length(t)-1

                dCB(n,i) = ((normintensity(n,i)+normintensity(n,i+1))/2)/(alpha); %this is the amount of photobleaching that occured in our measured value

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
        elseif parameter==3
               for i=1:length(t)-1

                dCB(n,i) = normintensity(n,i+1)/(alpha); %this is the amount of photobleaching that occured in our measured value

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
        end
    end
    
    gamma = nan(height(normintensity), 1);
    
    for n=1:height(normintensity)
        p = dCP(n,end)/dt;
        gamma(n,1) = p/Cnew(n, end-1);
    end

end